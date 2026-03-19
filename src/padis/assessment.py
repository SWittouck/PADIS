# terms/abbreviations: 
# - ext_region = extended region: DNA sequence around a gene to align 
# - hom_region = homologous region: the section of the extended region that
#   aligns to another extended region
# - term = terminus: first/last N nucleotides of aligned section
# - start: start position of a DNA segment, always smaller than "end" and with a
#   1-base offset (gff style)
# - end: end position of a DNA segment, always greater than "start" and with a
#   1-base offset (gff style)
# - up = upstream
# - down = downstream
# - rc = reverse complement

from Bio import Align
from concurrent.futures import ProcessPoolExecutor
import logging as lg
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta, Sequence
from random import sample
import statistics

from .input import read_acc_genes

#############
# top level #
#############

def assess_orthogroups(
        acc_genes_file: Path, assembly_files: dict[str, Path],
        acc_orthogroups_file: Path, summary_file: Path, max_length: int,
        threads: int, strategy: str
    ) -> None:
    """
    Assess accessory orthogroups for IS status.

    :param acc_genes_file: Path to accessory gene table with columns gene,
        genome, orthogroup, contig, strand, start, end and position.
    :param assembly_files: List of paths to assembly (.fna) files.
    :param acc_orthogroups_file: Path to output file for orthogroup stats.
    :param summary_file: Path to output file for summary stats.
    :max_length: Maximum length of detected insertion sequences.
    :threads: Number of threads to use.
    """

    if acc_orthogroups_file.exists():
        lg.info(
            "Existing accessory orthogroup table file found - skipping phase 2"
        )
        return() 
    
    lg.info("Reading accessory gene table")
    acc_genes = read_acc_genes(acc_genes_file)

    lg.info("Indexing genome sequences (creates .fai files)")
    for p in assembly_files.values():
        _ = Fasta(p)

    # # subset orthogroups for testing purposes
    # topn = acc_genes["orthogroup"].unique()[:20]
    # acc_genes = acc_genes[acc_genes["orthogroup"].isin(topn)]

    if threads == 1: 
        lg.info("Assessing flanking regions per orthogroup")
        acc_ogs = (
            acc_genes
            .groupby("orthogroup")
            .apply(lambda genes: process_orthogroup(
                genes, assembly_files, max_length, strategy
            ))
            .reset_index(drop = True)
        )
    else: 
        lg.info(
            f"Assessing flanking regions per orthogroup using {threads} "
            "threads"
        )
        tasks = [
            (orthogroup, genes, assembly_files, max_length, strategy) for
            orthogroup, genes in acc_genes.groupby("orthogroup", sort = False)
        ]
        with ProcessPoolExecutor(max_workers = threads) as executor: 
            results = list(executor.map(_worker, tasks))
        acc_ogs = pd.DataFrame(results)

    lg.info("Writing candidate insertion sequence orthogroups")
    acc_ogs.to_csv(acc_orthogroups_file, index = False)

    lg.info("Writing summary file")
    (
        acc_ogs["status"]
        .value_counts()
        .reset_index(name = "orthogroups")
        .to_csv(summary_file, index = False)
    )

###############
# lower level #
###############

# helper for running process_orthogroup in parallel
# needs to be top-level function
def _worker(task):
    orthogroup, genes, assembly_files, max_length, strategy = task
    genes.name = orthogroup
    result = process_orthogroup(genes, assembly_files, max_length, strategy)
    return(result)

def process_orthogroup(
        genes: pd.DataFrame, assembly_files: dict[str, Path], max_length: int,
        strategy: str
    ) -> pd.Series:

    orthogroup = genes.name
    lg.debug(f"Processing orthogroup {orthogroup}")

    result = pd.Series({
        "orthogroup": orthogroup,
        "status": "potential IS",
        "genes": genes["gene"].size,
        "genomes": genes["genome"].nunique(),
        "located": genes["position"].count(),
        "positions": genes["position"].nunique(),
        "length": pd.NA,
        "tir_score": pd.NA,
        "tir_nscore": pd.NA,
        "tir_length": pd.NA,
        "tir_offset_up": pd.NA,
        "tir_offset_down": pd.NA,
        "tir_random_score": pd.NA,
        "tir_random_nscore": pd.NA,
        "tir_random_length": pd.NA,
        "tir_up": "",
        "tir_down": ""
    })

    # identify homologous region
    try:
        hom_region = representative_homologous_region(
            genes, assembly_files, max_length, strategy
        )
    except RuntimeError as e:
        result["status"] = e.args[0]
        return(result)

    # determine length of homologous region
    result["length"] = len(hom_region)
    if result["length"] > max_length:
        result["status"] = "too long"
        return(result)

    # set up aligner
    aligner = Align.PairwiseAligner(scoring = "blastn")
    aligner.mode = "local"

    # assess tir of homologous region
    termseq_up = hom_region[:30].seq
    termseq_down = hom_region[-30:].seq
    tir_alignments = aligner.align(termseq_up, termseq_down, strand = "-")
    if not tir_alignments:
        result["tir_score"] = 0
        return(result)
    tir_alignment = tir_alignments[0]
    tirco = tir_alignment.coordinates
    result["tir_score"] = np.int64(tir_alignment.score)
    result["tir_offset_up"] = np.int64(tirco[0, 0])
    result["tir_offset_down"] = np.int64(30 - tirco[1, 0])
    result["tir_length"] = np.int64(tir_alignment.length)
    result["tir_up"] = tir_alignment[0]
    result["tir_down"] = tir_alignment[1]

    # assess randomized tir
    randomseq_down = ''.join(sample(termseq_down,len(termseq_down)))
    tir_alignments = aligner.align(termseq_up, randomseq_down, strand = "-")
    if tir_alignments:
        tir_alignment = tir_alignments[0]
        result["tir_random_score"] = np.int64(tir_alignment.score)
        result["tir_random_length"] = np.int64(tir_alignment.length)

    # normalize tir scores
    scores = []
    for i in range(100):
        randomseq_down = ''.join(sample(termseq_down,len(termseq_down)))
        tir_alignments = aligner.align(termseq_up, randomseq_down, strand = "-")
        if tir_alignments:
            scores.append(tir_alignments[0].score)
        else:
            scores.append(0)
    m = statistics.mean(scores)
    sd = statistics.stdev(scores)
    if sd == 0: return(result)
    result["tir_nscore"] = (result["tir_score"] - m) / sd
    result["tir_random_nscore"] = (result["tir_random_score"] - m) / sd

    return(result)

def representative_homologous_region(
        genes: pd.DataFrame, assembly_files: dict[str, Path], max_length: int,
        strategy: str
    ) -> pd.Series:

    # make copy to avoid modifying original gene table 
    genes = genes.copy()

    # set up aligner
    aligner = Align.PairwiseAligner(scoring = "blastn")
    aligner.mode = "local"

    if strategy == "pairwise_alignment":

        # check if at least two positions
        if genes["position"].nunique() < 2:
            raise RuntimeError("less than two positions determined")

        # identify two example extended regions
        genes = genes[genes["position"].notna()]
        gene1, ext_region1 = representative_extended_region(
            genes, assembly_files, max_length
        )
        genes = genes.loc[(genes["position"] != gene1.position)]
        gene2, ext_region2 = representative_extended_region(
            genes, assembly_files, max_length
        )

        # align extended regions -> identify homologous region
        strand = "+" if gene1.strand == gene2.strand else "-"
        alignment = aligner.align(
            ext_region1.seq, ext_region2.seq, strand = strand
        )[0]
        homco = alignment.coordinates
        hom_region1 = ext_region1[homco[0, 0]:homco[0, -1]]
        f, t = sorted(homco[1, [0, -1]])
        hom_region2 = ext_region2[f:t]

        # check whether gene lies within homologous region
        if (
            hom_region1.start > gene1.start
            or hom_region1.end < gene1.end
            or hom_region2.start > gene2.start
            or hom_region2.end < gene2.end
        ):
            raise RuntimeError("gene outside homologous region")

    elif strategy == "short_contigs":

        if genes["position"].count() == genes["gene"].size:
            raise RuntimeError("all genes located")

        genes = genes[genes["position"].isna()]
        genes["contig_length"] = [
            len(Fasta(assembly_files[gene.genome])[gene.contig]) for gene in
            genes.itertuples()
        ]
        genes = genes[genes["contig_length"] <= max_length]

        if genes.empty:
            raise RuntimeError("too long")

        gene_ix = genes["contig_length"].idxmax()
        gene1 = genes.loc[gene_ix]
        hom_region1 = Fasta(assembly_files[gene1.genome])[gene1.contig][:]

    else:

        lg.error(f"Strategy {strategy} unknown")
        import sys; sys.exit(0)

    if gene1.strand == "-":
        hom_region1 = hom_region1.reverse.complement

    return(hom_region1)

def representative_extended_region(
        genes: pd.DataFrame, assembly_files: dict[str, Path], max_length: int
    ) -> tuple[tuple, Sequence]:
    """
    Find the gene with the best surrounding region.

    The best region is the longest (within 2 * max_length - gene_length),
    excluding uncalled bases.

    :param assembly_files: Assembly path (fna file) for every genome
    :param genes: Coordinates of all genes 
    :param max_length: Maximum length of insertion sequence
    """

    best_gene = None
    best_region = None
    best_length = -1
    perfect_region = True

    for gene in genes.itertuples():
        region_start = min(gene.end - max_length, gene.start)
        region_end = max(gene.start + max_length, gene.end)
        if region_start < 1:
            perfect_region = False
            region_start = 1
        assembly = Fasta(assembly_files[gene.genome])
        # important: pyfaidx sequence slicing uses a 0-base offset but gff 
        # start/end positions use a 1-base offset!
        region = assembly[gene.contig][(region_start - 1):region_end]
        if region.end < region_end:
            perfect_region = False
        uncalled_bases = region.seq.count("N")
        # punish length for uncalled bases
        region_length = len(region.seq) - uncalled_bases
        if uncalled_bases > 0:
            perfect_region = False
        if region_length > best_length:
            best_length = region_length
            best_gene = gene
            best_region = region
        if perfect_region:
            return(gene, region)

    return(best_gene, best_region)
