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
        threads: int
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
                genes, assembly_files, max_length
            ))
            .reset_index(drop = True)
        )
    else: 
        lg.info(
            f"Assessing flanking regions per orthogroup using {threads} "
            "threads"
        )
        tasks = [(orthogroup, genes, assembly_files, max_length) for
            orthogroup, genes in acc_genes.groupby("orthogroup", sort = False)]
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

def process_orthogroup(
        genes: pd.DataFrame, assembly_files: dict[str, Path],
        max_length: int
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

    if result["genes"] == 1:
        result["status"] = "singleton"
        return(result)
    if result["located"] <= 1:
        result["status"] = "insufficient positions"
        return(result)
    if result["positions"] == 1:
        result["status"] = "no position variation"
        return(result)

    # make copy to avoid modifying original gene table 
    genes = genes.copy()

    # identify two example extended regions
    genes = genes[genes["position"].notna()]
    gene1, ext_region1 = best_region(assembly_files, genes, max_length)
    genes = genes.loc[(genes["position"] != gene1.position)]
    gene2, ext_region2 = best_region(assembly_files, genes, max_length)

    # align extended regions -> identify homologous region
    aligner = Align.PairwiseAligner(scoring = "blastn")
    aligner.mode = "local"
    strand = "+" if gene1.strand == gene2.strand else "-"
    alignment = aligner.align(
        ext_region1.seq, ext_region2.seq, strand = strand
    )[0]
    homco = alignment.coordinates
    hom_region1 = ext_region1[homco[0, 0]:homco[0, -1]]
    f, t = sorted(homco[1, [0, -1]])
    hom_region2 = ext_region2[f:t]

    # determine length of homologous region
    result["length"] = len(hom_region1)
    if result["length"] > max_length:
        result["status"] = "too long"
        return(result)

    # check whether gene lies within homologous region
    if (
        hom_region1.start > gene1.start
        or hom_region1.end < gene1.end
        or hom_region2.start > gene2.start
        or hom_region2.end < gene2.end
    ):
        result["status"] = "gene outside homologous region"
        return(result)

    # assess tir of homologous region
    if gene1.strand == "-":
        ext_region1 = ext_region1.reverse.complement
    termseq_up = hom_region1[:30].seq
    termseq_down = hom_region1[-30:].seq
    tir_alignment = aligner.align(termseq_up, termseq_down, strand = "-")[0]
    tir_co = tir_alignment.coordinates
    result["tir_score"] = np.int64(tir_alignment.score)
    result["tir_offset_up"] = np.int64(tir_co[0, 0])
    result["tir_offset_down"] = np.int64(30 - tir_co[1, 0])
    result["tir_length"] = np.int64(tir_alignment.length)
    result["tir_up"] = tir_alignment[0]
    result["tir_down"] = tir_alignment[1]

    # assess randomized tir
    randomseq_down = ''.join(sample(termseq_down,len(termseq_down)))
    tir_alignment = aligner.align(termseq_up, randomseq_down, strand = "-")[0]
    result["tir_random_score"] = np.int64(tir_alignment.score)
    result["tir_random_length"] = np.int64(tir_alignment.length)

    # normalize tir scores
    scores = []
    for i in range(100):
        randomseq_down = ''.join(sample(termseq_down,len(termseq_down)))
        tir_alignment = aligner.align(
            termseq_up, randomseq_down, strand = "-"
        )[0]
        scores.append(tir_alignment.score)
    m = statistics.mean(scores)
    sd = statistics.stdev(scores)
    result["tir_nscore"] = (result["tir_score"] - m) / sd
    result["tir_random_nscore"] = (result["tir_random_score"] - m) / sd

    return(result)

# helper for running process_orthogroup in parallel 
# needs to be top-level function
def _worker(task):
    orthogroup, genes, assembly_files, max_length = task
    genes.name = orthogroup
    result = process_orthogroup(genes, assembly_files, max_length)
    return(result)

def best_region(
        assembly_files: dict[str, Path], genes: pd.DataFrame, max_length: int
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
