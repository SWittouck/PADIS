# terms/abbreviations: 
# - ext_region = extended region: DNA region around a gene
# - hom_region = homologous region: the section of an extended region that
#   aligns to another extended region
# - term = terminus: first/last N nucleotides of aligned section
# - start: start position of a DNA segment, with a 1-base offset (gff style)
# - end: end position of a DNA segment, with a 1-base offset (gff style)
# - up = upstream
# - down = downstream
# - rc = reverse complement

# remarks:
# - The objects ext_region and hom_region are pyfaidx Sequence objects. They
#   have a start and end coordinate that indicate their position on their
#   contig. A Sequence object can be reverse complemented. In that case, the
#   start position will be larger than the end position.

from Bio import Align
from concurrent.futures import ProcessPoolExecutor
import logging as lg
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta, Sequence
from random import sample
import statistics
import subprocess
import tempfile

from .input import read_acc_genes

#############
# top level #
#############

def assess_orthogroups(
        acc_genes_file: Path, assembly_files: dict[str, Path],
        acc_orthogroups_file: Path, summary_file: Path, max_length: int,
        threads: int, strategy: str, ali_dir: Path = None
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
    
    if strategy == "pairwise_alignment":
        if not ali_dir:
            lg.error(
                "Pairwise alignment strategy needs a folder for alignments"
            )
            sys.exit(1)
        lg.info("Creating output folder for alignments")
        ali_dir.mkdir(exist_ok = True)

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
                genes, assembly_files, max_length, strategy, ali_dir
            ))
            .reset_index(drop = True)
        )
    else: 
        lg.info(
            f"Assessing flanking regions per orthogroup using {threads} "
            "threads"
        )
        tasks = [
            (
                orthogroup, genes, assembly_files, max_length, strategy, ali_dir
            ) for
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
    orthogroup, genes, assembly_files, max_length, strategy, ali_dir = task
    genes.name = orthogroup
    result = process_orthogroup(
        genes, assembly_files, max_length, strategy, ali_dir
    )
    return(result)

def process_orthogroup(
        genes: pd.DataFrame, assembly_files: dict[str, Path], max_length: int,
        strategy: str, ali_dir: Path
    ) -> pd.Series:

    term_length = 100

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
        "fdr_score": pd.NA,
        "fdr_length": pd.NA,
        "fdr_offset_up": pd.NA,
        "fdr_offset_down": pd.NA,
        "tir_up": "",
        "tir_down": "",
        "fdr_up": "",
        "fdr_down": ""
    })

    if strategy in ["pairwise_alignment", "multiple_alignment"]:
        if result["genes"] == 1:
            result["status"] = "singleton orthogroup"
            return(result)

    # identify homologous region
    try:
        hom_region = representative_homologous_region(
            genes, assembly_files, max_length, strategy, ali_dir
        )
    except RuntimeError as e:
        result["status"] = e.args[0]
        return(result)

    # determine length of homologous region
    result["length"] = len(hom_region)
    if result["length"] > max_length:
        result["status"] = "homologous region too long"
        return(result)

    # set up aligner
    aligner = Align.PairwiseAligner(scoring = "blastn")
    aligner.mode = "local"

    # assess tir of homologous region
    term_length = min(len(hom_region) // 2, term_length)
    termseq_up = hom_region[:term_length].seq
    termseq_down = hom_region[-term_length:].reverse.complement.seq
    aliresult = pairwise_alignment(aligner, termseq_up, termseq_down)
    result["tir_score"] = aliresult[0]
    result["tir_length"] = aliresult[1]
    result["tir_offset_up"] = aliresult[2]
    result["tir_offset_down"] = aliresult[4]
    result["tir_up"] = aliresult[6]
    result["tir_down"] = aliresult[7]

    # assess fdr of homologous region
    termseq_down = hom_region[-term_length:].seq
    aliresult = pairwise_alignment(aligner, termseq_up, termseq_down)
    result["fdr_score"] = aliresult[0]
    result["fdr_length"] = aliresult[1]
    result["fdr_offset_up"] = aliresult[2]
    result["fdr_offset_down"] = term_length - aliresult[5]
    result["fdr_up"] = aliresult[6]
    result["fdr_down"] = aliresult[7]

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
        strategy: str, ali_dir: Path
    ) -> pd.Series:

    orthogroup = genes.name

    # make copy to avoid modifying original gene table 
    genes = genes.copy()

    if strategy == "multiple_alignment":

        extended_regions = [
            representative_extended_region(
                genes, assembly_files, max_length
            )[1] for position, genes in genes.groupby("position")
        ]

        ali = multiple_alignment(extended_regions[0:10])

        ali_file = Path("~/padis_alignments") / f"{orthogroup}.aln"
        with open(ali_file, "w") as ali_handle:
            ali_handle.write(ali)

        hom_region1 = extended_regions[0]

        # to finish

    elif strategy == "pairwise_alignment":

        # identify two example extended regions
        genome, copy_number = (
            genes["genome"]
            .value_counts()
            .agg(["idxmax", "max"])
        )
        if genes["position"].nunique() > 1:
            genes = genes[genes["position"].notna()]
            gene1, ext_region1 = representative_extended_region(
                genes, assembly_files, max_length
            )
            genes = genes.loc[genes["position"] != gene1.position]
            gene2, ext_region2 = representative_extended_region(
                genes, assembly_files, max_length
            )
        elif copy_number > 1:
            genes = genes[genes["genome"] == genome]
            gene1, ext_region1 = representative_extended_region(
                genes, assembly_files, max_length
            )
            genes = genes.loc[genes["gene"] != gene1.gene]
            gene2, ext_region2 = representative_extended_region(
                genes, assembly_files, max_length
            )
        else:
            raise RuntimeError("no multipositionality evidence")

        # set up aligner
        aligner = Align.PairwiseAligner(scoring = "blastn", mode = "local")

        # align extended regions -> identify homologous region
        alignment_up = aligner.align(
            ext_region1.seq[:max_length], ext_region2.seq[:max_length],
            strand = "+"
        )[0]
        homco_up = alignment_up.coordinates
        alignment_down = aligner.align(
            ext_region1.seq[-max_length:], ext_region2.seq[-max_length:],
            strand = "+"
        )[0]
        homco_down = alignment_down.coordinates
        gene1_length = gene1.end - gene1.start + 1
        start = homco_up[0, 0]
        end = homco_down[0, -1] + max_length - gene1_length
        hom_region1 = ext_region1[start:end]

        # write alignments
        Align.write(alignment_up, ali_dir / f"{orthogroup}_up.aln", "fasta")
        Align.write(alignment_down, ali_dir / f"{orthogroup}_down.aln", "fasta")

        # hits_up = blastn_pairwise(
        #     ext_region1.seq[:max_length], ext_region2.seq[:max_length]
        # )
        # if hits_up.empty:
        #     raise RuntimeError("no homology detected")
        # start = hits_up["qstart"].max() - 1
        # hits_down = blastn_pairwise(
        #     ext_region1.seq[-max_length:], ext_region2.seq[-max_length:]
        # )
        # if hits_down.empty:
        #     raise RuntimeError("no homology detected")
        # gene1_length = gene1.end - gene1.start + 1
        # end = hits_down["qend"].min() + max_length - gene1_length
        # if start >= end:
        #     raise RuntimeError("gene not fully homologous")
        # hom_region1 = ext_region1[start:end]

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
            raise RuntimeError("homologous region too long")

        gene_ix = genes["contig_length"].idxmax()
        gene1 = genes.loc[gene_ix]
        hom_region1 = Fasta(assembly_files[gene1.genome])[gene1.contig][:]

    else:

        lg.error(f"Strategy {strategy} unknown")
        import sys; sys.exit(0)

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
        # example: a region of 100 nucleotides starting at position 500 in gff
        # coordinates would be 500 to 599 --> 500 to (500 + 100 - 1)
        region_start = min(gene.end - max_length + 1, gene.start)
        region_end = max(gene.start + max_length - 1, gene.end)
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
            break

    if best_gene.strand == "-":
        best_region = best_region.reverse.complement

    return(best_gene, best_region)

def multiple_alignment(seqs):
    """
    Align a list of sequences with MAFFT.
    """

    fasta = (
        "\n"
        .join(f">seq{i}\n{seq}" for i, seq in enumerate(seqs))
        .encode()
    )

    ali = subprocess.run(
        ["mafft", "-"], input = fasta, stdout = subprocess.PIPE,
        stderr = subprocess.PIPE
    )

    return(ali.stdout.decode())

def pairwise_alignment(aligner, seq1, seq2):
    alignments = aligner.align(seq1, seq2)
    if not alignments:
        return((0, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA, pd.NA))
    alignment = alignments[0]
    alico = alignment.coordinates
    score = np.int64(alignment.score)
    length = np.int64(alignment.length)
    ali1_start = np.int64(alico[0, 0])
    ali1_end = np.int64(alico[0, -1])
    ali2_start = np.int64(alico[1, 0])
    ali2_end = np.int64(alico[1, -1])
    ali1 = alignment[0]
    ali2 = alignment[1]
    return((
        score, length, ali1_start, ali1_end, ali2_start, ali2_end, ali1, ali2
    ))

def blastn_pairwise(query_seq, subject_seq):
    """
    Align two sequences with blastn.
    """

    with (
        tempfile.NamedTemporaryFile(
            mode = "w", suffix = ".fa", delete = True, delete_on_close = False
        ) as query_handle,
        tempfile.NamedTemporaryFile(
            mode = "w", suffix = ".fa", delete = True, delete_on_close = False
        ) as subject_handle
    ):

        query_handle.write(">query\n" + query_seq + "\n")
        subject_handle.write(">subject\n" + subject_seq + "\n")

        query_file = query_handle.name
        subject_file = subject_handle.name

        query_handle.close()
        subject_handle.close()

        cmd = [
            "blastn",
            "-query", query_file,
            "-subject", subject_file,
            "-outfmt", "6",
            "-strand", "plus"
        ]

        result = subprocess.run(
            cmd, capture_output = True, text = True, check = True
        )

        output = result.stdout.strip()

        colnames = [
            "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        ]
        hits = [line.split("\t") for line in output.splitlines()]
        hits = pd.DataFrame(hits, columns = colnames)

        numeric_cols = [
            "pident", "length", "mismatch", "gapopen", "qstart", "qend",
            "sstart", "send", "evalue", "bitscore"
        ]
        hits[numeric_cols] = hits[numeric_cols].apply(pd.to_numeric)

    return(hits)
