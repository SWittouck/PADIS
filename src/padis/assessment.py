# terms/abbreviations: 
# - reg = region: DNA sequence around a gene to align 
# - ali = aligned section: part of a region that aligns to another region
# - term = terminus: first/last N nucleotides of aligned section
# - start: start position of a DNA segment, always smaller than end
# - end: end position of a DNA segment, always greater than start
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

from .input import read_canisgenes

#############
# top level #
#############

def assess_canisorthogroups(
        canisgenes_file: Path, assemblies_files: dict[str, Path],
        canisorthogroups_file: Path, max_length: int, threads: int
        ) -> None:

    if canisorthogroups_file.exists():
        lg.info("Existing canisorthogroups file found - skipping phase 2")
        return() 
    
    lg.info("Reading candidate insertion sequence genes")
    canisgenes = read_canisgenes(canisgenes_file)

    lg.info("Indexing genome sequences (creates .fai files)")
    for p in assemblies_files.values():
        _ = Fasta(p)

    # # subset orthogroups for testing purposes
    # topn = canisgenes["orthogroup"].unique()[:20]
    # canisgenes = canisgenes[canisgenes["orthogroup"].isin(topn)]

    lg.info("Calculating basic orthogroup stats")
    canisogs1 = (
        canisgenes
        .groupby("orthogroup")
        .agg(
            genes = ("gene", "size"),
            genomes = ("genome", "nunique"),
            positions = ("position", "nunique")
        )
    )

    if threads == 1: 
        lg.info("Assessing flanking regions per orthogroup")
        canisogs2 = (
            canisgenes
            .groupby("orthogroup")
            .apply(
                lambda g: process_orthogroup(g, assemblies_files, max_length)
            )
        )
    else: 
        lg.info(
            f"Assessing flanking regions per orthogroup using {threads} "
            "threads")
        tasks = [(orthogroup, genes, assemblies_files, max_length) for
            orthogroup, genes in canisgenes.groupby("orthogroup", sort = False)]
        with ProcessPoolExecutor(max_workers = threads) as executor: 
            results = list(executor.map(_worker, tasks))
        canisogs2 = pd.DataFrame(results).set_index("orthogroup")

    lg.info("Merging orthogroup statistics")
    canisogs = (
        canisogs1
        .merge(canisogs2, on = "orthogroup")
        .reset_index()
    )

    lg.info("Writing candidate insertion sequence orthogroups")
    canisogs.to_csv(canisorthogroups_file, index = False)

###############
# lower level #
###############

def process_orthogroup(
        genes: pd.DataFrame, assemblies_files: dict[str, Path], 
        max_length: int
        ) -> pd.Series:

    orthogroup = genes.name
    lg.debug(f"Processing orthogroup {orthogroup}")

    result = pd.Series({
        "length": 0,
        "tir_score": 0,
        "tir_offset_up": 0,
        "tir_offset_down": 0,
        "tir_length": 0,
        "tir_up": "",
        "tir_down": "",
        "fdr_score": 0,
        "fdr_offset_up": 0,
        "fdr_offset_down": 0,
        "fdr_length": 0,
        "fdr_up": "",
        "fdr_down": "",
        "too_short": False,
        "too_long": False,
        "includes_region_boundary": False
    })

    # make copy to avoid modifying original gene table 
    genes = genes.copy()

    genes = genes[genes["position"].notna()]

    # open connections to assembly files
    # remark: indexing will not happen again since .fai files already exist
    genomes = {p.stem: Fasta(p) for p in assemblies_files.values()}

    gene1, reg1 = best_region(genomes, genes, max_length)
    genes = genes.loc[(genes["position"] != gene1.position)]
    gene2, reg2 = best_region(genomes, genes, max_length)

    # remark: this should normally not happen
    if not reg1 or not reg2:
        lg.warning(f"Unable to find two regions to align for {orthogroup}")
        return(result)

    aligner = Align.PairwiseAligner(scoring = "blastn")
    aligner.mode = "local"
    strand = "+" if gene1.strand == gene2.strand else "-"
    alignments = aligner.align(reg1.seq, reg2.seq, strand = strand)
    alignment = alignments[0]
    alico = alignment.coordinates
    ali1_start, ali1_end = [reg1.start + c for c in sorted(alico[0, [0, -1]])]
    ali2_start, ali2_end = [reg2.start + c for c in sorted(alico[1, [0, -1]])]
    result["length"] = ali1_end - ali1_start

    if (
        ali1_start > gene1.start
        or ali1_end < gene1.end
        or ali2_start > gene2.start
        or ali2_end < gene2.end
    ):
        result["too_short"] = True
        return(result)

    if (ali1_end - ali1_start) > max_length:
        result["too_long"] = True

    if (
        ali1_start == reg1.start
        or ali1_end == reg1.end
        or ali2_start == reg2.start
        or ali2_end == reg2.end
    ): 
        result["includes_region_boundary"] = True

    # switch to coordinates relative to region
    # alico[0, 0] is always smaller than alico[0, -1]
    # (but if gene.strand == "-", alico[1, 0] is larger than alico[1, -1])
    if gene1.strand == "+":
        ali_start = alico[0, 0]
        ali_end = alico[0, -1]
    else:
        reg1 = reg1.reverse.complement
        ali_start = len(reg1.seq) - alico[0, -1]
        ali_end = len(reg1.seq) - alico[0, 0]

    # assess tir
    term_up = reg1[ali_start:ali_end][:30]
    term_down = reg1[ali_start:ali_end][-30:]
    tir_alignments = aligner.align(term_up.seq, term_down.seq, strand = "-")
    tir_co = tir_alignments[0].coordinates
    result["tir_score"] = np.int64(tir_alignments[0].score)
    result["tir_offset_up"] = np.int64(tir_co[0, 0])
    result["tir_offset_down"] = np.int64(tir_co[1, 0] - 30)
    result["tir_length"] = np.int64(tir_co[0, -1] - tir_co[0, 0])
    result["tir_up"] = term_up[tir_co[0, 0]:tir_co[0, -1]]
    result["tir_down"] = term_down[tir_co[1, -1]:tir_co[1, 0]]
    result["tir_down"] = result["tir_down"].reverse.complement

    # assess fdr
    flank_up = reg1[(ali_start - 10):ali_start]
    flank_down = reg1[ali_end:(ali_end + 10)]
    try:
        fdr_alignments = aligner.align(flank_up.seq, flank_down.seq)
    except ValueError: 
        return(result)
    try:
        fdr_co = fdr_alignments[0].coordinates
    except IndexError: 
        return(result)
    result["fdr_score"] = np.int64(fdr_alignments[0].score)
    result["fdr_offset_up"] = np.int64(fdr_co[0, -1] - 10 - 1)
    result["fdr_offset_down"] = np.int64(fdr_co[1, 0] + 1)
    result["fdr_length"] = np.int64(fdr_co[0, -1] - fdr_co[0, 0])
    result["fdr_up"] = flank_up[fdr_co[0, 0]:fdr_co[0, -1]]
    result["fdr_down"] = flank_down[fdr_co[1, 0]:fdr_co[1, -1]]

    return(result)

# helper for running process_orthogroup in parallel 
# needs to be top-level function
def _worker(task):
    orthogroup, genes, assemblies_files, max_length = task
    genes.name = orthogroup
    result = process_orthogroup(genes, assemblies_files, max_length)
    result.loc["orthogroup"] = orthogroup
    return(result)

def best_region(
        genomes: dict[Fasta], genes: pd.DataFrame, max_length: int
        ) -> tuple[tuple, Sequence]:
    """
    Find the gene with the best surrounding region.

    The best region is the longest (within 2 * max_length - gene_length),
    excluding uncalled bases.

    :param genomes: Genome sequence (fna file) for every genome
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
        if region_start < 0:
            perfect_region = False
            region_start = 0
        region = genomes[gene.genome][gene.contig][region_start:region_end]
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
