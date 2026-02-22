from Bio import Align
import logging as lg
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta, Sequence

from .input import read_canisgenes

#############
# top level #
#############

def assess_canisorthogroups(
        canisgenes_file: Path, assemblies_files: dict[str, Path], 
        canisorthogroups_file: Path) -> None:
    
    flank_size = 3000

    if canisorthogroups_file.exists():
        lg.info("Existing canisorthogroups file found - skipping phase 2")
        return() 
    
    lg.info("Reading candidate insertion sequence genes")
    canisgenes = read_canisgenes(canisgenes_file)
    genomes = canisgenes["genome"].unique().tolist()

    lg.info("Indexing genome sequences")
    genomes = {p.stem: Fasta(p) for p in assemblies_files.values()}

    # # for testing purposes
    # topn = canisgenes["orthogroup"].unique()[:20]
    # canisgenes = canisgenes[canisgenes["orthogroup"].isin(topn)]

    lg.info("Calculating basic orthogroup stats")
    canisogs1 = (
        canisgenes
        .groupby("orthogroup")
        .agg(
            n_genes = ("gene", "size"),
            genomes = ("genome", "nunique"),
            positions = ("position", "nunique")
        )
    )

    lg.info("Aligning terminal regions where possible for every orthogroup")
    canisogs2 = (
        canisgenes
        .groupby("orthogroup")
        .apply(lambda g: process_orthogroup(g, genomes, flank_size))
    )

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

def process_orthogroup(genes, genomes, flank_size):
    
    FS = flank_size

    lg.debug(f"Processing orthogroup {genes.name}")

    genes = genes[genes["position"].notna()]

    gene1, seq1 = sequence_with_flanks(genomes, genes, FS)
    genes = genes.loc[(genes["position"] != gene1.position)]
    gene2, seq2 = sequence_with_flanks(genomes, genes, FS)

    if not seq1 or not seq2:
        return(pd.Series({
            "tir_score": 0,
            "fdr_score": 0,
            "comment": "All flanking regions too short"
        }))

    aligner = Align.PairwiseAligner(scoring = "blastn")
    aligner.mode = "local"
    alignments = aligner.align(seq1.seq, seq2.seq, strand = "+")
    alignment = alignments[0]
    alicoo = alignment.coordinates

    if (
        alicoo[0, 0] > FS
        or (alicoo[0, -1] < len(seq1) - FS)
        or alicoo[1, 0] > FS
        or alicoo[1, -1] < len(seq2) - FS
    ):
        return(pd.Series({
            "tir_score": 0,
            "fdr_score": 0,
            "comment": "Genes not fully inside the alignment"
        }))

    if (
        alicoo[0, 0] == 0 
        or alicoo[1, 0] == 0
        or alicoo[0, 1] == FS + len(seq1)
        or alicoo[1, 1] ==  FS + len(seq2)
    ): 
        return(pd.Series({
            "tir_score": 0,
            "fdr_score": 0,
            "comment": "Alignment extends beyond flanks"
        }))

    term_left = seq1[alicoo[0, 0]:alicoo[0, -1]][:30]
    term_right = seq1[alicoo[0, 0]:alicoo[0, -1]][-30:]
    term_right_rc = term_right.reverse.complement
    tir_alignments = aligner.align(term_left.seq, term_right_rc.seq)
    fdr_alignments = aligner.align(term_left.seq, term_right.seq)
    return(pd.Series({
        "tir_score": tir_alignments[0].score,
        "fdr_score": fdr_alignments[0].score,
        "comment": "Termini extraction successful"
    }))

def sequence_with_flanks(
        genomes: dict[Fasta], genes: pd.DataFrame, flanksize: int) -> (tuple, 
        Sequence):
    """
    Sample a gene and extract its sequencing including flanking regions.

    If no complete flanks can be extracted or if many uncalled bases are found,
    a different gene will be sampled. If no genes meet these criteria, None will
    be returned as sequence. 
    
    :param genomes: Genome sequence (fna file) for every genome
    :param genes: Coordinates of all genes 
    :param flanksize: Size of the flanking region to extract
    """

    gene = None
    seq = None

    for gene in genes.itertuples():
        start = gene.start - flanksize
        end = gene.end + flanksize
        if start < 0: continue
        seq = genomes[gene.genome][gene.contig][start:end]
        length = gene.end - gene.start
        if len(seq.seq) < length + 2 * flanksize: continue
        if seq.seq.count("N") > 10: continue
        if gene.strand == "-": seq = seq.reverse.complement
        break 

    return(gene, seq)
