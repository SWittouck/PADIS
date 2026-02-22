# abbrevations:
# canis: candidate insertion sequence
# posind: position indicator, e.g. "orthogroup1+" = upstream of orthogroup 1

import logging as lg
import pandas as pd
from pathlib import Path

from .input import read_annotation

#############
# top level #
#############

def identify_canisorthogroups(
        annotation_files: dict[str, Path], genes: pd.DataFrame,
        canisgenes_file: Path, intervals_file: Path = None) -> None:
    """
    Identify candidate insertion sequence genes. 
    
    :param annotation_files: Annotation files (gff). 
    :param genes: Pangenome in SCARAP format (columns gene, genome and 
        orthogroup)
    :param canisgenes_file: Path to output file for candidate insertion sequence 
        genes. 
    :param intervals_file: Path to output file for intervals (left and right
        position indicators). 
    """

    if canisgenes_file.exists():
        lg.info("Existing canisgenes file found - skipping phase 1")
        return() 
    
    lg.info("Processing pangenome")
    genomes = genes["genome"].unique()
    orthogroups = genes["orthogroup"].unique()
    lg.info(
        f"Pangenome contains {len(genes)} genes, {len(genomes)} genomes "
        f"and {len(orthogroups)} orthogroups")

    lg.info("Identifying single-copy core orthogroups")
    core, multicopy, zerocopy = determine_core(genes) 
    lg.info(f"Identified {len(core)} 95% single-copy core orthogroups")
    lg.info(
        f"Identified {len(multicopy)} instances of a core orthogroup copy "
        "number greater than one")
    lg.info(
        f"Identified {len(zerocopy)} instances of a core orthogroup copy "
        "number of zero")
    
    iscore = genes["orthogroup"].isin(core)
    coregenes, accgenes = genes[iscore], genes[~ iscore]
    coregene2orthogroup = dict(zip(coregenes["gene"], coregenes["orthogroup"]))

    lg.info("Processing gene annotation files")
    accgenes_chunks = []
    intervals_chunks = []
    for genome in genomes: 
        lg.debug(f"Processing genome {genome}")
        annotation = read_annotation(annotation_files[genome])
        # intervals needs to be separate from accgenes because we also want to
        # keep track of empty intervals
        accgenes_chunk, intervals_chunk = process_annotation(
            annotation, coregene2orthogroup)
        intervals_chunk["genome"] = genome
        accgenes_chunks.append(accgenes_chunk)
        intervals_chunks.append(intervals_chunk)
    accgenes2 = pd.concat(accgenes_chunks).set_index("gene")
    accgenes = accgenes.set_index("gene").join(accgenes2) 
    accgenes = accgenes.reset_index() 
    intervals = pd.concat(intervals_chunks)
    intervals = intervals[["genome", "interval", "posind1", "posind2"]]

    lg.info("Removing multi-copy core gene instances from intervals")
    intervals["tocheck"] = list(zip(
        intervals["genome"], intervals["posind1"].str[:-1]))
    intervals.loc[intervals["tocheck"].isin(multicopy), "posind1"] = None
    intervals["tocheck"] = list(zip(
        intervals["genome"], intervals["posind2"].str[:-1]))
    intervals.loc[intervals["tocheck"].isin(multicopy), "posind2"] = None
    intervals = intervals.drop(columns = ["tocheck"])

    lg.info("Removing singleton position indicator pairs")
    intervals["is_singleton"] = ~ intervals.duplicated(
        subset = ["posind1", "posind2"], keep = False)
    intervals.loc[intervals["is_singleton"], "posind1"] = None
    intervals.loc[intervals["is_singleton"], "posind2"] = None
    intervals = intervals.drop(columns = ["is_singleton"])

    lg.info(
        "Removing position indicator pairs likely created due to core gene "
        "absence")
    ind2pos = define_positions(intervals) # ind2pos = posind2position
    intervals["position"] = intervals["posind1"].map(ind2pos).astype("Int64")
    # we need to convert the zerocopy instances from (genome, orthogroup) to
    # (genome, position), where every orthogroup corresponds to two positions:
    # one upstream and one downstream
    zerocopy_positions = []
    for genome, orthogroup in zerocopy:
        try:
            zerocopy_positions.append((genome, ind2pos[orthogroup + "+"]))
            zerocopy_positions.append((genome, ind2pos[orthogroup + "-"]))
        except KeyError as e:
            lg.warning(
                f"Zerocopy position indicator {e.args[0]} not found in "
                "position indicator table")
    zerocopy_positions = set(zerocopy_positions)
    # every interval has a position identifier; we have to check those against
    # the zerocopy positions
    intervals["tocheck"] = list(zip(intervals["genome"], intervals["position"]))
    intervals["may_miss_core"] = intervals["tocheck"].isin(zerocopy_positions)
    intervals.loc[intervals["may_miss_core"], "posind1"] = None
    intervals.loc[intervals["may_miss_core"], "posind2"] = None
    intervals = intervals.drop(
        columns = ["tocheck", "may_miss_core", "position"])

    if intervals_file: 
        lg.info("Writing interval table")
        intervals.to_csv(intervals_file, index = False)

    lg.info("Defining genomic positions")
    ind2pos = define_positions(intervals)

    lg.info("Assigning genomic positions to accessory genes")
    intervals["position"] = intervals["posind1"].map(ind2pos).astype("Int64")
    intervals = intervals.set_index(["genome", "interval"])
    accgenes = accgenes.set_index(["genome", "interval"])
    accgenes = accgenes.join(intervals)
    accgenes = accgenes.reset_index().drop(columns = ["interval"])
    lg.info(
        f"Assigned {accgenes['position'].nunique()} unique genomic positions")
    n_withpos = accgenes['position'].count()
    n_acc = len(accgenes)
    lg.info(f"Assigned a position to {n_withpos} out of {n_acc} accessory "
        f"genes ({n_withpos / n_acc * 100:.1f}%)")

    lg.info("Selecting candidate insertion sequence genes")
    canisgenes = (
        accgenes
        .groupby("orthogroup")
        .filter(lambda g: g["position"].nunique() > 1)
    )
    lg.info(
        f"Selected {len(canisgenes)} candidate insertion sequence genes "
        f"belonging to {canisgenes['orthogroup'].nunique()} orthogroups"
    )

    lg.info("Writing candidate insertion sequence genes")
    canisgenes = canisgenes[[
        "gene", "genome", "orthogroup", "contig", "strand", "start", "end", 
        "position"]]
    canisgenes.to_csv(canisgenes_file, index = False)

###############
# lower level #
###############

def determine_core(genes: pd.DataFrame) -> (set[str], set[tuple[str, str]], 
        set[tuple[str, str]]):
    """
    Identify the 95% single-copy core orthogroups in a pangenome. 
    
    :param genes: Table with columns gene, genome and orthogroup. 
    :return: Single-copy core orthogroups. 
    :return: Multi-copy occurrences of core orthogroups in genomes. 
    :return: Zero-copy occurrences of core orthogroups in genomes. 
    """
    counts = (
        genes
        .groupby(["genome", "orthogroup"])
        .size()
        .unstack(fill_value = 0)
    )
    n_genomes = counts.shape[0]
    singlecopy_prevalence = (counts == 1).sum(axis = 0) / n_genomes 
    core = singlecopy_prevalence[singlecopy_prevalence >= 0.95].index.tolist() 
    counts_core = counts[core]
    multicopy = counts_core.stack().loc[lambda x: x >= 2].index.tolist()
    zerocopy = counts_core.stack().loc[lambda x: x == 0].index.tolist()
    core = set(core)
    multicopy = set(multicopy)
    zerocopy = set(zerocopy)
    return(core, multicopy, zerocopy)

def process_annotation(
        annotation: pd.DataFrame, coregene2orthogroup: dict[str, str]) -> (
        pd.DataFrame, pd.DataFrame):
    """
    Process gene annotation of a genome (helper of identify_canisgenes). 

    Implementation note: an interval has the form (interval_id, posind1,
    posind2) where posind1 and posind2 are sorted in lexicographical order. 
    
    :param annotation: Table with columns of gff file. 
    :param coregene2orthogroup: Dictionary with orthogroups of core genes.
    :return: Table with accessory genes; includes interval column.
    :return: Table with intervals; has columns interval, posind1, posind2.
    """

    # reconstruct gene id that prodigal uses in ffn/faa files but not in gff 
    # files
    annotation["gene"] = (
        annotation["seqid"] 
        + "_" 
        + annotation["attr"].str.extract(r"^[^_]+_([^;]+)")[0]
    )
    annotation = annotation.set_index("gene", drop = False)

    annotation["interval"] = pd.Series(dtype = "Int64")
    intervals = []
    
    contig = annotation["seqid"].tolist()[0]
    interval = 0
    members = []
    left_posind = None
    right_posind = None
    
    # itertuples is faster than iterrows (tested)
    for row in annotation.itertuples():
        
        # if new contig: resolve interval 
        if row.seqid != contig:

            annotation.loc[members, "interval"] = interval
            intervals.append((interval, left_posind, None))
            
            contig = row.seqid
            interval += 1
            members = []
            left_posind = None
            right_posind = None

        # if core gene: create position indicators and resolve interval
        orthogroup = coregene2orthogroup.get(row.gene, None)
        if orthogroup: 

            # define position indicator for interval to the left
            # downstream means: interval is downstream of core orthogroup
            # --> core orthogroup is on minus strand
            downstream = row.strand == "-"
            right_posind = orthogroup + ("+" if downstream else "-")

            # resolve interval to the left of core gene
            annotation.loc[members, "interval"] = interval
            if left_posind:
                if left_posind < right_posind:
                    intervals.append((interval, left_posind, right_posind))
                else:
                    intervals.append((interval, right_posind, left_posind))
            else:
                intervals.append((interval, right_posind, None))

            interval += 1
            members = []
            left_posind = None
            right_posind = None

            # define position indicator for interval to the right
            downstream = row.strand == "+"
            left_posind = orthogroup + ("+" if downstream else "-")

        # if accessory gene: add to interval
        else: 

            members.append(row.gene)

    # resolve last interval
    annotation.loc[members, "interval"] = interval
    intervals.append((interval, left_posind, None))

    annotation = annotation.reset_index(drop = True)
    genes = annotation.rename(columns = {"seqid": "contig"})
    genes = genes[["gene", "contig", "strand", "start", "end", "interval"]]
    genes = genes[genes["interval"].notnull()]

    intervals = pd.DataFrame.from_records(
        intervals, columns = ["interval", "posind1", "posind2"])

    return(genes, intervals)

def define_positions(intervals: pd.DataFrame) -> dict[str, int]:
    """
    Define genomic positions from intervals. 
    
    :param intervals: Table with columns "posind1" and "posind2". 
    :return: Dictionary with position for every position indicator. 
    """
    # posinds = [o + "+" for o in core] + [o + "-" for o in core]
    posinds = (
        intervals[["posind1", "posind2"]]
        .stack()
        .dropna()
        .unique()
        .tolist()
    )
    posinds = pd.DataFrame(
        {"posind": posinds, "position": range(len(posinds))})
    posinds = posinds.set_index("posind")
    for row in intervals.itertuples():
        if pd.isna(row.posind1): continue
        if pd.isna(row.posind2): continue
        position1 = posinds.at[row.posind1, "position"] 
        position2 = posinds.at[row.posind2, "position"] 
        if position1 != position2:
            posinds.loc[posinds.position == position2, "position"] = \
                position1
    posinds = posinds.reset_index()
    posind2position = dict(zip(posinds["posind"], posinds["position"]))
    return(posind2position)
