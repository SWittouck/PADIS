import gzip
import pandas as pd
from pathlib import Path
from typing import TextIO

def open_smart(file: Path, mode: str = "r") -> TextIO:
    """
    Open a file that may be compressed or not.
    """
    if file.suffix == ".gz":
        return(gzip.open(file, mode = mode))
    else:
        return(file.open(mode = mode))

def read_pangenome(file: Path) -> pd.DataFrame:
    """
    Read a pangenome file in SCARAP format (tsv without column names and with  
    the columns gene, genome and orthogroup). 
    """
    with open_smart(file) as handle:
        pangenome = pd.read_csv(
            file, sep = "\t", names = ["gene", "genome", "orthogroup"])
    return(pangenome)

def read_annotation(file: Path) -> pd.DataFrame:
    """
    Read a gene annotation file in gff format. 
    """
    colnames = ["seqid", "source", "type", "start", "end", "score", "strand", 
        "phase", "attr"]
    types = {"seqid": "str", "source": "str", "type": "str", "start": "int", 
        "end": "int", "score": "float", "strand": "str", "phase": "int", 
        "attr": "str"}
    with open_smart(file) as handle:
        genes = pd.read_csv(
            file, sep = "\t", names = colnames, dtype = types, comment = "#")
    return(genes)

def read_canisgenes(file: Path) -> pd.DataFrame:
    """
    Read a table with candidate insertion sequence genes.
    """
    types = {"start": "Int64", "end": "Int64", "position": "Int64"}
    canisgenes = pd.read_csv(file, dtype = types)
    return(canisgenes)

def read_files(path: Path) -> list[Path]:
    """
    Read list of file paths, either from a file containing the paths or from a 
    directory containing the files themselves. 
    """
    if path.is_dir():
        files = [f for f in path.iterdir()]
    else:
        files = [Path(f.strip()) for f in open(path)]
    return(files)
