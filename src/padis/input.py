import sys
import gzip
import logging as lg
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
    if not file.is_file(): 
        lg.error(f"File not found: {file}")
        sys.exit(1)
    if file.stat().st_size == 0: 
        lg.error("File is empty: {file}")
        sys.exit(1)
    colnames = ["gene", "genome", "orthogroup"]
    pangenome = pd.read_csv(file, sep = "\t", names = colnames)
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
    genes = pd.read_csv(
            file, sep = "\t", names = colnames, dtype = types, comment = "#")
    return(genes)

def read_acc_genes(file: Path) -> pd.DataFrame:
    """
    Read a table with candidate insertion sequence genes.
    """
    types = {"start": "Int64", "end": "Int64", "position": "Int64"}
    acc_genes = pd.read_csv(file, dtype = types)
    return(acc_genes)

def read_files(path: Path) -> dict[str, Path]:
    """
    Read list of file paths, either from a file containing the paths or from a 
    directory containing the files themselves. 
    """
    if not path.exists(): 
        lg.error(f"File or directory not found: {path}")
        sys.exit(1)
    if path.is_dir():
        files = [f for f in path.iterdir()]
    else:
        files = [Path(f.strip()) for f in open(path)]
    if not files: 
        lg.error(f"No files found in {path}")
        sys.exit(1)
    files = {filename_from_path(p): p for p in files}
    return(files)

def filename_from_path(path: Path) -> str:
    """
    Extract filename from path of potentially compressed file. 
    
    :param path: A path. 
    :return: The filename without compression extension (if present) and without
        filetype extension. 
    """
    if path.suffix == ".gz":
        path = path.with_suffix("")
    filename = path.stem
    return(filename)
