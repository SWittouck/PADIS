import argparse
from importlib.metadata import metadata
from pathlib import Path
import sys

from .padis import run_padis

def main() -> None:

    meta = metadata("padis")
    args = parse_arguments(meta) 

    intro = f"Hi, this is {meta["Name"]} version {meta["Version"]}\n"
    print(intro)

    assemblies_path = Path(args.assemblies)
    annotations_path = Path(args.annotations)
    pangenome_file = Path(args.pangenome)
    output_dir = Path(args.outputdir)

    run_padis(
        assemblies_path, annotations_path, pangenome_file, output_dir, 
        args.max_length, args.write_intervals, args.debug, args.threads)

def parse_arguments(meta) -> None:

    parser = argparse.ArgumentParser(
        prog = meta["Name"].lower(), description = meta["Summary"])

    parser.add_argument(
        "assemblies", 
        help = "assemblies in fasta format, either as a file with paths or "
        "as a directory (files should no be compressed)")
    parser.add_argument(
        "annotations", 
        help = "gene coordinates in gff format, either as a file with paths "
        "or as a directory (files can be gzipped)")
    parser.add_argument(
        "pangenome", 
        help = "file with pangenome in SCARAP format: no header; tab-separated; " \
        "columns should be gene, genome and orthogroup (can be gzipped)")
    parser.add_argument(
        "outputdir", help = "directory for output")
    parser.add_argument(
        "-l", "--max_length", type = int, default = 3000,
        help = "maximum length of the insertion sequences in base pairs "
        "[default: 3,000]")
    parser.add_argument(
        "-i", "--write_intervals", action = "store_true",
        help = "write left and right position indicators of intervals")
    parser.add_argument(
        "-b", "--debug", action = "store_true", 
        help = "be extra verbose for debugging and continue in output folder " \
        "if already exists (log file will be appended)")
    parser.add_argument(
        "-t", "--threads", type = int, default = 1,
        help = "number of threads [default: 1]")
    parser.add_argument(
        "-v", "--version", action = "version", 
        version = f"""{meta["Name"]} version {meta["Version"]}""")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args() 

    return(args)

if __name__ == "__main__":
    main() 
