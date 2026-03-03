from importlib.metadata import metadata
import logging as lg
from pathlib import Path
import sys

from .assignment import assign_positions
from .assessment import assess_orthogroups
from .input import read_files, read_pangenome

def run_padis(
        assemblies_path: Path, annotation_path: Path, pangenome_file: Path, 
        output_dir: Path, max_length: int = 3000, write_intervals: bool = False, 
        debug: bool = False, threads = 1
    ) -> None:

    args = locals()
    
    lg.basicConfig(
        level = lg.DEBUG if debug else lg.INFO,
        format = "[%(asctime)s] %(levelname)s: %(message)s",
        datefmt = "%d/%m %H:%M:%S"
    )
    
    lg.info("Processing arguments")
    if not debug and output_dir.is_dir(): 
        lg.error("Output folder already exists")
        sys.exit(1)

    # why input file reading already here? --> checking presence of assembly and
    # annotation files for all genomes in pangenome should happen at the start
    assembly_files = read_files(assemblies_path)
    annotation_files = read_files(annotation_path)
    genes = read_pangenome(pangenome_file)
    genomes = genes["genome"].unique()
    try: 
        assembly_files = {g: assembly_files[g] for g in genomes}
    except KeyError as e: 
        lg.error(f"Assembly file for {e.args[0]} not supplied")
        sys.exit(1)
    try: 
        annotation_files = {g: annotation_files[g] for g in genomes}
    except KeyError as e: 
        lg.error(f"Annotation file for {e.args[0]} not supplied")
        sys.exit(1)
    for assembly_file in assembly_files.values():
        if assembly_file.suffix in [".gz", ".bz2", ".xz"]: 
            lg.error("Assemblies should not be compressed")
            sys.exit(1)

    # output dir and log file creation
    # --> should happen after argument checking
    if not output_dir.exists():
        lg.info("Creating output folder")
        output_dir.mkdir(exist_ok = False)
    acc_genes_file = output_dir / "accessory_genes.csv"
    intervals_file = output_dir / "intervals.csv" if write_intervals else None
    acc_orthogroups_file = output_dir / "accessory_orthogroups.csv"
    summary_file = output_dir / "summary.csv"

    # log file initiation
    # --> should happen after argument checking
    log_file = output_dir / "log.txt" 
    log_fileh = lg.FileHandler(log_file, encoding = "utf-8")
    log_fileh.setFormatter(lg.getLogger().handlers[0].formatter)
    lg.getLogger().addHandler(log_fileh)
    meta = metadata('padis')
    lg.info(f"This is {meta['Name']} version {meta['Version']}")
    if debug: lg.info("Running in debugging mode")

    lg.info("Arguments:")
    for arg, value in args.items(): lg.info(f"  {arg}: {value}")

    lg.info("Starting phase 1: position assignment")
    assign_positions(annotation_files, genes, acc_genes_file, intervals_file)

    lg.info("Starting phase 2: orthogroup assessment")
    assess_orthogroups(
        acc_genes_file, assembly_files, acc_orthogroups_file, summary_file, 
        max_length, threads
    )

    lg.info("PADIS out\n")
