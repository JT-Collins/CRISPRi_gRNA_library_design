#!/usr/bin/env python3

"""
CRISPRi gRNA Library Designer (Deterministic, K-mer Accelerated)

Run:
    python design_library.py config.txt
"""

import os
import sys
import pandas as pd
import configparser
import re

from crispri_design import (
    design_guides,
    reverse_complement,
    build_kmer_index,
    CAS_SYSTEMS
)


# =========================
# CONFIG PARSING
# =========================

def parse_config(config_file):
    """
    Parse config file into dictionary.
    """

    if not os.path.exists(config_file):
        raise FileNotFoundError(f"Config file not found: {config_file}")

    config = configparser.ConfigParser()
    config.read(config_file)

    if len(config.sections()) != 1:
        raise ValueError("Config must contain exactly one section")

    section = config.sections()[0]

    return {k: v for k, v in config.items(section)}


# =========================
# FASTA PARSING
# =========================

import re

def parse_header_metadata(header_line):
    """
    Extract metadata from NCBI-style CDS FASTA headers.

    Example header:
    >lcl|... [gene=dnaA] [locus_tag=...] [protein=...] [protein_id=...]

    Returns dict with keys:
        gene, locus_tag, protein, protein_id
    Missing fields are set to None.
    """

    meta = {
        "gene": None,
        "locus_tag": None,
        "protein": None,
        "protein_id": None
    }

    # extract [key=value] terms
    for key in meta.keys():
        match = re.search(rf"\[{key}=([^\]]+)\]", header_line)
        if match:
            meta[key] = match.group(1)

    return meta


def read_fasta(filepath):
    """
    Read FASTA and return:
    sequences: {gene_id: sequence}
    metadata:  {gene_id: {gene, locus_tag, protein, protein_id}}
    """

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    sequences = {}
    metadata = {}

    current_id = None
    current_meta = None

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                header = line[1:]

                # fallback ID: first token (existing behavior)
                current_id = header.split()[0]

                sequences[current_id] = ""
                current_meta = parse_header_metadata(header)
                metadata[current_id] = current_meta

            else:
                sequences[current_id] += line.upper()

    return sequences, metadata


# =========================
# GENOME LOADING
# =========================

def read_genome(filepath):
    """
    Load genome FASTA into a single concatenated string.
    """

    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Genome file not found: {filepath}")

    genome = []

    with open(filepath) as f:
        for line in f:
            if not line.startswith(">"):
                genome.append(line.strip().upper())

    genome_seq = "".join(genome)

    if len(genome_seq) == 0:
        raise ValueError("Genome file contains no sequence")

    return genome_seq


# =========================
# MAIN PIPELINE
# =========================

def design_library(config):

    # -------------------------
    # PARAMETERS
    # -------------------------

    cas_type = config["crispri_system"].lower()
    fasta_path = config["target_fasta"]

    max_guides = int(config["grna_number"])
    min_gc = float(config["min_gc_content"])
    max_gc = float(config["max_gc_content"])
    orf_cutoff = float(config["orf_cutoff"])

    prefix = config["prefix"]

    if cas_type not in CAS_SYSTEMS:
        raise ValueError(f"Unsupported CAS system: {cas_type}")

    # -------------------------
    # OUTPUT DIRECTORY
    # -------------------------

    print(f"[INFO] Creating output directory: {prefix}")
    os.makedirs(prefix, exist_ok=True)

    # -------------------------
    # LOAD TARGETS
    # -------------------------

    print(f"[INFO] Loading FASTA: {fasta_path}")
    sequences, metadata = read_fasta(fasta_path)
    print(f"[INFO] Loaded {len(sequences)} sequences")

    # -------------------------
    # LOAD GENOME + INDEX
    # -------------------------

    genome_seq = None
    genome_rc = None
    kmer_index = None

    if "genome" in config and config["genome"] != "":

        genome_path = config["genome"]

        print(f"[INFO] Loading genome: {genome_path}")
        genome_seq = read_genome(genome_path)

        print(f"[INFO] Genome length: {len(genome_seq):,} bp")

        print("[INFO] Computing reverse complement")
        genome_rc = reverse_complement(genome_seq)

        k = CAS_SYSTEMS[cas_type]["seed_len"]

        print(f"[INFO] Building k-mer index (k={k})")
        kmer_index = build_kmer_index(genome_seq, k)

        print(f"[INFO] K-mer index size: {len(kmer_index):,}")

    else:
        print("[INFO] No genome provided → skipping off-target filtering")

    # -------------------------
    # DESIGN LOOP
    # -------------------------

    all_guides = []

    for gene_id, seq in sequences.items():

        print(f"[INFO] Designing guides for {gene_id}")

        df = design_guides(
            seq,
            cas_type,
            genome_seq=genome_seq,
            genome_rc=genome_rc,
            kmer_index=kmer_index,
            max_guides=max_guides,
            min_gc=min_gc,
            max_gc=max_gc,
            orf_cutoff=orf_cutoff
        )

        print(f"[INFO] {len(df)} guides retained")

        if df.empty:
            print(f"[WARNING] No valid guides for {gene_id}")
            continue

        meta = metadata.get(gene_id, {})

        df["gene_id"] = gene_id

        # prefer parsed metadata if present; fallback to gene_id
        df["gene"] = meta.get("gene") or gene_id
        df["locus_tag"] = meta.get("locus_tag") or gene_id
        df["protein"] = meta.get("protein") or gene_id
        df["protein_id"] = meta.get("protein_id") or gene_id
        all_guides.append(df)

    # -------------------------
    # CHECK OUTPUT
    # -------------------------

    if not all_guides:
        raise RuntimeError("No guides generated for any gene")

    library_df = pd.concat(all_guides, ignore_index=True)

    # -------------------------
    # WRITE OUTPUT
    # -------------------------

    output_csv = os.path.join(prefix, "gRNA_library.csv")

    print(f"[INFO] Writing output: {output_csv}")
    library_df.to_csv(output_csv, index=False)

    print("\n✅ SUCCESS: CRISPRi library design complete")


# =========================
# ENTRY POINT
# =========================

def main():

    if len(sys.argv) != 2:
        print("Usage: python design_library.py config.txt")
        sys.exit(1)

    config_file = sys.argv[1]

    print("[INFO] Parsing config")
    config = parse_config(config_file)

    print("[INFO] Running design pipeline")
    design_library(config)


if __name__ == "__main__":
    main()