#!/usr/bin/env python3

import argparse
import csv
import itertools
import subprocess
import sys

# =========================
# CONSTANTS
# =========================

REPEAT = "aatttctactcttgtagat"
BOTTOM_SUFFIX = "aaatt"
TOP_PREFIX = "tagat"
TAIL = "a"

JUNCTIONS = [
    "GCTC", "GAGT", "AACG", "TACA",
    "TTCT", "AGAA", "TGGC", "CCCT"
]

MAX_SPACERS = len(JUNCTIONS)

# =========================
# UTILITIES
# =========================

def reverse_complement(seq):
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1]

# check for BsmBI sites when making an array
def validate_no_bsmbi_sites(spacers):
    forbidden = ["gagacg", "cgtctc"]
    failures = []

    for name, seq in spacers:
        seq_lower = seq.lower()
        for site in forbidden:
            if site in seq_lower:
                failures.append((name, seq, site.upper()))

    if failures:
        msg = ["The following spacers are incompatible with multiplex design due to BsmBI sites:"]
        for name, seq, site in failures:
            msg.append(f"  - {name}: contains {site} in {seq}")
        raise ValueError("\n".join(msg))

# =========================
# INPUT PARSER
# =========================

def parse_input(file_path):
    sequences = []
    with open(file_path) as f:
        name = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif line:
                sequences.append((name, line))
    return sequences

# =========================
# SINGLE DESIGN
# =========================

def design_single(spacer, name):

    top    = TOP_PREFIX + spacer + TAIL
    bottom = BOTTOM_SUFFIX + reverse_complement(spacer) + TAIL   # DO NOT include RC(TOP_PREFIX)

    return [
        {"name": f"{name}_top",    "mode": "single", "sequence": top,    "junction": "NONE"},
        {"name": f"{name}_bottom", "mode": "single", "sequence": bottom, "junction": "NONE"},
    ]


# =========================
# MULTIPLEX DESIGN
# =========================

def design_multiplex_crates(spacers):

    n = len(spacers)

    if n < 2:
        raise ValueError("Multiplex mode requires at least 2 spacers.")
    if n > MAX_SPACERS:
        raise ValueError(
            f"Too many spacers: {n} provided but only {MAX_SPACERS} "
            f"junctions are available."
        )

    results = []

    for i, (name, spacer) in enumerate(spacers):

        has_left  = (i > 0)
        has_right = (i < n - 1)

        lci = i - 1
        rci = i

        left_typeA  = has_left  and (lci % 2 == 0)
        left_typeB  = has_left  and (lci % 2 == 1)
        right_typeA = has_right and (rci % 2 == 0)
        right_typeB = has_right and (rci % 2 == 1)

        # ── TOP ──────────────────────────

        if not has_left:
            top = TOP_PREFIX
        elif left_typeA:
            top = REPEAT
        else:
            top = JUNCTIONS[lci] + REPEAT

        top += spacer

        if right_typeA:
            top += JUNCTIONS[rci].lower()
        elif not has_right:
            top += TAIL
        # right_typeB → nothing

        # ── BOTTOM ───────────────────────

        if not has_right:
            bottom = BOTTOM_SUFFIX
        elif right_typeB:
            bottom = reverse_complement(JUNCTIONS[rci]).lower()
        else:
            bottom = ""

        bottom += reverse_complement(spacer)

        if has_left:
            bottom += reverse_complement(REPEAT)

        if not has_left:
            bottom += TAIL
        elif left_typeA:
            bottom += reverse_complement(JUNCTIONS[lci]).lower()
        # left_typeB → nothing

        if i == 0:
            junction_label = JUNCTIONS[0]
        elif i < n - 1:
            junction_label = f"{JUNCTIONS[i-1]}|{JUNCTIONS[i]}"
        else:
            junction_label = JUNCTIONS[i-1]

        results.append({
            "name": f"{name}_top",
            "mode": "multiplex_unit",
            "sequence": top,
            "junction": junction_label
        })

        results.append({
            "name": f"{name}_bottom",
            "mode": "multiplex_unit",
            "sequence": bottom,
            "junction": junction_label
        })

    return results

# =========================
# RNAFOLD
# =========================

def score_structure(seq):
    p = subprocess.Popen(
        ["RNAfold", "--noPS"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, _ = p.communicate(seq)

    try:
        mfe = float(stdout.strip().split("\n")[-1].split("(")[-1].replace(")", ""))
    except:
        mfe = 0.0

    return mfe

# =========================
# BUILD RNA ARRAY (FOR OPTIMIZATION ONLY)
# =========================

def build_full_array(spacers):
    seq = ""
    for i, (_, spacer) in enumerate(spacers):
        if i == 0:
            seq += spacer
        else:
            seq += REPEAT + spacer
    return seq

# =========================
# OPTIMIZATION
# =========================

def optimize_spacer_order(spacers):

    best_score = None
    best_order = None
    rows = []

    print("\n🔬 RNAfold evaluation:\n")

    for idx, perm in enumerate(itertools.permutations(spacers), 1):

        seq = build_full_array(perm)
        mfe = score_structure(seq)

        label = " -> ".join(n for n, _ in perm)
        print(f"{idx:02d} | {label} | MFE: {mfe}")

        rows.append({"order": label, "mfe": mfe})

        if best_score is None or mfe > best_score:
            best_score = mfe
            best_order = perm

    print("\n✅ BEST ORDER:")
    print(" -> ".join(n for n, _ in best_order))
    print(f"MFE: {best_score}\n")

    with open("rnafold_permutations.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["order", "mfe"])
        writer.writeheader()
        writer.writerows(rows)

    return list(best_order)

# =========================
# VISUAL ALIGNMENT
# =========================

def display_alignment(results):
    """
    Display top oligos and correctly offset bottom oligos.
    Bottom strands are reversed and shifted by 4 nt to account for plasmid overhangs.
    """

    top_oligos = [r["sequence"] for r in results if r["name"].endswith("_top")]
    bottom_oligos = [r["sequence"] for r in results if r["name"].endswith("_bottom")]

    # Join with delimiter
    top_line = " | ".join(top_oligos)

    # Reverse ONLY (no complement)
    bottom_reversed = [seq[::-1] for seq in bottom_oligos]
    bottom_line_raw = " | ".join(bottom_reversed)

    # Apply 4 nt offset
    offset = 4
    bottom_line = (" " * offset) + bottom_line_raw

    # Pad top to match length if needed
    max_len = max(len(top_line), len(bottom_line))
    top_line = top_line.ljust(max_len)
    bottom_line = bottom_line.ljust(max_len)

    # Build complementarity line
    match_line = []
    for t, b in zip(top_line, bottom_line):

        if t == " " or b == " ":
            match_line.append(" ")
        elif t == "|":
            match_line.append("|")
        elif reverse_complement(t) == b:
            match_line.append(":")
        else:
            match_line.append(" ")

    match_line = "".join(match_line)

    print("\n🔬 Oligo Alignment:\n")
    print("TOP    : ", top_line)
    print("MATCH  : ", match_line)
    print("BOTTOM : ", bottom_line)
    print("\nLegend: ':' = complementary basepair\n")





# =========================
# OUTPUT
# =========================

def write_csv(results, output_file):
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["name", "mode", "sequence", "junction"]
        )
        writer.writeheader()
        for r in results:
            writer.writerow(r)

# =========================
# MAIN
# =========================

def main():
    parser = argparse.ArgumentParser(
        description=(
                "Design dCas12 CRISPRi oligonucleotides for single or "
                f"multiplex (2–{MAX_SPACERS} spacers) array assembly."
            )
        )

    parser.add_argument("input",
        help="FASTA file (or plain text, one spacer per line) of protospacer sequences.")
    parser.add_argument("--mode", choices=["single", "multiplex"], required=True,
        help="'single': independent cloning oligos.  'multiplex': repeat-spacer array oligos.")
    parser.add_argument("--output", default="oligos.csv",
        help="Output CSV filename (default: oligos.csv).")
    parser.add_argument("--optimize", action="store_true")

    args = parser.parse_args()

    spacers = parse_input(args.input)
    results = []

    try:
        if args.mode == "multiplex":
            validate_no_bsmbi_sites(spacers)

    except ValueError as e:
        print("\n[ERROR] Invalid spacer(s) detected:", file=sys.stderr)
        print(e, file=sys.stderr)
        print("\nFix the sequences and rerun.", file=sys.stderr)
        sys.exit(1)
    
# SINGLE MODE
    if args.mode == "single":
        for name, spacer in spacers:
            results.extend(design_single(spacer, name))

    # MULTIPLEX MODE
    elif args.mode == "multiplex":

        if args.optimize:
            spacers = optimize_spacer_order(spacers)

        results = design_multiplex_crates(spacers)

    write_csv(results, args.output)
    
    if args.mode == "multiplex":
        display_alignment(results)



if __name__ == "__main__":
    main()
