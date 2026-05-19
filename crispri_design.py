#!/usr/bin/env python3

import re
import pandas as pd


# =========================
# CAS SYSTEM DEFINITIONS
# =========================

CAS_SYSTEMS = {
    "sp-dcas9": {
        "pam_regex": r"[ACGT]GG",
        "pam_side": "3prime",
        "guide_len": 20,
        "seed_len": 10,
        "seed_side": "3prime"
    },
    "fn-dcas12a": {
        "pam_regex": r"TTT[ACG]",
        "pam_side": "5prime",
        "guide_len": 23,
        "seed_len": 12,
        "seed_side": "5prime"
    },
    "lb-dcas12a": {
        "pam_regex": r"TTT[ACG]",
        "pam_side": "5prime",
        "guide_len": 23,
        "seed_len": 12,
        "seed_side": "5prime"
    }
}


# =========================
# UTILITIES
# =========================

def gc_content(seq):
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)


def reverse_complement(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def count_mismatches(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


# =========================
# GUIDE DISCOVERY
# =========================

def find_guides(sequence, cas_type):

    config = CAS_SYSTEMS[cas_type]
    seq = sequence.upper()
    guides = []

    for match in re.finditer(config["pam_regex"], seq):
        pam_start = match.start()
        pam_seq = match.group()

        if config["pam_side"] == "3prime":
            g_start = pam_start - config["guide_len"]
            g_end = pam_start
        else:
            g_start = pam_start + len(pam_seq)
            g_end = g_start + config["guide_len"]

        if g_start < 0 or g_end > len(seq):
            continue

        guides.append({
            "guide": seq[g_start:g_end],
            "pam": pam_seq,
            "start": g_start,
            "strand": "+"
        })

    return guides


# =========================
# KMER INDEX
# =========================

def build_kmer_index(genome_seq, k):

    index = {}

    for i in range(len(genome_seq) - k + 1):
        kmer = genome_seq[i:i+k]
        if kmer not in index:
            index[kmer] = []
        index[kmer].append(i)

    return index


# =========================
# OFF-TARGET FILTER (FAST, STRAND-AWARE)
# =========================

def passes_offtarget_filter(
    guide,
    genome_seq,
    genome_rc,
    kmer_index,
    cas_type,
    max_mismatches=2,
    seed_mismatch_threshold=1
):

    config = CAS_SYSTEMS[cas_type]

    g_len = config["guide_len"]
    seed_len = config["seed_len"]
    seed_side = config["seed_side"]

    g = guide.upper()

    if seed_side == "3prime":
        guide_seed = g[-seed_len:]
    else:
        guide_seed = g[:seed_len]

    candidate_positions = kmer_index.get(guide_seed, [])

    if len(candidate_positions) == 0:
        return False

    exact_match_count = 0

    for pos in candidate_positions:

        for genome in (genome_seq, genome_rc):

            if pos + g_len > len(genome):
                continue

            window = genome[pos:pos+g_len]

            mismatches = 0
            for a, b in zip(g, window):
                if a != b:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break

            if mismatches == 0:
                exact_match_count += 1
                if exact_match_count > 1:
                    return False

            elif mismatches <= max_mismatches:

                if seed_side == "3prime":
                    window_seed = window[-seed_len:]
                else:
                    window_seed = window[:seed_len]

                seed_mm = count_mismatches(guide_seed, window_seed)

                if seed_mm <= seed_mismatch_threshold:
                    return False

    return True


# =========================
# GC FILTER
# =========================

def filter_guides(guides, min_gc, max_gc):

    filtered = []

    for g in guides:
        gc = gc_content(g["guide"])

        if min_gc <= gc <= max_gc:
            g["GC_content"] = round(gc)
            filtered.append(g)

    return filtered


# =========================
# SELECTION
# =========================

def select_guides(guides, max_guides, gene_len, orf_cutoff):

    guides = sorted(guides, key=lambda x: x["start"])
    cutoff = int(gene_len * orf_cutoff)

    early = [g for g in guides if g["start"] <= cutoff]
    late = [g for g in guides if g["start"] > cutoff]

    selected = []

    for g in early:
        selected.append(g)
        if len(selected) >= max_guides:
            return selected

    if late:
        spacing = max(1, len(late) // max(1, (max_guides - len(selected))))

        for i in range(0, len(late), spacing):
            selected.append(late[i])
            if len(selected) >= max_guides:
                break

    return selected


# =========================
# MAIN PIPELINE
# =========================

def design_guides(sequence,
                  cas_type,
                  genome_seq=None,
                  genome_rc=None,
                  kmer_index=None,
                  max_guides=10,
                  min_gc=20,
                  max_gc=80,
                  orf_cutoff=0.6):

    guides = find_guides(sequence, cas_type)
    guides = filter_guides(guides, min_gc, max_gc)

    if genome_seq is not None:
        filtered = []
        for g in guides:
            if passes_offtarget_filter(
                g["guide"],
                genome_seq,
                genome_rc,
                kmer_index,
                cas_type
            ):
                filtered.append(g)
        guides = filtered

    if len(guides) == 0:
        return pd.DataFrame()

    guides = select_guides(guides, max_guides, len(sequence), orf_cutoff)

    return pd.DataFrame(guides)