#!/usr/bin/env python3
"""
compute_conservation.py

Inputs:
 - MSA (A3M or FASTA alignment). Recommended: A3M produced by hhblits with query=spike_protein.fasta
 - Reference FASTA (canonical spike): data/raw/spike_protein.fasta
 - Mutation table: data/processed/final_dataset_clean.csv

Outputs:
 - data/processed/conservation_features.csv  (per UniProt position)
 - data/processed/final_dataset_with_conservation.csv (mutation-level merge)
"""

import os
import math
import argparse
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")  # 20 standard AAs


def read_a3m(path):
    """Read A3M and remove lowercase insertion letters.
    Returns dict(header -> sequence string).
    """
    seqs = {}
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line=line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq_lines)
                header = line[1:].split()[0]
                seq_lines = []
            else:
                # remove lowercase letters (insertions), keep uppercase and '-' chars
                cleaned = "".join([c for c in line if (not c.islower())])
                seq_lines.append(cleaned)
        if header is not None:
            seqs[header] = "".join(seq_lines)
    # Verify alignment lengths
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        raise ValueError(f"Sequences in A3M have varying lengths: {sorted(lengths)[:10]}")
    return seqs


def read_msa(msa_path):
    """Auto-detect a3m or fasta MSA. Return (headers, sequences list)."""
    if msa_path.endswith(".a3m") or msa_path.endswith(".A3M"):
        seqdict = read_a3m(msa_path)
        headers = list(seqdict.keys())
        seqs = [seqdict[h] for h in headers]
        return headers, seqs
    else:
        # try AlignIO; common fasta MSA formats should work
        aln = AlignIO.read(msa_path, "fasta")
        headers = [rec.id for rec in aln]
        seqs = [str(rec.seq) for rec in aln]
        L = {len(s) for s in seqs}
        if len(L) != 1:
            raise ValueError("MSA sequences are not all same length")
        return headers, seqs


def henikoff_weights(msa_seqs):
    """Compute Henikoff sequence weights (returns normalized weights summing to 1)."""
    N = len(msa_seqs)
    L = len(msa_seqs[0])
    weights = np.zeros(N, dtype=float)
    for j in range(L):
        col = [msa_seqs[i][j] for i in range(N)]
        # count residues excluding gaps
        counts = {}
        for a in col:
            if a == '-':
                continue
            counts[a] = counts.get(a, 0) + 1
        rj = len(counts)  # number of different residues at col
        if rj == 0:
            continue
        for i, a in enumerate(col):
            if a == '-':
                continue
            nj = counts[a]
            weights[i] += 1.0 / (rj * nj)
    # normalize sum to 1
    total = weights.sum()
    if total == 0:
        # fallback uniform
        return np.ones(N, dtype=float) / N
    weights = weights / total
    return weights


def compute_conservation(msa_headers, msa_seqs, ref_record, rbd_start=None, rbd_end=None, pseudocount=1.0):
    """
    Compute per-reference-position conservation stats.
    ref_record: Biopython SeqRecord with canonical spike sequence (1-based positions).
    If rbd_start/rbd_end provided, will only return positions in that range.
    """
    N = len(msa_seqs)
    L = len(msa_seqs[0])
    # pick reference sequence from MSA: try to find the sequence identical to ref_record (allowing gaps)
    ref_seq = str(ref_record.seq)
    # Find MSA entry that best matches reference (prefer exact uppercase match)
    best_idx = None
    for i, h in enumerate(msa_headers):
        # quick test: compare letters ignoring gaps
        s = msa_seqs[i].replace("-", "")
        if s == ref_seq:
            best_idx = i
            break
    if best_idx is None:
        # fallback: pick first entry
        best_idx = 0

    msa_ref = msa_seqs[best_idx]
    # mapping MSA columns -> UniProt pos (1-based). None if column is an insertion in ref.
    mapping = []
    up_pos = 0
    for j, ch in enumerate(msa_ref):
        if ch != "-":
            up_pos += 1
            mapping.append(up_pos)
        else:
            mapping.append(None)

    # compute Henikoff weights and Neff
    weights = henikoff_weights(msa_seqs)
    neff = 1.0 / np.sum(weights ** 2)

    # background frequencies (simple empirical across all non-gap residues)
    bg_counts = Counter()
    total_non_gaps = 0
    for j in range(L):
        for i in range(N):
            a = msa_seqs[i][j]
            if a == "-":
                continue
            if a not in AA_ORDER:
                continue
            bg_counts[a] += 1
            total_non_gaps += 1
    if total_non_gaps == 0:
        raise ValueError("No valid residues in MSA")
    bg_freq = {aa: bg_counts[aa] / total_non_gaps for aa in AA_ORDER}

    rows = []
    import math
    for j in range(L):
        up = mapping[j]
        if up is None:
            continue
        if rbd_start and (up < rbd_start or up > rbd_end):
            continue
        col = [msa_seqs[i][j] for i in range(N)]
        # counts ignoring gaps and non-standard AAs
        counts = Counter([a for a in col if a in AA_ORDER])
        num_non_gaps = sum(counts.values())
        if num_non_gaps == 0:
            # no data
            p = {aa: 0.0 for aa in AA_ORDER}
            entropy = float("nan")
            conservation = float("nan")
        else:
            # pseudocount Dirichlet smoothing
            denom = num_non_gaps + pseudocount * len(AA_ORDER)
            p = {aa: (counts.get(aa, 0) + pseudocount) / denom for aa in AA_ORDER}
            # Shannon entropy
            entropy = -sum((p[aa] * math.log2(p[aa]) if p[aa] > 0 else 0.0) for aa in AA_ORDER)
            conservation = 1.0 - (entropy / math.log2(len(AA_ORDER)))  # normalized 0..1

        # PSSM log-odds (log2)
        pssm = {aa: (math.log2((p.get(aa, 1e-12) + 1e-12) / (bg_freq.get(aa, 1e-12) + 1e-12))) for aa in AA_ORDER}

        row = {
            "uniprot_pos": up,
            "num_obs": num_non_gaps,
            "neff": float(neff),
            "entropy": float(entropy),
            "conservation": float(conservation),
        }
        # include freq and pssm for all 20 AAs
        for aa in AA_ORDER:
            row[f"freq_{aa}"] = float(counts.get(aa, 0) / num_non_gaps) if num_non_gaps > 0 else 0.0
            row[f"pssm_{aa}"] = float(pssm[aa])
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("uniprot_pos").reset_index(drop=True)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--msa", required=True, help="MSA file (A3M or fasta) with query=spike included")
    parser.add_argument("--ref", default="data/raw/spike_protein.fasta", help="Reference fasta (canonical spike)")
    parser.add_argument("--mut_table", default="data/processed/final_dataset_with_context.csv", help="Mutation-level table with columns: pos, wt, mut, mut_id")
    parser.add_argument("--out", default="data/processed/conservation_features.csv", help="Output per-position conservation CSV")
    parser.add_argument("--merged_out", default="data/processed/final_dataset_with_conservation.csv", help="Mutation-level merged output")
    parser.add_argument("--rbd_start", type=int, default=331, help="Canonical RBD start (UniProt pos)")
    parser.add_argument("--rbd_end", type=int, default=524, help="Canonical RBD end (UniProt pos)")
    args = parser.parse_args()

    # read reference
    ref_record = SeqIO.read(args.ref, "fasta")
    # load MSA
    headers, seqs = read_msa(args.msa)
    print(f"MSA loaded: {len(seqs)} sequences, length {len(seqs[0])}")
    # compute conservation table
    cons_df = compute_conservation(headers, seqs, ref_record, rbd_start=args.rbd_start, rbd_end=args.rbd_end)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    cons_df.to_csv(args.out, index=False)
    print(f"Wrote conservation per-position to {args.out} (rows={len(cons_df)})")

    # Merge with mutation table
    mut_df = pd.read_csv(args.mut_table)
    # ensure pos column exists
    if "pos" not in mut_df.columns:
        if "site_SARS2" in mut_df.columns:
            mut_df = mut_df.rename(columns={"site_SARS2":"pos"})
        else:
            raise ValueError("mut_table must have a 'pos' column (UniProt position)")

    # merge pos -> uniprot_pos
    merged = mut_df.merge(cons_df, left_on="pos", right_on="uniprot_pos", how="left", suffixes=("","_cons"))
    # compute pssm_wt and pssm_mut, freq_wt, freq_mut
    def get_pssm(row, aa):
        col = f"pssm_{aa}"
        if col in row and not pd.isna(row[col]):
            return row[col]
        return float("nan")

    def get_freq(row, aa):
        col = f"freq_{aa}"
        if col in row and not pd.isna(row[col]):
            return row[col]
        return float("nan")

    merged["pssm_wt"] = merged.apply(lambda r: get_pssm(r, r["wt"]) if ("wt" in r and pd.notna(r["wt"])) else float("nan"), axis=1)
    merged["pssm_mut"] = merged.apply(lambda r: get_pssm(r, r["mut"]) if ("mut" in r and pd.notna(r["mut"])) else float("nan"), axis=1)
    merged["freq_wt"] = merged.apply(lambda r: get_freq(r, r["wt"]) if ("wt" in r and pd.notna(r["wt"])) else float("nan"), axis=1)
    merged["freq_mut"] = merged.apply(lambda r: get_freq(r, r["mut"]) if ("mut" in r and pd.notna(r["mut"])) else float("nan"), axis=1)

    merged.to_csv(args.merged_out, index=False)
    print(f"Wrote merged mutation table to {args.merged_out} (rows={len(merged)})")


if __name__ == "__main__":
    main()
