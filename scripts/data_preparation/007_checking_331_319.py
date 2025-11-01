# # # # # #!/usr/bin/env python3
# # # # # """
# # # # # build_uniprot_to_pdb_map.py

# # # # # Inputs (edit paths below):
# # # # #  - uniprot_fasta: full spike FASTA (UniProt P0DTC2)
# # # # #  - fragment_seq_file: RBD fragment fasta/txt you pasted (or will be extracted from UniProt positions)
# # # # #  - pdb_file: repaired PDB (e.g. data/external/6m0j_Repair.pdb)
# # # # #  - mutations_input: plain text file with one mutation per line like 'N331A' (or CSV, simple mode)

# # # # # Outputs:
# # # # #  - data/processed/uniprot_pdb_mapping.csv
# # # # #  - data/processed/individual_list.txt   (FoldX-ready)
# # # # #  - data/processed/mapping_log.txt
# # # # # """

# # # # # from pathlib import Path
# # # # # from Bio.PDB import PDBParser, Polypeptide
# # # # # from Bio import SeqIO
# # # # # import re
# # # # # import csv

# # # # # # ---------- CONFIG: change these file paths ----------
# # # # # uniprot_fasta = Path("data/raw/spike_protein.fasta")   # full UniProt FASTA (P0DTC2)
# # # # # fragment_seq_file = Path("data/raw/rbd_region.txt")  # or .txt (the fragment you pasted)
# # # # # pdb_file = Path("data/external/6m0j_Repair.pdb")
# # # # # mutations_input = Path("data/external/mutations_input.txt")   # lines like N331A or CSV (simple)
# # # # # out_map = Path("data/processed/uniprot_pdb_mapping.csv")
# # # # # out_foldx = Path("data/processed/individual_list.txt")
# # # # # out_log = Path("data/processed/mapping_log.txt")
# # # # # # ----------------------------------------------------

# # # # # out_map.parent.mkdir(parents=True, exist_ok=True)

# # # # # def read_fasta_or_text(path: Path) -> str:
# # # # #     txt = path.read_text().strip()
# # # # #     if txt.startswith(">"):
# # # # #         # parse first FASTA record
# # # # #         seq = str(next(SeqIO.parse(str(path), "fasta")).seq)
# # # # #     else:
# # # # #         seq = "".join(txt.split())  # remove newline/spaces
# # # # #     return seq.upper()

# # # # # # read sequences
# # # # # if not uniprot_fasta.exists():
# # # # #     raise FileNotFoundError(f"UniProt FASTA not found at {uniprot_fasta}. Please provide full spike FASTA.")
# # # # # full_uniprot = read_fasta_or_text(uniprot_fasta)
# # # # # frag_seq = read_fasta_or_text(fragment_seq_file)

# # # # # # parse PDB and build chain sequences + residue number lists
# # # # # parser = PDBParser(QUIET=True)
# # # # # structure = parser.get_structure("s", str(pdb_file))
# # # # # model = structure[0]

# # # # # chain_info = {}  # chain -> {seq:..., res_list: [(resnum, icode), ...]}
# # # # # for chain in model:
# # # # #     seq = []
# # # # #     res_list = []
# # # # #     for res in chain:
# # # # #         if not Polypeptide.is_aa(res, standard=True):
# # # # #             continue
# # # # #         three = res.get_resname()
# # # # #         try:
# # # # #             one = Polypeptide.three_to_one(three)
# # # # #         except Exception:
# # # # #             one = "X"
# # # # #         seq.append(one)
# # # # #         # store residue number and insertion code
# # # # #         res_list.append((res.get_id()[1], res.get_id()[2]))
# # # # #     chain_info[chain.id] = {"seq": "".join(seq), "res_list": res_list}

# # # # # # search for fragment in all chains (exact substring)
# # # # # matches = []
# # # # # for chain_id, info in chain_info.items():
# # # # #     idx = info["seq"].find(frag_seq)
# # # # #     if idx != -1:
# # # # #         matches.append({"chain": chain_id, "seq_idx": idx, "pdb_start_resnum": info["res_list"][idx][0]})
# # # # # # if none found, try searching for smaller unique seed (first 20 aa)
# # # # # if not matches:
# # # # #     seed = frag_seq[:20]
# # # # #     for chain_id, info in chain_info.items():
# # # # #         idx = info["seq"].find(seed)
# # # # #         if idx != -1:
# # # # #             matches.append({"chain": chain_id, "seq_idx": idx, "pdb_start_resnum": info["res_list"][idx][0], "note": "seed-match"})
# # # # # # Prepare logging + mapping table
# # # # # log_lines = []
# # # # # mapping_rows = []  # each: uni_pos, uni_aa, pdb_chain, pdb_resnum, pdb_aa, status
# # # # # foldx_lines = []
# # # # # if not matches:
# # # # #     log_lines.append("ERROR: fragment sequence not found in any PDB chain. Cannot auto-map.")
# # # # # else:
# # # # #     log_lines.append(f"Found {len(matches)} fragment match(es):")
# # # # #     for m in matches:
# # # # #         log_lines.append(str(m))

# # # # #     # For each match, create mapping for fragment positions
# # # # #     # We'll try to infer the UniProt offset by aligning the fragment to full UniProt:
# # # # #     # find where the fragment occurs in the full UniProt sequence
# # # # #     frag_positions_in_uniprot = []
# # # # #     uni_idx = full_uniprot.find(frag_seq)
# # # # #     if uni_idx != -1:
# # # # #         # uni_idx is 0-based index in UniProt sequence. UniProt residue numbers are 1-based.
# # # # #         frag_start_uni = uni_idx + 1
# # # # #         log_lines.append(f"Fragment found in UniProt at 1-based position {frag_start_uni}")
# # # # #     else:
# # # # #         # fallback: ask user assumption: commonly RBD starts at UniProt 319. We'll attempt to find best match by scanning
# # # # #         # but avoid interactive ask; try to find best local alignment by substring matches
# # # # #         # Try to find a 10-mer from fragment in UniProt to infer start
# # # # #         found = False
# # # # #         for k in range(0, len(frag_seq)-9):
# # # # #             sub = frag_seq[k:k+10]
# # # # #             p = full_uniprot.find(sub)
# # # # #             if p != -1:
# # # # #                 frag_start_uni = p - k + 1
# # # # #                 log_lines.append(f"Estimated fragment UniProt start by 10-mer match: {frag_start_uni}")
# # # # #                 found = True
# # # # #                 break
# # # # #         if not found:
# # # # #             log_lines.append("Could not locate fragment in UniProt sequence automatically. Please ensure fragment_seq matches UniProt content.")
# # # # #             frag_start_uni = None

# # # # #     # build the mapping using the first match (common case)
# # # # #     match = matches[0]
# # # # #     chain_id = match["chain"]
# # # # #     seq_idx = match["seq_idx"]
# # # # #     res_list = chain_info[chain_id]["res_list"]
# # # # #     pdb_seq = chain_info[chain_id]["seq"]

# # # # #     if frag_start_uni is None:
# # # # #         log_lines.append("Aborting detailed mapping due to inability to locate fragment in UniProt.")
# # # # #     else:
# # # # #         # For each residue in fragment, assign mapping
# # # # #         for i in range(len(frag_seq)):
# # # # #             uni_pos = frag_start_uni + i
# # # # #             uni_aa = frag_seq[i]
# # # # #             res_list_idx = seq_idx + i
# # # # #             if res_list_idx < 0 or res_list_idx >= len(res_list):
# # # # #                 # PDB has no residue at this index
# # # # #                 mapping_rows.append((uni_pos, uni_aa, chain_id, None, None, "MISSING_IN_PDB"))
# # # # #                 continue
# # # # #             pdb_resnum, icode = res_list[res_list_idx]
# # # # #             pdb_aa = pdb_seq[res_list_idx]
# # # # #             status = "OK" if pdb_aa == uni_aa else "AA_MISMATCH"
# # # # #             mapping_rows.append((uni_pos, uni_aa, chain_id, pdb_resnum, pdb_aa, status))

# # # # # # write mapping CSV
# # # # # with open(out_map, "w", newline="") as fh:
# # # # #     writer = csv.writer(fh)
# # # # #     writer.writerow(["uni_pos", "uni_aa", "pdb_chain", "pdb_resnum", "pdb_aa", "status"])
# # # # #     for row in mapping_rows:
# # # # #         writer.writerow(row)

# # # # # # Now read mutation input (lines like N331A) and convert using mapping table
# # # # # # load mapping into dict: (uni_pos) -> (pdb_chain, pdb_resnum, pdb_aa, status)
# # # # # mapdict = {r[0]: {"uni_aa": r[1], "chain": r[2], "pdb_resnum": r[3], "pdb_aa": r[4], "status": r[5]} for r in mapping_rows}

# # # # # if mutations_input.exists():
# # # # #     for line in mutations_input.read_text().splitlines():
# # # # #         s = line.strip().rstrip(";")
# # # # #         if not s: continue
# # # # #         # accept formats like N331A or N 331 A
# # # # #         m = re.match(r'^([A-Za-z])\s*([0-9]+)\s*([A-Za-z\*])$', s)
# # # # #         if not m:
# # # # #             log_lines.append(f"Unrecognized mutation format (skipping): {s}")
# # # # #             continue
# # # # #         orig, pos_str, new = m.group(1).upper(), int(m.group(2)), m.group(3).upper()
# # # # #         if pos_str not in mapdict and pos_str==str(pos_str): pass
# # # # #         entry = mapdict.get(pos_str) if False else mapdict.get(pos_str)  # placeholder (unused)
# # # # #         entry = mapdict.get(int(pos_str))
# # # # #         if entry is None:
# # # # #             # not in fragment mapping
# # # # #             log_lines.append(f"SKIP: mutation {orig}{pos_str}{new} -> residue not within mapped fragment (or mapping missing).")
# # # # #             continue
# # # # #         if entry["status"] != "OK":
# # # # #             log_lines.append(f"SKIP: mutation {orig}{pos_str}{new} -> mapping status {entry['status']} (pdb:{entry['pdb_aa']}{entry['pdb_resnum']})")
# # # # #             continue
# # # # #         # everything ok: build FoldX mutation string
# # # # #         foldx_line = f"{orig}{entry['pdb_resnum']}{new};"
# # # # #         foldx_lines.append(foldx_line)
# # # # #         log_lines.append(f"MAPPED: {orig}{pos_str}{new} -> {foldx_line} (chain {entry['chain']})")
# # # # # else:
# # # # #     log_lines.append(f"No mutation input file found at {mutations_input}; skipping FoldX mutation generation.")

# # # # # # write foldx mutation file and log
# # # # # out_foldx.write_text("\n".join(foldx_lines) + ("\n" if foldx_lines else ""))
# # # # # out_log.write_text("\n".join(log_lines) + ("\n" if log_lines else ""))

# # # # # print("Wrote mapping to:", out_map)
# # # # # print("Wrote foldx list to:", out_foldx)
# # # # # print("Wrote log to:", out_log)




# # # # # from Bio.PDB import PDBParser

# # # # # pdb_file = "data/external/6m0j_Repair.pdb"
# # # # # parser = PDBParser(QUIET=True)
# # # # # structure = parser.get_structure("spike", pdb_file)
# # # # # model = structure[0]

# # # # # for chain in model:
# # # # #     print(f"Chain {chain.id}")
# # # # #     for res in chain.get_residues():
# # # # #         print(res.get_resname(), res.get_id())
# # # # #         break  # only first residue per chain




# # # # # from Bio.PDB import PDBParser
# # # # # from Bio.PDB.Polypeptide import is_aa

# # # # # pdb_file = "data/external/6m0j_Repair.pdb"
# # # # # parser = PDBParser(QUIET=True)
# # # # # structure = parser.get_structure("spike", pdb_file)
# # # # # model = structure[0]

# # # # # # Chain E = Spike RBD, Chain A = ACE2
# # # # # spike_chain = model["E"]
# # # # # ace2_chain = model["A"]

# # # # # print("Spike chain (E):")
# # # # # for res in spike_chain:
# # # # #     if is_aa(res, standard=True):
# # # # #         print(res.get_resname(), res.get_id()[1])
# # # # #         break

# # # # # print("\nACE2 chain (A):")
# # # # # for res in ace2_chain:
# # # # #     if is_aa(res, standard=True):
# # # # #         print(res.get_resname(), res.get_id()[1])
# # # # #         break




# # # # from pathlib import Path

# # # # def clean_pdb(infile, outfile):
# # # #     with open(infile) as f, open(outfile, "w") as out:
# # # #         for line in f:
# # # #             if line.startswith("ATOM") or line.startswith("HETATM"):
# # # #                 out.write(line)
# # # #             elif line.startswith("CONECT"):  # keep CONECT if present
# # # #                 out.write(line)

# # # # # Use Path to keep it OS-agnostic
# # # # infile = Path(r"C:\Users\avant\OneDrive\Desktop\VSCode\C++\SpikeMutations\data\external\6m0j_Repair.pdb")
# # # # outfile = Path(r"C:\Users\avant\OneDrive\Desktop\VSCode\C++\SpikeMutations\data\external\6m0j_clean.pdb")

# # # # clean_pdb(infile, outfile)
# # # # print("Cleaned PDB written to:", outfile)

# # # # from pdbfixer import PDBFixer
# # # # from openmm.app import PDBFile

# # # # fixer = PDBFixer(filename="data/external/6m0j_Repair.pdb")
# # # # # remove heterogens (water, ions, ligands) if you want
# # # # fixer.removeHeterogens(keepWater=False)
# # # # PDBFile.writeFile(fixer.topology, fixer.positions, open("data/external/6m0j_clean.pdb", "w"))

# # # print("-"*100)
# # # from Bio.PDB import PDBParser

# # # parser = PDBParser(QUIET=True)
# # # structure = parser.get_structure("spike", "data/external/6m0j_clean.pdb")

# # # # for model in structure:
# # # #     print("Model", model.id, "chains:", [c.id for c in model])
# # # # print("-"*100)
# # # model = structure[0]
# # # # Detect the chain containing RBD residues (331–531)
# # # rbd_chain = None
# # # for chain in model:
# # #     res_nums = [res.id[1] for res in chain if res.id[0] == " "]
# # #     if 331 in res_nums and 531 in res_nums:
# # #         rbd_chain = chain
# # #         break

# # # if rbd_chain is None:
# # #     raise ValueError("❌ Could not find chain with RBD (331–531)")
# # # else:
# # #     print(f"✅ RBD chain found: {rbd_chain.id}")

# # # """
# # # Build structural context features for SARS-CoV-2 RBD from 6M0J.pdb

# # # Outputs:
# # #     data/processed/context_features.csv
# # # """

# # # import os
# # # import pandas as pd
# # # from Bio.PDB import PDBParser, DSSP
# # # from Bio.PDB.Polypeptide import is_aa
# # # from Bio.SeqUtils import seq1

# # # # -----------------------------
# # # # Config
# # # # -----------------------------
# # # pdb_file = "data/external/6m0j_clean.pdb"
# # # output_csv = "data/processed/context_features.csv"

# # # # UniProt RBD region (319–541)
# # # RBD_START, RBD_END = 319, 541
# # # RBM_RESIDUES = range(438, 507)  # Receptor binding motif residues

# # # # Offsets:
# # # # Chain E starts at THR333 in PDB, but UniProt RBD starts at ASN331
# # # OFFSET = +14  # UniProt -> PDB numbering shift

# # # # -----------------------------
# # # # Load structure
# # # # -----------------------------
# # # parser = PDBParser(QUIET=True)
# # # structure = parser.get_structure("complex", pdb_file)
# # # model = structure[0]

# # # spike_chain = model["E"]  # Spike RBD
# # # ace2_chain = model["A"]   # ACE2 receptor

# # # # -----------------------------
# # # # Collect RBD residues
# # # # -----------------------------
# # # rbd_residues = [
# # #     res for res in spike_chain if is_aa(res, standard=True)
# # #     and (RBD_START + OFFSET) <= res.get_id()[1] <= (RBD_END + OFFSET)
# # # ]

# # # # -----------------------------
# # # # DSSP: secondary structure + RSA
# # # # -----------------------------
# # # dssp = DSSP(model, pdb_file, dssp="mkdssp")

# # # residue_numbers, aa_list, ss_list, rsa_list = [], [], [], []
# # # for key in dssp.keys():
# # #     chain_id, res_id = key[0], key[1][1]
# # #     if chain_id != spike_chain.id:
# # #         continue
# # #     if not ((RBD_START + OFFSET) <= res_id <= (RBD_END + OFFSET)):
# # #         continue
# # #     residue_numbers.append(res_id - OFFSET)  # back to UniProt numbering
# # #     aa_list.append(dssp[key][1])
# # #     ss_list.append(dssp[key][2])
# # #     rsa_list.append(dssp[key][3])

# # # # -----------------------------
# # # # Distance to ACE2
# # # # -----------------------------
# # # ace2_atoms = [atom for atom in ace2_chain.get_atoms()]
# # # dist_to_ace2 = []
# # # for res in rbd_residues:
# # #     atoms = list(res.get_atoms())
# # #     if len(atoms) == 0:
# # #         dist_to_ace2.append(None)
# # #         continue
# # #     min_dist = min(atom - ace_atom for atom in atoms for ace_atom in ace2_atoms)
# # #     dist_to_ace2.append(min_dist)

# # # # -----------------------------
# # # # Distance to RBM
# # # # -----------------------------
# # # dist_to_rbm = []
# # # for res in rbd_residues:
# # #     res_id_uniprot = res.get_id()[1] - OFFSET
# # #     if res_id_uniprot in RBM_RESIDUES:
# # #         dist_to_rbm.append(0)
# # #     else:
# # #         dist_to_rbm.append(min(abs(res_id_uniprot - r) for r in RBM_RESIDUES))

# # # # -----------------------------
# # # # Merge into DataFrame
# # # # -----------------------------
# # # df = pd.DataFrame({
# # #     "residue": residue_numbers,         # UniProt numbering
# # #     "aa": aa_list,
# # #     "secondary_structure": ss_list,
# # #     "rsa": rsa_list,
# # #     "dist_to_ACE2": dist_to_ace2,
# # #     "dist_to_RBM": dist_to_rbm,
# # # })

# # # # -----------------------------
# # # # Save
# # # # -----------------------------
# # # os.makedirs(os.path.dirname(output_csv), exist_ok=True)
# # # df.to_csv(output_csv, index=False)

# # # print(f"Context features saved to {output_csv}")
# # # print(df.head())


# # from Bio.PDB import PDBParser, DSSP

# # # Path to your cleaned PDB
# # pdb_file = "data/external/6m0j_clean.pdb"

# # # Load structure
# # parser = PDBParser(QUIET=True)
# # structure = parser.get_structure("spike", pdb_file)
# # model = structure[0]

# # print("-" * 100)
# # print("Model", model.id, "chains:", [c.id for c in model])
# # print("-" * 100)

# # # -------------------------------------------------------------------
# # # Detect the RBD chain dynamically (residues 331–531 must exist there)
# # # -------------------------------------------------------------------
# # rbd_chain = None
# # for chain in model:
# #     res_nums = [res.id[1] for res in chain if res.id[0] == " "]
# #     if 331 in res_nums and 531 in res_nums:
# #         rbd_chain = chain
# #         break

# # if rbd_chain is None:
# #     raise ValueError("❌ Could not find chain with RBD (331–531)")
# # else:
# #     print(f"✅ RBD chain found: {rbd_chain.id}")

# # # Use this chain for downstream analysis
# # spike_chain = rbd_chain

# # # -------------------------------------------------------------------
# # # Example: get sequence of residues 331–531
# # # -------------------------------------------------------------------
# # rbd_residues = [
# #     res for res in spike_chain 
# #     if res.id[0] == " " and 331 <= res.id[1] <= 531
# # ]

# # rbd_seq = "".join(res.resname for res in rbd_residues)
# # print(f"Extracted {len(rbd_residues)} residues for RBD (331–531).")

# # # -------------------------------------------------------------------
# # # Run DSSP on the model to annotate secondary structure
# # # -------------------------------------------------------------------
# # try:
# #     dssp = DSSP(model, pdb_file, dssp="mkdssp")  # requires DSSP installed
# #     print("✅ DSSP run successfully")
# # except Exception as e:
# #     print("❌ DSSP failed:", e)

# # # -------------------------------------------------------------------
# # # Example: print out residue numbers + DSSP assignment
# # # -------------------------------------------------------------------
# # if "dssp" in locals():
# #     for res in rbd_residues[:10]:  # just first 10 for demo
# #         key = (spike_chain.id, res.id)
# #         if key in dssp:
# #             aa, ss, acc = dssp[key][1], dssp[key][2], dssp[key][3]
# #             print(f"Residue {res.id[1]} {aa}: SS={ss}, ASA={acc}")


# from Bio.PDB import PDBParser, DSSP
# import csv

# # Path to your cleaned PDB
# pdb_file = "data/external/6m0j_clean.pdb"

# # Load structure
# parser = PDBParser(QUIET=True)
# structure = parser.get_structure("spike", pdb_file)
# model = structure[0]

# print("-" * 100)
# print("Model", model.id, "chains:", [c.id for c in model])
# print("-" * 100)

# # -------------------------------------------------------------------
# # Detect the RBD chain dynamically (residues 331–531 must exist there)
# # -------------------------------------------------------------------
# rbd_chain = None
# for chain in model:
#     res_nums = [res.id[1] for res in chain if res.id[0] == " "]
#     if 331 in res_nums and 531 in res_nums:
#         rbd_chain = chain
#         break

# if rbd_chain is None:
#     raise ValueError("❌ Could not find chain with RBD (331–531)")
# else:
#     print(f"✅ RBD chain found: {rbd_chain.id}")

# # Use this chain for downstream analysis
# spike_chain = rbd_chain

# # -------------------------------------------------------------------
# # Collect RBD residues (331–531)
# # -------------------------------------------------------------------
# rbd_residues = [
#     res for res in spike_chain 
#     if res.id[0] == " " and 331 <= res.id[1] <= 531
# ]

# print(f"Extracted {len(rbd_residues)} residues for RBD (331–531).")

# # -------------------------------------------------------------------
# # Run DSSP on the model to annotate secondary structure
# # -------------------------------------------------------------------
# try:
#     dssp = DSSP(model, pdb_file, dssp="mkdssp")  # requires DSSP installed
#     print("✅ DSSP run successfully")
# except Exception as e:
#     print("❌ DSSP failed:", e)
#     dssp = None

# # -------------------------------------------------------------------
# # Save results to CSV
# # -------------------------------------------------------------------
# if dssp:
#     out_csv = "data/processed/rbd_dssp_features.csv"
#     with open(out_csv, "w", newline="") as f:
#         writer = csv.writer(f)
#         writer.writerow(["position", "amino_acid", "SS", "ASA"])
        
#         for res in rbd_residues:
#             key = (spike_chain.id, res.id)
#             if key in dssp:
#                 aa, ss, acc = dssp[key][1], dssp[key][2], dssp[key][3]
#                 writer.writerow([res.id[1], aa, ss, acc])

#     print(f"DSSP features saved to {out_csv}")


"""
Build structural context features for SARS-CoV-2 RBD (319–541) from 6M0J.pdb

Outputs:
    data/processed/context_features.csv
"""

import os
import pandas as pd
from Bio.PDB import PDBParser, DSSP
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

# -----------------------------
# Config
# -----------------------------
pdb_file = "data/external/6m0j_clean.pdb"   # use cleaned/repaired file
output_csv = "data/processed/context_features.csv"

# SARS-CoV-2 RBD residue range
RBD_START, RBD_END = 319, 541
ACE2_CHAIN = "A"  # ACE2 receptor chain in 6M0J
RBM_RESIDUES = range(438, 507)  # receptor binding motif

# -----------------------------
# Utility
# -----------------------------
def try_three_to_one(resname: str) -> str:
    """Convert 3-letter AA code to 1-letter, fallback to X."""
    try:
        return seq1(resname)
    except Exception:
        return "X"

# -----------------------------
# Load structure
# -----------------------------
parser = PDBParser(QUIET=True)
structure = parser.get_structure("complex", pdb_file)
model = structure[0]

# Identify RBD chain automatically
rbd_chain = None
for chain in model:
    res_ids = [res.get_id()[1] for res in chain if is_aa(res, standard=True)]
    if any(RBD_START <= r <= RBD_END for r in res_ids):
        rbd_chain = chain
        break

if rbd_chain is None:
    raise ValueError("❌ Could not find RBD chain (319–541) in PDB")
print(f"✅ RBD chain detected: {rbd_chain.id}")

# Get ACE2 chain
ace2_chain = model[ACE2_CHAIN]

# Collect RBD residues
rbd_residues = [
    res for res in rbd_chain
    if is_aa(res, standard=True) and RBD_START <= res.get_id()[1] <= RBD_END
]

# -----------------------------
# DSSP: secondary structure + RSA
# -----------------------------
dssp = DSSP(model, pdb_file, dssp="mkdssp")

residue_numbers, aa_list, ss_list, rsa_list = [], [], [], []
for key in dssp.keys():
    chain_id, res_id = key[0], key[1][1]
    if chain_id != rbd_chain.id:
        continue
    if not (RBD_START <= res_id <= RBD_END):
        continue
    residue_numbers.append(res_id)
    aa_list.append(dssp[key][1])
    ss_list.append(dssp[key][2])
    rsa_list.append(dssp[key][3])

# -----------------------------
# Distance to ACE2
# -----------------------------
ace2_atoms = [atom for atom in ace2_chain.get_atoms()]
dist_to_ace2 = []
for res in rbd_residues:
    atoms = list(res.get_atoms())
    if not atoms:
        dist_to_ace2.append(None)
        continue
    min_dist = min(atom - ace_atom for atom in atoms for ace_atom in ace2_atoms)
    dist_to_ace2.append(min_dist)

# -----------------------------
# Distance to RBM (sequence distance)
# -----------------------------
dist_to_rbm = []
for res in rbd_residues:
    res_id = res.get_id()[1]
    if res_id in RBM_RESIDUES:
        dist_to_rbm.append(0)
    else:
        dist_to_rbm.append(min(abs(res_id - r) for r in RBM_RESIDUES))

# -----------------------------
# Merge into DataFrame
# -----------------------------
df = pd.DataFrame({
    "residue": residue_numbers,
    "aa": aa_list,
    "secondary_structure": ss_list,
    "rsa": rsa_list,
})

dist_df = pd.DataFrame({
    "residue": [res.get_id()[1] for res in rbd_residues],
    "dist_to_ACE2": dist_to_ace2,
    "dist_to_RBM": dist_to_rbm,
})

df = pd.merge(df, dist_df, on="residue", how="inner")

# -----------------------------
# Add one-hot SS and normalized RSA
# -----------------------------
ss_types = ["H", "E", "T", "C"]  # Helix, Strand, Turn, Coil
for ss_type in ss_types:
    df[f"ss_{ss_type}"] = (df["secondary_structure"] == ss_type).astype(int)

# Normalize RSA between 0 and 1
df["rsa_norm"] = (df["rsa"] - df["rsa"].min()) / (df["rsa"].max() - df["rsa"].min())

# -----------------------------
# Save
# -----------------------------
os.makedirs(os.path.dirname(output_csv), exist_ok=True)
df.to_csv(output_csv, index=False)

print(f"✅ Context features saved to {output_csv}")
print(df.head())
