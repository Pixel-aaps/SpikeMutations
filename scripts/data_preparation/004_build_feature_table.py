# # """
# # Build structural context features for SARS-CoV-2 RBD from 6M0J.pdb

# # Outputs:
# #     data/processed/context_features.csv
# # """

# # import os
# # import pandas as pd
# # from Bio.PDB import PDBParser, DSSP, NeighborSearch
# # from Bio.PDB.Polypeptide import is_aa
# # from Bio.SeqUtils import seq1

# # # -----------------------------
# # # Config
# # # -----------------------------
# # pdb_file = "data/external/6M0J.pdb"
# # output_csv = "data/processed/context_features.csv"

# # # SARS-CoV-2 RBD is residues 319–541 in the spike protein
# # RBD_START, RBD_END = 319, 541
# # ACE2_CHAIN = "A"  # ACE2 receptor chain
# # RBM_RESIDUES = range(438, 507)  # Receptor binding motif residues

# # # -----------------------------
# # # Utility
# # # -----------------------------
# # def try_three_to_one(resname: str) -> str:
# #     """Convert 3-letter AA code to 1-letter, fallback to 'X'."""
# #     try:
# #         return seq1(resname)
# #     except Exception:
# #         return "X"

# # # -----------------------------
# # # Load structure
# # # -----------------------------
# # parser = PDBParser(QUIET=True)
# # structure = parser.get_structure("complex", pdb_file)
# # model = structure[0]

# # # Identify RBD chain automatically by residue range (319–541)
# # rbd_chain = None
# # for chain in model:
# #     res_ids = [res.get_id()[1] for res in chain if is_aa(res, standard=True)]
# #     if any(RBD_START <= r <= RBD_END for r in res_ids):
# #         rbd_chain = chain
# #         break

# # if rbd_chain is None:
# #     raise ValueError("Could not find RBD chain (319–541) in 6M0J.pdb")

# # print(f"RBD chain detected: {rbd_chain.id}")

# # # Get ACE2 chain
# # ace2_chain = model[ACE2_CHAIN]

# # # Filter RBD residues (standard AA only, 319–541)
# # rbd_residues = [res for res in rbd_chain if is_aa(res, standard=True)]
# # rbd_residues = [res for res in rbd_residues if RBD_START <= res.get_id()[1] <= RBD_END]

# # # -----------------------------
# # # DSSP: secondary structure + RSA
# # # -----------------------------
# # dssp = DSSP(model, pdb_file, dssp="mkdssp")

# # residue_numbers, aa_list, ss_list, rsa_list = [], [], [], []
# # for key in dssp.keys():
# #     chain_id, res_id = key[0], key[1][1]
# #     if chain_id != rbd_chain.id:
# #         continue
# #     if not (RBD_START <= res_id <= RBD_END):
# #         continue
# #     residue_numbers.append(res_id)
# #     aa_list.append(dssp[key][1])
# #     ss_list.append(dssp[key][2])
# #     rsa_list.append(dssp[key][3])

# # # -----------------------------
# # # Distance to ACE2
# # # -----------------------------
# # ace2_atoms = [atom for atom in ace2_chain.get_atoms()]
# # dist_to_ace2 = []
# # for res in rbd_residues:
# #     atoms = list(res.get_atoms())
# #     if len(atoms) == 0:
# #         dist_to_ace2.append(None)
# #         continue
# #     min_dist = min(atom - ace_atom for atom in atoms for ace_atom in ace2_atoms)
# #     dist_to_ace2.append(min_dist)

# # # -----------------------------
# # # Distance to RBM
# # # -----------------------------
# # dist_to_rbm = []
# # for res in rbd_residues:
# #     res_id = res.get_id()[1]
# #     if res_id in RBM_RESIDUES:
# #         dist_to_rbm.append(0)
# #     else:
# #         dist_to_rbm.append(min(abs(res_id - r) for r in RBM_RESIDUES))

# # # -----------------------------
# # # Merge into DataFrame
# # # -----------------------------
# # df = pd.DataFrame({
# #     "residue": residue_numbers,
# #     "aa": aa_list,
# #     "secondary_structure": ss_list,
# #     "rsa": rsa_list,
# # })

# # # Align DSSP and distance lists by residue ID
# # dist_df = pd.DataFrame({
# #     "residue": [res.get_id()[1] for res in rbd_residues],
# #     "dist_to_ACE2": dist_to_ace2,
# #     "dist_to_RBM": dist_to_rbm,
# # })

# # df = pd.merge(df, dist_df, on="residue", how="inner")

# # # -----------------------------
# # # Save
# # # -----------------------------
# # os.makedirs(os.path.dirname(output_csv), exist_ok=True)
# # df.to_csv(output_csv, index=False)

# # print(f"Context features saved to {output_csv}")
# # print(df.head())


# """
# Build structural context features for SARS-CoV-2 RBD from 6M0J.pdb

# Outputs:
#     data/processed/context_features.csv
# """

# import os
# import pandas as pd
# from Bio.PDB import PDBParser, DSSP, NeighborSearch
# from Bio.PDB.Polypeptide import is_aa
# from Bio.SeqUtils import seq1

# # -----------------------------
# # Config
# # -----------------------------
# pdb_file = "data/external/6m0j_Repair(1).pdb"
# output_csv = "data/processed/context_features.csv"

# # SARS-CoV-2 RBD is residues 319–541 in the spike protein
# RBD_START, RBD_END = 319, 541
# ACE2_CHAIN = "A"  # ACE2 receptor chain
# RBM_RESIDUES = range(438, 507)  # Receptor binding motif residues

# # -----------------------------
# # Utility
# # -----------------------------
# def try_three_to_one(resname: str) -> str:
#     """Convert 3-letter AA code to 1-letter, fallback to 'X'."""
#     try:
#         return seq1(resname)
#     except Exception:
#         return "X"

# # -----------------------------
# # Load structure
# # -----------------------------
# parser = PDBParser(QUIET=True)
# structure = parser.get_structure("complex", pdb_file)
# model = structure[0]

# # Identify RBD chain automatically by residue range (319–541)
# rbd_chain = None
# for chain in model:
#     res_ids = [res.get_id()[1] for res in chain if is_aa(res, standard=True)]
#     if any(RBD_START <= r <= RBD_END for r in res_ids):
#         rbd_chain = chain
#         break

# if rbd_chain is None:
#     raise ValueError("Could not find RBD chain (319–541) in 6M0J.pdb")

# print(f"RBD chain detected: {rbd_chain.id}")

# # Get ACE2 chain
# ace2_chain = model[ACE2_CHAIN]

# # Filter RBD residues (standard AA only, 319–541)
# rbd_residues = [res for res in rbd_chain if is_aa(res, standard=True)]
# rbd_residues = [res for res in rbd_residues if RBD_START <= res.get_id()[1] <= RBD_END]

# # -----------------------------
# # DSSP: secondary structure + RSA
# # -----------------------------
# dssp = DSSP(model, pdb_file, dssp="mkdssp")

# residue_numbers, aa_list, ss_list, rsa_list = [], [], [], []
# for key in dssp.keys():
#     chain_id, res_id = key[0], key[1][1]
#     if chain_id != rbd_chain.id:
#         continue
#     if not (RBD_START <= res_id <= RBD_END):
#         continue
#     residue_numbers.append(res_id)
#     aa_list.append(dssp[key][1])
#     ss_list.append(dssp[key][2])
#     rsa_list.append(dssp[key][3])

# # -----------------------------
# # Distance to ACE2
# # -----------------------------
# ace2_atoms = [atom for atom in ace2_chain.get_atoms()]
# dist_to_ace2 = []
# for res in rbd_residues:
#     atoms = list(res.get_atoms())
#     if len(atoms) == 0:
#         dist_to_ace2.append(None)
#         continue
#     min_dist = min(atom - ace_atom for atom in atoms for ace_atom in ace2_atoms)
#     dist_to_ace2.append(min_dist)

# # -----------------------------
# # Distance to RBM
# # -----------------------------
# dist_to_rbm = []
# for res in rbd_residues:
#     res_id = res.get_id()[1]
#     if res_id in RBM_RESIDUES:
#         dist_to_rbm.append(0)
#     else:
#         dist_to_rbm.append(min(abs(res_id - r) for r in RBM_RESIDUES))

# # -----------------------------
# # Merge into DataFrame
# # -----------------------------
# df = pd.DataFrame({
#     "residue": residue_numbers,
#     "aa": aa_list,
#     "secondary_structure": ss_list,
#     "rsa": rsa_list,
# })

# # Align DSSP and distance lists by residue ID
# dist_df = pd.DataFrame({
#     "residue": [res.get_id()[1] for res in rbd_residues],
#     "dist_to_ACE2": dist_to_ace2,
#     "dist_to_RBM": dist_to_rbm,
# })

# df = pd.merge(df, dist_df, on="residue", how="inner")

# # -----------------------------
# # Save
# # -----------------------------
# os.makedirs(os.path.dirname(output_csv), exist_ok=True)
# df.to_csv(output_csv, index=False)

# print(f"Context features saved to {output_csv}")
# print(df.head())

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
output_csv = "data/processed/context_features1.csv"

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
# Save
# -----------------------------
os.makedirs(os.path.dirname(output_csv), exist_ok=True)
df.to_csv(output_csv, index=False)

print(f"✅ Context features saved to {output_csv}")
print(df.head())
