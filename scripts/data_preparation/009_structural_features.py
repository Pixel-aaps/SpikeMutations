# structure_features.py
import os
import pandas as pd
from Bio.PDB import PDBParser, DSSP, NeighborSearch, Selection
import numpy as np

def compute_dssp_features(pdb_file, chain_id, output_csv="structure_features.csv"):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = structure[0]
    dssp = DSSP(model, pdb_file, acc_array='Sander')  # requires mkdssp
    rows = []
    for key in dssp.keys():
        (chain, res_id) = key
        if chain != chain_id:
            continue
        aa = dssp[key][1]
        ss = dssp[key][2]
        asa = dssp[key][3]  # absolute ASA
        rows.append({
            "chain": chain,
            "resid": res_id[1],
            "aa": aa,
            "sec_struct": ss,
            "asa": asa
        })
    df = pd.DataFrame(rows)

    # relative ASA
    max_acc = {
      'A': 129,'R': 274,'N': 195,'D': 193,'C': 167,'Q': 225,'E': 223,'G': 104,'H': 224,'I': 197,
      'L': 201,'K': 236,'M': 224,'F': 240,'P': 159,'S': 155,'T': 172,'W': 285,'Y': 263,'V': 174
    }
    df['max_acc'] = df['aa'].map(max_acc)
    df['rsa'] = df['asa'] / df['max_acc']
    df.to_csv(output_csv, index=False)
    return df

def neighbor_counts(pdb_file, chain_id, radius=8.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)
    model = structure[0]
    atoms = Selection.unfold_entities(model, "A")
    ns = NeighborSearch(atoms)
    rows = []
    for chain in model:
        if chain.id != chain_id: continue
        for res in chain:
            if res.id[0] != " ": continue
            ca = res.get('CA')
            if ca is None: 
                rows.append({"resid": res.id[1], "neighbors": None})
                continue
            neighbors = ns.search(ca.coord, radius, level="R")
            rows.append({"resid": res.id[1], "neighbors": len(neighbors)-1})
    dfn = pd.DataFrame(rows)
    return dfn

# ---------- Run Section ----------
if __name__ == "__main__":
    pdb_path = "data/external/fold_x/6m0j_Repair.pdb"   # <-- Linux style
    df_dssp = compute_dssp_features(pdb_path, chain_id="E", output_csv="structure_features.csv")
    df_neigh = neighbor_counts(pdb_path, chain_id="E")

    merged = df_dssp.merge(df_neigh, on="resid")
    merged.to_csv("merged_structure_features.csv", index=False)
