import pandas as pd
import os

# ==== 1. Load BLOSUM62 Matrix ====
def load_blosum62(file):
    blosum = {}
    with open(file) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    headers = lines[0].split()
    for line in lines[1:]:
        parts = line.split()
        row_aa = parts[0]
        scores = list(map(int, parts[1:]))
        for j, col_aa in enumerate(headers):
            blosum[(row_aa, col_aa)] = scores[j]
    return blosum

# ==== 2. Load Grantham Distances (auto-detect format) ====
def load_grantham(file):
    df = pd.read_csv(file)
    if {"wt", "mut", "distance"}.issubset(df.columns):  
        # long format
        return {(row["wt"], row["mut"]): row["distance"] for _, row in df.iterrows()}
    else:
        # matrix format
        df = pd.read_csv(file, index_col=0)
        return {(i, j): df.loc[i, j] for i in df.index for j in df.columns}

# ==== 3. Physicochemical properties ====
charge = {'A': 0,'R': 1,'N': 0,'D': -1,'C': 0,'E': -1,'Q': 0,'G': 0,
          'H': 1,'I': 0,'L': 0,'K': 1,'M': 0,'F': 0,'P': 0,'S': 0,
          'T': 0,'W': 0,'Y': 0,'V': 0}

hydropathy = {'A': 1.8,'R': -4.5,'N': -3.5,'D': -3.5,'C': 2.5,'E': -3.5,'Q': -3.5,'G': -0.4,
              'H': -3.2,'I': 4.5,'L': 3.8,'K': -3.9,'M': 1.9,'F': 2.8,'P': -1.6,'S': -0.8,
              'T': -0.7,'W': -0.9,'Y': -1.3,'V': 4.2}

volume = {'A':  88,'R': 173,'N': 114,'D': 111,'C': 108,'E': 138,'Q': 143,'G': 60,
          'H': 153,'I': 166,'L': 166,'K': 168,'M': 162,'F': 189,'P': 112,'S': 89,
          'T': 116,'W': 227,'Y': 193,'V': 140}

polarity = {'A':  8.1,'R': 10.5,'N': 11.6,'D': 13.0,'C': 5.5,'E': 12.3,'Q': 10.5,'G': 9.0,
            'H': 10.4,'I': 5.2,'L': 4.9,'K': 11.3,'M': 5.7,'F': 5.2,'P': 8.0,'S': 9.2,
            'T': 8.6,'W': 5.4,'Y': 6.2,'V': 5.9}

# ==== 4. Load input files ====
blosum62 = load_blosum62("data/external/blosum62.txt")
grantham = load_grantham("data/external/grantham_distances.csv")
labels = pd.read_csv("data/processed/labels.csv")

# ==== 5. Compute features ====
rows = []
for _, row in labels.iterrows():
    wt, mut = row["wildtype"], row["mutant"]

    bscore = blosum62.get((wt, mut), 0)
    gscore = grantham.get((wt, mut), 0)

    is_stop = int(mut == "*" or wt == "*")

    if is_stop:
        delta_charge = 0
        delta_hydro = 0
        delta_volume = 0
        delta_polar = 0
    else:
        delta_charge = charge.get(mut, 0) - charge.get(wt, 0)
        delta_hydro  = hydropathy.get(mut, 0) - hydropathy.get(wt, 0)
        delta_volume = volume.get(mut, 0) - volume.get(wt, 0)
        delta_polar  = polarity.get(mut, 0) - polarity.get(wt, 0)

    rows.append({
        "mut_id": row["mut_id"],
        "pos": row["site_RBD"],
        "wt": wt, "mut": mut,
        "blosum62": bscore,
        "grantham": gscore,
        "delta_charge": delta_charge,
        "delta_hydro": delta_hydro,
        "delta_volume": delta_volume,
        "delta_polarity": delta_polar,
        "is_stop": is_stop
    })



seq_features = pd.DataFrame(rows)


# ==== 6. Save ====
os.makedirs("data/processed", exist_ok=True)
seq_features.to_csv("data/processed/seq_features.csv", index=False)

print("Saved data/processed/seq_features.csv")
