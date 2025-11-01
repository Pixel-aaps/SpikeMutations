# scripts/01_download_dms.py

import pandas as pd

raw_path = "data/raw/single_mut_effects.csv"
cleaned_path = "data/interim/dms_core.csv"

# Load raw CSV
df = pd.read_csv(raw_path)

# Keep key columns
df_clean = df[[
    "site_RBD",
    "site_SARS2",
    "wildtype",
    "mutant",
    "mutation",        
    "mutation_RBD",    
    "bind_lib1", "bind_lib2", "bind_avg",
    "expr_lib1", "expr_lib2", "expr_avg"
]].rename(columns={
    "bind_avg": "ace2_score",
    "expr_avg": "expr_score"
})

# Add an identifier column for uniqueness
df_clean["mut_id"] = (
    df_clean["wildtype"] + df_clean["site_SARS2"].astype(str) + df_clean["mutant"]
)

# Save
df_clean.to_csv(cleaned_path, index=False)
print(f"Saved cleaned DMS to {cleaned_path}")
