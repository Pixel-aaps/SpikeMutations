from rdkit import Chem

# --- Load the SDF file ---
sdf_path = "data/raw/drug_ligands.sdf"
supplier = Chem.SDMolSupplier(sdf_path)

# --- Extract and clean names ---
ligand_names = []

for mol in supplier:
    if mol is None:
        continue
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else "unknown"
    # Remove CID if present, lowercase it
    cleaned_name = name.split("(")[0].strip().lower()
    ligand_names.append(cleaned_name)

# --- Print result ---
print("Ligands found in SDF:")
for name in ligand_names:
    print(f"  - {name}")

print(f"\nTotal ligands found: {len(ligand_names)}")
