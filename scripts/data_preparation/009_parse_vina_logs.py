# robust_parse_vina_logs.py
import glob, os, re, pandas as pd

out_folder = "vina_wt_out"
log_paths = sorted(glob.glob(os.path.join(out_folder, "*.log")))
pdbqt_paths = sorted(glob.glob(os.path.join(out_folder, "*.pdbqt")))

def parse_pdbqt_for_scores(pdbqt_path):
    """Parse .pdbqt produced by vina for REMARK VINA RESULT lines."""
    best = None
    with open(pdbqt_path, "r", errors="ignore") as fh:
        for line in fh:
            if "REMARK" in line.upper() and "VINA" in line.upper() and "RESULT" in line.upper():
                m = re.search(r"(-?\d+\.\d+)", line)
                if m:
                    val = float(m.group(1))
                    if best is None or val < best:
                        best = val
    return best

def parse_log_for_scores(log_path):
    """Parse the vina log. Try several patterns and return best (most negative) affinity."""
    best = None
    with open(log_path, "r", errors="ignore") as fh:
        lines = fh.readlines()
    # Try lines explicitly containing VINA RESULT
    for line in lines:
        if "REMARK" in line.upper() and "VINA" in line.upper() and "RESULT" in line.upper():
            m = re.search(r"(-?\d+\.\d+)", line)
            if m:
                val = float(m.group(1))
                if best is None or val < best:
                    best = val
    # Try lines that look like "Affinity: -7.3 (kcal/mol)" or "Score: -7.3"
    if best is None:
        for line in lines:
            if "AFFINITY" in line.upper() or "SCORE" in line.upper() or "kcal/mol" in line.lower():
                m = re.search(r"(-?\d+\.\d+)", line)
                if m:
                    val = float(m.group(1))
                    if best is None or val < best:
                        best = val
    # As a final fallback, search any floats and pick the most negative
    if best is None:
        for line in lines:
            m = re.findall(r"(-?\d+\.\d+)", line)
            if m:
                vals = [float(x) for x in m]
                mn = min(vals)
                if best is None or mn < best:
                    best = mn
    return best, lines

rows = []
debug_report = []

for log in log_paths:
    basefile = os.path.basename(log)
    # e.g. "Apilimod_wt.log" or "Apilimod(10173277)_3D_h_wt.log"
    # remove common suffixes .log / _wt etc
    name = os.path.splitext(basefile)[0]
    name = re.sub(r"_?wt$","", name, flags=re.IGNORECASE)
    # try to find a matching pdbqt (out) file
    # common patterns: <name>_wt_out.pdbqt, <name>_out.pdbqt, <name>.out.pdbqt, <name>_out.pdbqt
    candidates = []
    patterns = [f"{name}*_out.pdbqt", f"{name}*.out.pdbqt", f"{name}*.pdbqt"]
    for p in patterns:
        candidates.extend(glob.glob(os.path.join(out_folder, p)))
    candidates = sorted(set(candidates))
    score = None

    # 1) prefer parsing pdbqt (contains REMARK VINA RESULT)
    for pdbqt in candidates:
        s = parse_pdbqt_for_scores(pdbqt)
        if s is not None:
            score = s
            break

    # 2) if not found, parse the log file robustly
    if score is None:
        s, lines = parse_log_for_scores(log)
        score = s
        # if score still None, add debug snippet for this log
        if score is None:
            snippet = "".join(lines[:30]) if lines else ""
            debug_report.append((name, log, len(lines), snippet[:800]))

    rows.append((name, score))

# Also try to catch orphan pdbqt files that had no .log counterpart
seen_names = {r[0] for r in rows}
for pdbqt in pdbqt_paths:
    base = os.path.splitext(os.path.basename(pdbqt))[0]
    # strip suffixes like _wt_out, _out, .out
    base_norm = re.sub(r"(_wt_out|_out|\.out|_out\.pdbqt)$", "", base, flags=re.IGNORECASE)
    if base_norm in seen_names:
        continue
    s = parse_pdbqt_for_scores(pdbqt)
    rows.append((base_norm, s))

df = pd.DataFrame(rows, columns=["LigandName","DockingScore_WT"]).sort_values("LigandName").reset_index(drop=True)
df.to_csv(os.path.join(out_folder,"wt_scores.csv"), index=False)

# Print summary
print(df)
if debug_report:
    print("\nSome logs had no parsable scores. Showing up to 5 debug snippets:")
    for i, (name, logpath, nlines, snippet) in enumerate(debug_report[:5]):
        print(f"\n--- {i+1}. {name} | file: {logpath} | lines: {nlines} ---\n{snippet}\n--- END SNIPPET ---\n")

# Also print any zero-size or tiny log/pdbqt files to help debugging
small_files = []
for p in glob.glob(os.path.join(out_folder,"*")):
    try:
        sz = os.path.getsize(p)
    except:
        sz = 0
    if sz < 200:  # 200 bytes cut-off for "probably empty"
        small_files.append((p, sz))
if small_files:
    print("Warning: small/empty files detected (likely failed runs):")
    for p,sz in small_files:
        print(p, sz)
