#!/usr/bin/env bash
set -euo pipefail

# =========================
# User-configurable paths
# =========================
: "${EGG_DIR:=00_input/01_eggNOG_functional_annotation/01_bacterial_protein_genes}"   # where Re-annotated_XX.csv +Ref_NCBI_mapped_XX.csv live
: "${FAA_DIR:=00_input/01_eggNOG_functional_annotation/01_bacterial_protein_genes}"  # where Re-annotated_XX.faa +Ref_NCBI_mapped_XX.faa live (set properly if different)
: "${OBO_FILE:=00_input/go-basic.obo}"

: "${OUT_TSV:=04_results/01_final_tables/01_Bacterial_genome_functional_summary_table.tsv}"
: "${OUT_CSV:=04_results/01_final_tables/01_Bacterial_genome_functional_summary_table.csv}"

: "${NEW_LABEL:=Re-annotation Genomes (prodigal-only)}"
: "${OLD_LABEL:=Reference Genomes (NCBI-only)}"

# Optional species mapping file (TSV): index<TAB>species
# Example:
# 01    Blautia producta
# 02    Bacteroides thetaiotaomicron
: "${SPECIES_TSV:=}"

export EGG_DIR FAA_DIR OBO_FILE OUT_TSV OUT_CSV NEW_LABEL OLD_LABEL SPECIES_TSV

python3 - <<'PY'
import os, re, sys
from pathlib import Path
import pandas as pd

EGG_DIR = Path(os.environ.get("EGG_DIR", "00_input/01_eggNOG_functional_annotation/01_bacterial_protein_genes"))
FAA_DIR = Path(os.environ.get("FAA_DIR", str(EGG_DIR)))
OBO_FILE = Path(os.environ.get("OBO_FILE", "00_input/go-basic.obo"))

OUT_TSV = Path(os.environ.get("OUT_TSV", "04_results/01_final_tables/01_Bacterial_genome_functional_summary_table.tsv"))
OUT_CSV = Path(os.environ.get("OUT_CSV", "04_results/01_final_tables/01_Bacterial_genome_functional_summary_table.csv"))

NEW_LABEL = os.environ.get("NEW_LABEL", "Re-annotation Genomes (prodigal-only)")
OLD_LABEL = os.environ.get("OLD_LABEL", "Reference Genomes (NCBI-only)")

SPECIES_TSV = os.environ.get("SPECIES_TSV", "").strip()

def die(msg: str):
    raise SystemExit(msg)

# ----------------------------
# Optional species mapping
# ----------------------------
def load_species_map(path: str) -> dict[str, str]:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        die(f"[ERROR] SPECIES_TSV provided but not found: {p}")
    m = {}
    with p.open("r", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            idx = parts[0].strip()
            name = parts[1].strip()
            if idx and name:
                m[idx] = name
    return m

species_map = load_species_map(SPECIES_TSV)

def label_for(idx: str) -> str:
    return species_map.get(idx, f"IDX_{idx}")

# ----------------------------
# Helpers: files + parsing
# ----------------------------
def find_eggnog_file(prefix: str) -> Path:
    p = EGG_DIR / f"{prefix}.csv"
    if p.exists():
        return p
    die(f"[ERROR] Cannot find eggNOG CSV: {p}")

def find_faa_file(prefix: str) -> Path:
    p = FAA_DIR / f"{prefix}.faa"
    if p.exists():
        return p
    die(f"[ERROR] Cannot find FAA: {p}")

def count_faa_headers(faa_path: Path) -> int:
    n = 0
    with faa_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for ln in fh:
            if ln.startswith(">"):
                n += 1
    return n

def sniff_sep(path: Path) -> str:
    # eggNOG tables may be TAB-separated even if .csv
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if not ln.strip():
                continue
            if "\t" in ln:
                return "\t"
            if "," in ln:
                return ","
            break
    return "\t"

def load_eggnog(path: Path) -> pd.DataFrame:
    sep = sniff_sep(path)
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    return df

def must_have_cols(df: pd.DataFrame, cols: list[str], file_path: Path):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        die(f"[ERROR] Missing columns {missing} in {file_path}. Found: {list(df.columns)}")

# ----------------------------
# GO namespace map from OBO
# ----------------------------
def load_go_namespace_map(obo_path: Path) -> dict[str, str]:
    if not obo_path.exists():
        die(f"[ERROR] GO OBO file not found: {obo_path}")

    ns_map_full_to_short = {
        "biological_process": "BP",
        "molecular_function": "MF",
        "cellular_component": "CC",
    }

    go2ns = {}
    current_id = None
    current_ns = None

    with obo_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                current_id = None
                current_ns = None
                continue
            if not line or line.startswith("!"):
                continue
            if line.startswith("id: GO:"):
                current_id = line.split("id: ", 1)[1].strip()
                continue
            if line.startswith("namespace:"):
                current_ns = line.split("namespace:", 1)[1].strip()
                continue
            if current_id and current_ns:
                short = ns_map_full_to_short.get(current_ns)
                if short:
                    go2ns[current_id] = short
                current_id = None
                current_ns = None
    return go2ns

GO2NS = load_go_namespace_map(OBO_FILE)

GO_PATTERN = re.compile(r"GO:\d+")
KO_PATTERN = re.compile(r"(?:^|[,\s])(?:ko:)?K\d{5}(?=$|[,\s])")

def split_tokens(cell: str) -> list[str]:
    cell = (cell or "").strip()
    if not cell or cell == "-":
        return []
    cell = re.sub(r"\s+", " ", cell)
    return [t for t in re.split(r"[,;| ]+", cell) if t]

def letters_only(x: str) -> str:
    return re.sub(r"[^A-Za-z]", "", (x or "").upper())

def compute_metrics(eggnog_path: Path, faa_path: Path) -> dict[str, int]:
    df = load_eggnog(eggnog_path)
    must_have_cols(df, ["#query", "COG_category", "GOs", "KEGG_ko"], eggnog_path)

    total = count_faa_headers(faa_path)

    # ---------- COG ----------
    cog = df["COG_category"].astype(str).str.strip()
    no_hit_mask = (cog == "") | (cog == "-")
    cog_letters = cog.map(letters_only)

    unknown_mask = (~no_hit_mask) & (cog_letters != "") & (cog_letters.str.replace("S", "", regex=False) == "")
    known_mask   = (~no_hit_mask) & (cog_letters.str.replace("S", "", regex=False) != "")

    cog_nohit   = int(no_hit_mask.sum())
    cog_unknown = int(unknown_mask.sum())
    cog_known   = int(known_mask.sum())

    present = cog_known + cog_unknown + cog_nohit
    missing = total - present
    if missing < 0:
        missing = 0

    # ---------- KEGG KO ----------
    ko_col = df["KEGG_ko"].astype(str).str.strip()
    has_any_ko = (ko_col != "") & (ko_col != "-") & ko_col.str.contains(KO_PATTERN)
    kegg_annotated = int(has_any_ko.sum())

    kos = ko_col[has_any_ko].map(split_tokens)
    flat = []
    for toks in kos:
        for t in toks:
            t = re.sub(r"^ko:", "", t)
            if re.fullmatch(r"K\d{5}", t):
                flat.append(t)
    unique_ko = len(set(flat))

    # ---------- GO BP/MF/CC (unique genes per namespace) ----------
    go_col = df["GOs"].astype(str).str.strip()
    gene_col = df["#query"].astype(str).str.strip()

    bp_genes, mf_genes, cc_genes = set(), set(), set()

    for gene, cell in zip(gene_col, go_col):
        if not gene or gene == "-":
            continue
        if not cell or cell == "-":
            continue
        ids = set(GO_PATTERN.findall(cell))
        if not ids:
            continue
        namespaces = {GO2NS.get(go_id) for go_id in ids}
        if "BP" in namespaces:
            bp_genes.add(gene)
        if "MF" in namespaces:
            mf_genes.add(gene)
        if "CC" in namespaces:
            cc_genes.add(gene)

    return {
        "Total Proteins": total,
        "COG_Known Function": cog_known,
        "COG_Unknown Function": cog_unknown,
        "No Hit": cog_nohit,
        "Missing": missing,
        "KEGG Annotated": kegg_annotated,
        "Unique KO": unique_ko,
        "GO_BP": len(bp_genes),
        "GO_MF": len(mf_genes),
        "GO_CC": len(cc_genes),
    }

def pct_change(new: int, old: int) -> str:
    if old == 0:
        return "NA"
    return f"{((new - old) / old) * 100:.2f}"

# ----------------------------
# Detect available indices (paired CSVs + paired FAAs)
# ----------------------------
def detect_indices() -> list[str]:
    new_csv = set()
    old_csv = set()
    for p in EGG_DIR.glob("Re-annotated_*.csv"):
        m = re.search(r"Re-annotated_(\d+)\.csv$", p.name)
        if m:
            new_csv.add(m.group(1))
    for p in EGG_DIR.glob("Ref_NCBI_mapped_*.csv"):
        m = re.search(r"Ref_NCBI_mapped_(\d+)\.csv$", p.name)
        if m:
            old_csv.add(m.group(1))

    new_faa = set()
    old_faa = set()
    for p in FAA_DIR.glob("Re-annotated_*.faa"):
        m = re.search(r"Re-annotated_(\d+)\.faa$", p.name)
        if m:
            new_faa.add(m.group(1))
    for p in FAA_DIR.glob("Ref_NCBI_mapped_*.faa"):
        m = re.search(r"Ref_NCBI_mapped_(\d+)\.faa$", p.name)
        if m:
            old_faa.add(m.group(1))

    common = sorted((new_csv & old_csv) & (new_faa & old_faa), key=lambda x: int(x))
    if not common:
        die(
            f"[ERROR] No paired inputs detected.\n"
            f"  EGG_DIR={EGG_DIR}\n"
            f"  FAA_DIR={FAA_DIR}\n"
            f"Need pairs like:\n"
            f"  Re-annotated_01.csv +Ref_NCBI_mapped_01.csv\n"
            f"  Re-annotated_01.faa +Ref_NCBI_mapped_01.faa\n"
        )
    return common

indices = detect_indices()

metrics_order = [
    "Total Proteins",
    "COG_Known Function",
    "COG_Unknown Function",
    "No Hit",
    "Missing",
    "KEGG Annotated",
    "Unique KO",
    "GO_BP",
    "GO_MF",
    "GO_CC",
]

rows = []

for idx in indices:
    new_id = f"Re-annotated_{idx}"
    old_id = f"Ref_NCBI_mapped_{idx}"

    species = label_for(idx)

    new_egg = find_eggnog_file(new_id)
    old_egg = find_eggnog_file(old_id)

    new_faa = find_faa_file(new_id)
    old_faa = find_faa_file(old_id)

    new_m = compute_metrics(new_egg, new_faa)
    old_m = compute_metrics(old_egg, old_faa)

    for m in metrics_order:
        new_v = int(new_m.get(m, 0))
        old_v = int(old_m.get(m, 0))
        diff = new_v - old_v
        pc = pct_change(new_v, old_v)

        rows.append({
            "Index": idx,
            "Species": species,   # optional label (IDX_01 if no mapping)
            "Metric": m,
            NEW_LABEL: new_v,
            OLD_LABEL: old_v,
            "Δ (Difference)": diff,
            "% Change": pc
        })

out = pd.DataFrame(rows)

OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
out.to_csv(OUT_TSV, sep="\t", index=False)
out.to_csv(OUT_CSV, index=False)

print("[OK] Wrote:")
print(" -", OUT_TSV)
print(" -", OUT_CSV)
print("[INFO] Indices used:", ",".join(indices))
if SPECIES_TSV:
    print("[INFO] Species mapping:", SPECIES_TSV)
else:
    print("[INFO] No species mapping provided (Species=IDX_<nn>)")
PY

echo "[DONE] Table created:"
echo "  $OUT_TSV"
echo "  $OUT_CSV"
