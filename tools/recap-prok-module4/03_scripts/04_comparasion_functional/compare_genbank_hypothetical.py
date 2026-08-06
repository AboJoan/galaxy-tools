#!/usr/bin/env python3
import argparse
import re
import csv
from pathlib import Path
from typing import Dict, Tuple, Optional, List

CDS_START = re.compile(r"^\s{5}CDS\s")
QUAL_PROD = re.compile(r'^\s{21}/product=')

# ----------------------------
# Core GenBank parsing
# ----------------------------
def parse_product_value(line: str) -> str:
    s = line.strip()
    if "/product=" not in s:
        return ""
    val = s.split("/product=", 1)[1].strip()
    if val.startswith('"') and val.endswith('"'):
        val = val[1:-1]
    return val


def count_metrics_gb(gb_path: Path) -> Dict[str, int]:
    cds_total = 0
    hypo = 0

    in_cds = False
    saw_product = False
    product_val = ""

    def finish_cds():
        nonlocal in_cds, saw_product, product_val, hypo
        if not in_cds:
            return
        if not saw_product:
            hypo += 1
        else:
            if "hypothetical protein" in product_val.lower():
                hypo += 1
        in_cds = False
        saw_product = False
        product_val = ""

    with gb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if CDS_START.match(line):
                finish_cds()
                in_cds = True
                cds_total += 1
                continue

            if in_cds:
                if line.startswith("     ") and not CDS_START.match(line) and not line.startswith("                     "):
                    finish_cds()
                    continue
                if QUAL_PROD.match(line):
                    saw_product = True
                    product_val = parse_product_value(line)

    finish_cds()
    functional = max(0, cds_total - hypo)
    return {"cds_total": cds_total, "hypo": hypo, "functional": functional}


def fmt_n_pct(n: int, total: int) -> str:
    if total <= 0:
        return f"{n} (0.00%)"
    pct = (n / total) * 100.0
    return f"{n} ({pct:.2f}%)"


def organism_from_gb(gb_path: Path) -> str:
    """
    Extract ORGANISM name if present.
    If missing, returns gb_path.stem (fallback).
    Note: Some GenBanks may have ORGANISM '.' which is not useful.
    """
    org_re = re.compile(r"^ {2}ORGANISM\s+(.+?)\s*$")
    with gb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = org_re.match(line)
            if m:
                return m.group(1).strip()
            if line.startswith("FEATURES"):
                break
    return gb_path.stem

# ----------------------------
# File listing / pairing
# ----------------------------
def list_gb_files(dirp: Path) -> List[Path]:
    exts = (".gb", ".gbk", ".gbff", ".genbank")
    files: List[Path] = []
    for ext in exts:
        files.extend(dirp.glob(f"*{ext}"))
    return sorted(set(files))


def parse_key(stem: str, mode: str) -> Optional[str]:
    """
    Returns pairing key:
      bac:  "01", "3" ... normalized to 2-digit "01","03"
      pls:  normal plasmid -> "02"
            putative -> "03_01"
    """
    if mode == "bac":
        m = re.match(r"^Re-annotated_(\d+)$", stem)
        if not m:
            return None
        idx = int(m.group(1))
        return f"{idx:02d}"

    # plasmids
    m = re.match(r"^PLRe-annotated_(\d+)$", stem)
    if m:
        idx = int(m.group(1))
        return f"{idx:02d}"

    m2 = re.match(r"^Putative_PLRe-annotated_(\d+)_(\d+)$", stem)
    if m2:
        idx = int(m2.group(1))
        sub = int(m2.group(2))
        return f"{idx:02d}_{sub:02d}"

    return None


def ref_key(stem: str, mode: str) -> Optional[str]:
    if mode == "bac":
        m = re.match(r"^Ref_NCBI_mapped_(\d+)$", stem)
        if not m:
            return None
        idx = int(m.group(1))
        return f"{idx:02d}"

    m = re.match(r"^PLRef_NCBI_mapped_(\d+)$", stem)
    if m:
        idx = int(m.group(1))
        return f"{idx:02d}"
  
    return None


def auto_pairs(reann_dir: Path, ref_dir: Path, mode: str) -> Dict[str, Tuple[Path, Path]]:
    reann_files = list_gb_files(reann_dir)
    ref_files = list_gb_files(ref_dir)

    ref_map: Dict[str, Path] = {}
    for f in ref_files:
        k = ref_key(f.stem, mode)
        if k:
            ref_map[k] = f

    pairs: Dict[str, Tuple[Path, Path]] = {}
    for f in reann_files:
        k = parse_key(f.stem, mode)
        if not k:
            continue
        if k in ref_map:
            pairs[k] = (f, ref_map[k])

    return pairs

# ----------------------------
# NEW: Fill Species from metadata Host via replicon accession in GenBank
# ----------------------------
def load_metadata_map(tsv_path: Path) -> Dict[str, str]:
    """
    Map replicon accessions (and their no-version forms) -> Host
    We will use:
      Replicon_accession
      Replicon_accession_norm
    If Host is empty or '.', skip row.
    """
    m: Dict[str, str] = {}

    with tsv_path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            host = (r.get("Host") or "").strip()
            if not host or host == ".":
                continue

            acc = (r.get("Replicon_accession") or "").strip()
            accn = (r.get("Replicon_accession_norm") or "").strip()

            if acc:
                m[acc] = host
                m[re.sub(r"\.\d+$", "", acc)] = host

            if accn:
                m[accn] = host
                m[re.sub(r"\.\d+$", "", accn)] = host

    return m


def replicon_from_genbank(gb_path: Path) -> str:
    """
    Extract replicon accession from GenBank header:
    Prefer ACCESSION; fallback to VERSION.
    """
    acc = ""
    with gb_path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ACCESSION"):
                parts = line.split()
                if len(parts) >= 2:
                    acc = parts[1].strip()
                    break
            if line.startswith("VERSION") and not acc:
                parts = line.split()
                if len(parts) >= 2:
                    acc = parts[1].strip()
            if line.startswith("FEATURES"):
                break
    return acc.strip()


def choose_species_label(re_gb: Path, meta_map: Dict[str, str]) -> str:
    """
    Priority:
      1) ORGANISM from GenBank (if not '.' and not empty)
      2) Host from metadata via replicon accession in GenBank header
      3) file stem (safe fallback)
    """
    org = organism_from_gb(re_gb).strip()
    if org and org != ".":
        return org

    acc = replicon_from_genbank(re_gb)
    acc0 = re.sub(r"\.\d+$", "", acc) if acc else ""

    if acc and acc in meta_map:
        return meta_map[acc]
    if acc0 and acc0 in meta_map:
        return meta_map[acc0]

    return re_gb.stem

# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Compare hypothetical proteins between Re-annotated vs Reference GenBanks."
    )
    ap.add_argument("--mode", choices=["bac", "pls"], required=True)
    ap.add_argument("--reann_dir", required=True)
    ap.add_argument("--ref_dir", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--label_reann", default="Re-annotated")
    ap.add_argument("--label_ref", default="NCBI Reference genome")

    # optional: metadata to fill Species when ORGANISM is missing
    ap.add_argument(
        "--metadata_tsv",
        default="00_input/03_genbank/01_metadata/03_replicon_metadata.tsv",
        help="Metadata TSV with Host + Replicon_accession(_norm) to fill Species when ORGANISM is missing."
    )

    args = ap.parse_args()

    reann_dir = Path(args.reann_dir)
    ref_dir = Path(args.ref_dir)
    out_prefix = Path(args.out_prefix)

    pairs = auto_pairs(reann_dir, ref_dir, args.mode)
    if not pairs:
        raise SystemExit(
            f"[ERROR] No pairs found in mode={args.mode}.\n"
            f"  reann_dir={reann_dir}\n  ref_dir={ref_dir}\n"
            f"  Expected examples:\n"
            f"    bac: Re-annotated_01.gb + Ref_NCBI_mapped_01.gb\n"
            f"    pls: PLRe-annotated_02.gb + PLRef_NCBI_mapped_02.gb\n"
            f"         PLRe-annotated_031.gb + PLRef_NCBI_mapped_031.gb\n"
        )

    # Load metadata map once (safe even if file exists but not used)
    meta_map: Dict[str, str] = {}
    meta_path = Path(args.metadata_tsv)
    if meta_path.exists():
        meta_map = load_metadata_map(meta_path)

    metrics_order = [
        "Protein coding CDSs",
        "Hypothetical proteins",
        "Functionally annotated proteins",
    ]

    rows = []
    for key in sorted(pairs.keys()):
        re_gb, rf_gb = pairs[key]

        re_m = count_metrics_gb(re_gb)
        rf_m = count_metrics_gb(rf_gb)

        # >>> ONLY CHANGE: choose Species label robustly
        species = choose_species_label(re_gb, meta_map)

        for metric in metrics_order:
            if metric == "Protein coding CDSs":
                re_val = str(re_m["cds_total"])
                rf_val = str(rf_m["cds_total"])
            elif metric == "Hypothetical proteins":
                re_val = fmt_n_pct(re_m["hypo"], re_m["cds_total"])
                rf_val = fmt_n_pct(rf_m["hypo"], rf_m["cds_total"])
            else:
                re_val = fmt_n_pct(re_m["functional"], re_m["cds_total"])
                rf_val = fmt_n_pct(rf_m["functional"], rf_m["cds_total"])

            rows.append([species, metric, re_val, rf_val])

    header = ["Species", "Metric", args.label_reann, args.label_ref]

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    csv_path = out_prefix.with_suffix(".csv")
    tsv_path = out_prefix.with_suffix(".tsv")

    with csv_path.open("w", newline="", encoding="utf-8") as out:
        w = csv.writer(out)
        w.writerow(header)
        w.writerows(rows)

    with tsv_path.open("w", newline="", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    print(f"[OK] Wrote: {csv_path}")
    print(f"[OK] Wrote: {tsv_path}")
    print(f"[OK] Pairs compared: {len(pairs)}  mode={args.mode}")


if __name__ == "__main__":
    main()
