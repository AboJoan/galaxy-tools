#!/usr/bin/env python3
import re
import json
import argparse
from pathlib import Path
from collections import Counter, defaultdict


def extract_recorddata(js_path: Path):
    txt = js_path.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"var\s+recordData\s*=\s*(\[\{.*?\}\]);", txt, flags=re.S)
    if not m:
        return None
    try:
        return json.loads(m.group(1))
    except json.JSONDecodeError:
        return None


def classify(prod: str) -> str:
    p = (prod or "").lower()
    if "nrps" in p:
        return "NRPS"
    if "t1pks" in p or "t2pks" in p or "t3pks" in p or "pks" in p:
        return "PKS"
    if "ripp" in p:
        return "RiPP"
    if "terpene" in p:
        return "Terpene"
    if "phosphonate" in p:
        return "Phosphonate"
    if "autoinducer" in p or "auto-inducer" in p:
        return "Auto-inducer"
    if "deazapurine" in p:
        return "Deazapurine"
    return "Other"


def clean_genome_name(js_name: str) -> str:
    """
    Convert current Module 4 antiSMASH JS names into clean genome labels.

    Examples:
      01_reannotated_regions.js -> Re-annotated_01
      08_reannotated_regions.js -> Re-annotated_08
      regions.js                -> regions
    """
    name = Path(js_name).name

    m = re.match(r"^(.+?)_reannotated_regions\.js$", name)
    if m:
        return f"Re-annotated_{m.group(1)}"

    m = re.match(r"^(.+?)_regions\.js$", name)
    if m:
        return m.group(1)

    return name.replace(".js", "")


def load_name_map(path: str | None):
    if not path:
        return {}
    mp = Path(path)
    if not mp.exists():
        raise SystemExit(f"[ERROR] name_map not found: {mp}")
    mapping = {}
    for ln in mp.read_text(encoding="utf-8", errors="replace").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if len(parts) >= 2:
            mapping[parts[0]] = parts[1]
    return mapping


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_dir", required=True)
    ap.add_argument("--out", dest="out_tsv", required=True, help="Base output TSV long table.")
    ap.add_argument("--out_csv", dest="out_csv", default=None)
    ap.add_argument("--name_map", dest="name_map", default=None)
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_long_tsv = Path(args.out_tsv)
    out_long_csv = Path(args.out_csv) if args.out_csv else out_long_tsv.with_suffix(".csv")

    out_wide_tsv = out_long_tsv.with_name(out_long_tsv.stem.replace("counts", "counts_wide") + ".tsv")
    out_wide_csv = out_long_tsv.with_name(out_long_tsv.stem.replace("counts", "counts_wide") + ".csv")

    if not in_dir.exists():
        raise SystemExit(f"[ERROR] input dir not found: {in_dir}")

    out_long_tsv.parent.mkdir(parents=True, exist_ok=True)

    mapping = load_name_map(args.name_map)

    js_files = sorted(in_dir.glob("*.js"))
    if not js_files:
        raise SystemExit(f"[ERROR] no .js files found in: {in_dir}")

    counts = Counter()
    all_classes = set()
    all_genomes = set()

    for js_path in js_files:
        default_name = clean_genome_name(js_path.name)
        genome = mapping.get(js_path.name, mapping.get(default_name, default_name))
        all_genomes.add(genome)

        data = extract_recorddata(js_path)
        if not data:
            print(f"[WARN] could not parse recordData from: {js_path}")
            continue

        for rec in data:
            for region in rec.get("regions", []):
                products = region.get("products", []) or []
                if not products:
                    cls = "Other"
                    counts[(genome, cls)] += 1
                    all_classes.add(cls)
                    continue

                for prod in products:
                    cls = classify(prod)
                    counts[(genome, cls)] += 1
                    all_classes.add(cls)

    with out_long_tsv.open("w", encoding="utf-8") as f:
        f.write("Genomes\tCanonical_Class\tBGC_count\n")
        for (genome, cls), n in sorted(counts.items()):
            f.write(f"{genome}\t{cls}\t{n}\n")

    with out_long_csv.open("w", encoding="utf-8") as f:
        f.write("Genomes,Canonical_Class,BGC_count\n")
        for (genome, cls), n in sorted(counts.items()):
            f.write(f"{genome},{cls},{n}\n")

    classes_sorted = sorted(all_classes)
    matrix = defaultdict(lambda: defaultdict(int))
    totals = defaultdict(int)

    for (genome, cls), n in counts.items():
        matrix[genome][cls] = n
        totals[genome] += n

    with out_wide_tsv.open("w", encoding="utf-8") as f:
        header = ["Genomes"] + classes_sorted + ["Total_BGCs"]
        f.write("\t".join(header) + "\n")
        for genome in sorted(all_genomes):
            row = [genome] + [str(matrix[genome].get(cls, 0)) for cls in classes_sorted] + [str(totals.get(genome, 0))]
            f.write("\t".join(row) + "\n")

    with out_wide_csv.open("w", encoding="utf-8") as f:
        header = ["Genomes"] + classes_sorted + ["Total_BGCs"]
        f.write(",".join(header) + "\n")
        for genome in sorted(all_genomes):
            row = [genome] + [str(matrix[genome].get(cls, 0)) for cls in classes_sorted] + [str(totals.get(genome, 0))]
            f.write(",".join(row) + "\n")

    print("[OK] wrote:", out_long_tsv)
    print("[OK] wrote:", out_long_csv)
    print("[OK] wrote:", out_wide_tsv)
    print("[OK] wrote:", out_wide_csv)


if __name__ == "__main__":
    main()
