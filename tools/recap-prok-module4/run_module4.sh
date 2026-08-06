#!/usr/bin/env bash
set -euo pipefail

INPUT_ARCHIVE="$1"
THREADS="$2"
TOOL_DIR="$3"

EGGNOG_ENV="/home/ataa/miniconda3/envs/eggnog"
MODULE4_ENV="/home/ataa/miniconda3/envs/module4"
ANTISMASH_ENV="/home/ataa/miniconda3/envs/antismash"
EGGNOG_DB="/home/ataa/eggnog_db"

RUN_DIR="$(pwd)"
WORKDIR="$RUN_DIR/module4_work_01"

rm -rf "$WORKDIR"
mkdir -p "$WORKDIR"

echo "[INFO] Checking input archive..."
file "$INPUT_ARCHIVE" || true

if tar -tzf "$INPUT_ARCHIVE" >/dev/null 2>&1; then
    echo "[INFO] Detected gzip-compressed tar archive"
    tar -xzf "$INPUT_ARCHIVE" -C "$WORKDIR"
elif tar -tf "$INPUT_ARCHIVE" >/dev/null 2>&1; then
    echo "[INFO] Detected uncompressed tar archive"
    tar -xf "$INPUT_ARCHIVE" -C "$WORKDIR"
else
    echo "[FATAL] Input is neither tar.gz nor tar"
    file "$INPUT_ARCHIVE" || true
    exit 1
fi

cd "$WORKDIR"

ln -sfn "$TOOL_DIR/03_scripts" 03_scripts
cp "$TOOL_DIR/Snakefile.master" .
cp "$TOOL_DIR/Snakefile.hypothetical" .


mkdir -p 04_results/00_logs
mkdir -p 00_input/03_genbank/03_eggNOG/01_bacteria
mkdir -p 00_input/03_genbank/03_eggNOG/02_plasmid

export EGGNOG_DATA_DIR="$EGGNOG_DB"

convert_csv () {
python3 - "$1" "$2" <<'PY'
import sys, csv
from pathlib import Path

ann = Path(sys.argv[1])
out = Path(sys.argv[2])
out.parent.mkdir(parents=True, exist_ok=True)

header = None
rows = []

for line in ann.read_text(errors="replace").splitlines():
    if not line or line.startswith("##"):
        continue
    if line.startswith("#query"):
        header = line.split("\t")
        continue
    if line.startswith("#"):
        continue
    rows.append(line.split("\t"))

if header is None:
    raise SystemExit(f"[ERROR] no #query header in {ann}")

header = ["GOs" if x == "GO" else x for x in header]

with out.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(header)
    for r in rows:
        w.writerow(r + [""] * (len(header) - len(r)))

print("[OK] wrote", out)
PY
}

valid_faa () {
    local f="$1"

    [ -s "$f" ] || {
        echo "[WARN] Skip empty FAA: $f"
        return 1
    }

    local first
    first="$(grep -m1 -v '^[[:space:]]*$' "$f" || true)"

    [[ "$first" == ">"* ]] || {
        echo "[WARN] Skip invalid FAA: $f"
        return 1
    }

    return 0
}

run_emapper () {
    local faa="$1"
    local csv="$2"
    local base
    local outdir

    base="$(basename "$faa" .faa)"
    outdir="02_work/01_emapper/$base"
    mkdir -p "$outdir"

    echo "[INFO] eggNOG: $faa"

    EGGNOG_DATA_DIR="$EGGNOG_DB" conda run -p "$EGGNOG_ENV" emapper.py \
      -i "$faa" \
      --itype proteins \
      -m diamond \
      --data_dir "$EGGNOG_DB" \
      --cpu "$THREADS" \
      --sensmode fast \
      --dmnd_iterate no \
      --block_size 0.5 \
      --index_chunks 4 \
      --evalue 0.001 \
      -o "$base" \
      --output_dir "$outdir" \
      --override

    convert_csv "$outdir/${base}.emapper.annotations" "$csv"
}

echo "[INFO] Running eggNOG-mapper..."

for faa in 00_input/01_eggNOG_functional_annotation/01_bacterial_protein_genes/*.faa; do
    [ -f "$faa" ] || continue
    valid_faa "$faa" || continue

    base="$(basename "$faa" .faa)"
    out1="00_input/01_eggNOG_functional_annotation/01_bacterial_protein_genes/${base}.csv"

    run_emapper "$faa" "$out1"
    test -s "$out1"

    if [[ "$base" == Re-annotated_* ]]; then
        out2="00_input/03_genbank/03_eggNOG/01_bacteria/${base}.csv"
        mkdir -p "$(dirname "$out2")"
        cp -f "$out1" "$out2"
        test -s "$out2"
    fi
done

for faa in 00_input/01_eggNOG_functional_annotation/02_plasmid_protein_genes/*.faa; do
    [ -f "$faa" ] || continue
    valid_faa "$faa" || continue

    base="$(basename "$faa" .faa)"
    out1="00_input/01_eggNOG_functional_annotation/02_plasmid_protein_genes/${base}.csv"

    run_emapper "$faa" "$out1"
    test -s "$out1"

    if [[ "$base" == PLRe-annotated_* ]]; then
        out2="00_input/03_genbank/03_eggNOG/02_plasmid/${base}.csv"
        mkdir -p "$(dirname "$out2")"
        cp -f "$out1" "$out2"
        test -s "$out2"
    fi
done

echo "[INFO] Running Snakemake master_all..."

conda run -p "$MODULE4_ENV" snakemake \
  --snakefile Snakefile.master \
  --cores "$THREADS" \
  --printshellcmds \
  --rerun-incomplete \
  master_all \
  2>&1 | tee 04_results/00_logs/master_snakemake.log

echo "[INFO] Running antiSMASH outside Snakemake..."

mkdir -p \
  00_input/05_BGCs \
  00_input/05_antismash \
  04_results/00_logs \
  04_results/01_final_tables \
  04_results/02_figures/03_BGCs

for gb in 00_input/04_hypothetical_genes/01_Re-annotated_bacterial_genbank/Re-annotated_*.gb; do
    [ -f "$gb" ] || continue

    id="$(basename "$gb" .gb | sed 's/Re-annotated_//')"
    outdir="00_input/05_antismash/reannotated_${id}"
    log="04_results/00_logs/antismash_reannotated_${id}.log"

    echo "[INFO] antiSMASH: $gb"

    rm -rf "$outdir"

    env -i \
      HOME="/home/ataa" \
      USER="ataa" \
      LOGNAME="ataa" \
      TMPDIR="/tmp" \
      PYTHONNOUSERSITE=1 \
      CONDA_PREFIX="$ANTISMASH_ENV" \
      PATH="$ANTISMASH_ENV/bin:/usr/bin:/bin" \
      "$ANTISMASH_ENV/bin/antismash" "$gb" \
        --taxon bacteria \
        --genefinding-tool prodigal \
        --cpus 2 \
        --output-dir "$outdir" \
        > "$log" 2>&1 || {
            echo "[FATAL] antiSMASH failed for $gb"
            cat "$log"
            exit 1
        }

    if [ -s "$outdir/regions.js" ]; then
        cp -f "$outdir/regions.js" "00_input/05_BGCs/${id}_reannotated_regions.js"
    else
        echo "[WARN] No regions.js for $id" | tee -a "$log"
        touch "00_input/05_BGCs/${id}_reannotated_regions.js"
    fi
done

python3 03_scripts/05_antiSMASH/antismash_simple_table.py \
  --in 00_input/05_BGCs \
  --out 04_results/01_final_tables/04_01_BGC_counts.tsv \
  --out_csv 04_results/01_final_tables/04_01_BGC_counts.csv

export BGC_TSV="04_results/01_final_tables/04_01_BGC_counts.tsv"
export OUTDIR="04_results/02_figures/03_BGCs"

conda run -p /home/ataa/miniconda3/envs/module4 \
Rscript 03_scripts/05_antiSMASH/antiSMASH_top_hits_per_species.R

echo "[INFO] Packaging outputs..."

tar -czf "$RUN_DIR/module4_results.tar.gz" 04_results
tar -czf "$RUN_DIR/module4_work.tar.gz" 02_work 00_input
tar -czf "$RUN_DIR/module4_logs.tar.gz" 04_results/00_logs
tar -czf "$RUN_DIR/module4_completed_input.tar.gz" 00_input

echo "[INFO] Final archives created:"
ls -lh "$RUN_DIR"/module4_*.tar.gz

echo "[INFO] DONE"
