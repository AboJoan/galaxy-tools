import matplotlib
matplotlib.use("Agg")
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# =========================================================
# Genome metadata
# =========================================================

GENOME_NAMES = {
    "01": {"full_name": "Blautia producta", "short_label": "B. producta"},
    "02": {"full_name": "Bacteroides thetaiotaomicron", "short_label": "B. thetaiotaomicron"},
    "03": {"full_name": "Clostridium butyricum", "short_label": "C. butyricum"},
    "04": {"full_name": "Escherichia coli", "short_label": "E. coli"},
    "05": {"full_name": "Lactiplantibacillus plantarum", "short_label": "L. plantarum"},
    "06": {"full_name": "Thomasclavelia ramosa", "short_label": "T. ramosa"},
    "07": {"full_name": "Anaerostipes caccae", "short_label": "A. caccae"},
    "08": {"full_name": "Bifidobacterium longum", "short_label": "B. longum"},
    "09": {"full_name": "Bacteroides thetaiotaomicron plasmid", "short_label": "B. thetaiotaomicron plasmid"},
    "10": {"full_name": "Clostridium butyricum plasmid 01", "short_label": "C. butyricum plasmid 01"},
    "11": {"full_name": "Clostridium butyricum plasmids 02 and 03", "short_label": "C. butyricum plasmids 02/03"},
    "12": {"full_name": "Lactiplantibacillus plantarum plasmid", "short_label": "L. plantarum plasmid"},
    "13": {"full_name": "Bifidobacterium longum plasmid", "short_label": "B. longum plasmid"},
}

# =========================================================
# Input / output
# =========================================================

folder = Path("02_work/03_pairs_classified")
files = sorted(folder.glob("*_pairs_classified.tsv"))

outdir = Path("04_results/02_figures/03_cds_mapping_interpretation")
outtable = Path("04_results/01_final_tables")

outdir.mkdir(parents=True, exist_ok=True)
outtable.mkdir(parents=True, exist_ok=True)

print(f"Found {len(files)} input files")

# =========================================================
# Helper functions
# =========================================================

def cds_length(start, end):
    return abs(int(end) - int(start)) + 1


def classify_direction(row):
    """
    Strand-aware interpretation of CDS boundary changes.
    Reference = NCBI CDS mapped/projected coordinates.
    Re-annotation = Prodigal-based CDS coordinates.
    """

    strand = str(row["ref_strand"]).strip()

    ref_start = int(row["ref_start"])
    ref_end = int(row["ref_end"])
    re_start = int(row["reann_start"])
    re_end = int(row["reann_end"])

    ref_len = cds_length(ref_start, ref_end)
    re_len = cds_length(re_start, re_end)

    if ref_len == re_len:
        return "same"

    if strand == "+":
        if re_start < ref_start:
            return "extension_after_reannotation"
        if re_start > ref_start:
            return "truncation_after_reannotation"

    if strand == "-":
        if re_end > ref_end:
            return "extension_after_reannotation"
        if re_end < ref_end:
            return "truncation_after_reannotation"

    if re_len > ref_len:
        return "extension_after_reannotation"

    if re_len < ref_len:
        return "truncation_after_reannotation"

    return "other"


def safe_percent(count, total):
    return round((count / total) * 100, 2) if total else 0


# =========================================================
# Main analysis
# =========================================================

summaries = []

for file in files:
    print(f"Processing: {file.name}")

    df = pd.read_csv(file, sep="\t")
    genome_id = file.stem.split("_")[0]

    info = GENOME_NAMES.get(
        genome_id,
        {
            "full_name": f"Genome {genome_id}",
            "short_label": f"Genome {genome_id}",
        }
    )

    sequence_type = "bacterial chromosome" if int(genome_id) <= 8 else "plasmid"
    total_pairs = len(df)

    exact_match = df["class"].astype(str).str.contains(
        "MATCH_ALL", case=False, na=False
    ).sum()

    locus_shift = df["class"].astype(str).str.contains(
        "LOC", case=False, na=False
    ).sum()

    mapped = df.dropna(
        subset=["ref_start", "ref_end", "reann_start", "reann_end"]
    ).copy()

    if not mapped.empty:
        mapped["ref_len"] = mapped.apply(
            lambda r: cds_length(r["ref_start"], r["ref_end"]),
            axis=1
        )

        mapped["reann_len"] = mapped.apply(
            lambda r: cds_length(r["reann_start"], r["reann_end"]),
            axis=1
        )

        mapped["length_diff"] = mapped["reann_len"] - mapped["ref_len"]
        mapped["direction_class"] = mapped.apply(classify_direction, axis=1)

        extensions = (
            mapped["direction_class"] == "extension_after_reannotation"
        ).sum()

        truncations = (
            mapped["direction_class"] == "truncation_after_reannotation"
        ).sum()

        mean_length_change = round(mapped["length_diff"].mean(), 2)
    else:
        extensions = 0
        truncations = 0
        mean_length_change = 0

    summaries.append({
        "sort_id": int(genome_id),
        "Genome": info["short_label"],
        "ST_{a}": sequence_type,

        "EM_{b}": exact_match,
        "EM_(%)_{c}": safe_percent(exact_match, total_pairs),

        "EXT_{d}": extensions,
        "EXT_(%)_{e}": safe_percent(extensions, total_pairs),

        "TRUNC_{f}": truncations,
        "TRUNC_(%)_{g}": safe_percent(truncations, total_pairs),

        "LS_{h}": locus_shift,
        "Mean_∆CDS_(bp)_{i}": mean_length_change,
    })


summary_df = pd.DataFrame(summaries)
summary_df = summary_df.sort_values("sort_id")

final_df = summary_df[
    [
        "Genome",
        "ST_{a}",
        "EM_{b}",
        "EM_(%)_{c}",
        "EXT_{d}",
        "EXT_(%)_{e}",
        "TRUNC_{f}",
        "TRUNC_(%)_{g}",
        "LS_{h}",
        "Mean_∆CDS_(bp)_{i}",
    ]
].copy()

final_df.to_csv(
    outtable / "CDS_mapping_summary_NCBI_reference_vs_Prodigal_reannotation.tsv",
    sep="\t",
    index=False
)

bacteria_df = summary_df[
    summary_df["ST_{a}"] == "bacterial chromosome"
].copy()

plasmid_df = summary_df[
    summary_df["ST_{a}"] == "plasmid"
].copy()

# =========================================================
# Plot functions
# =========================================================

def add_segment_label(ax, x, y_center, percent, count, total, color="black"):
    ax.text(
        x,
        y_center,
        f"{percent:.1f}%\n(n={count})",
        ha="center",
        va="center",
        fontsize=8,
        color=color,
        fontweight="bold"
    )


def plot_stacked_percent(df, title, filename, xlabel):
    fig, ax = plt.subplots(figsize=(13, 6))

    x = list(range(len(df)))

    exact = df["EM_(%)_{c}"]
    ext = df["EXT_(%)_{e}"]
    trunc = df["TRUNC_(%)_{g}"]

    ax.bar(x, exact, label="Exact matches")
    ax.bar(x, ext, bottom=exact, label="CDS extensions after re-annotation")
    ax.bar(x, trunc, bottom=exact + ext, label="CDS truncations after re-annotation")

    for i, row in df.reset_index(drop=True).iterrows():
        total = row["EM_{b}"] + row["EXT_{d}"] + row["TRUNC_{f}"]

        if row["EM_(%)_{c}"] >= 5:
            add_segment_label(
                ax,
                i,
                row["EM_(%)_{c}"] / 2,
                row["EM_(%)_{c}"],
                row["EM_{b}"],
                total,
                color="white"
            )

        if row["EXT_(%)_{e}"] >= 5:
            add_segment_label(
                ax,
                i,
                row["EM_(%)_{c}"] + row["EXT_(%)_{e}"] / 2,
                row["EXT_(%)_{e}"],
                row["EXT_{d}"],
                total,
                color="black"
            )

        if row["TRUNC_(%)_{g}"] >= 5:
            add_segment_label(
                ax,
                i,
                row["EM_(%)_{c}"] + row["EXT_(%)_{e}"] + row["TRUNC_(%)_{g}"] / 2,
                row["TRUNC_(%)_{g}"],
                row["TRUNC_{f}"],
                total,
                color="black"
            )

    ax.set_ylim(0, 100)
    ax.set_ylabel("CDS proportion of mapped CDS pairs (%)")
    ax.set_xlabel(xlabel)
    ax.set_title(title)

    ax.set_xticks(x)
    ax.set_xticklabels(df["Genome"], rotation=45, ha="right")

    ax.legend(
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )

    plt.tight_layout()

    plt.savefig(outdir / f"{filename}.png", dpi=300, bbox_inches="tight")
    plt.savefig(outdir / f"{filename}.pdf", bbox_inches="tight")
    plt.close()


def plot_mean_length(df, title, filename, xlabel):
    fig, ax = plt.subplots(figsize=(12, 5))

    x = list(range(len(df)))

    ax.bar(x, df["Mean_∆CDS_(bp)_{i}"])
    ax.axhline(0, linewidth=1)

    for i, row in df.reset_index(drop=True).iterrows():
        value = row["Mean_∆CDS_(bp)_{i}"]
        ax.text(
            i,
            value,
            f"{value:.2f} bp",
            ha="center",
            va="bottom" if value >= 0 else "top",
            fontsize=8
        )

    ax.set_ylabel("Mean CDS length change (bp)\nProdigal re-annotation minus NCBI reference")
    ax.set_xlabel(xlabel)
    ax.set_title(title)

    ax.set_xticks(x)
    ax.set_xticklabels(df["Genome"], rotation=45, ha="right")

    plt.tight_layout()

    plt.savefig(outdir / f"{filename}.png", dpi=300, bbox_inches="tight")
    plt.savefig(outdir / f"{filename}.pdf", bbox_inches="tight")
    plt.close()

# =========================================================
# Final figures
# =========================================================

plot_stacked_percent(
    bacteria_df,
    "CDS loci comparison: NCBI reference annotations vs Prodigal-based re-annotations in bacterial chromosomes",
    "Figure_1A_bacterial_chromosomes_NCBI_vs_Prodigal_CDS_mapping",
    "Bacterial chromosome"
)

plot_stacked_percent(
    plasmid_df,
    "CDS loci comparison: NCBI reference annotations vs Prodigal-based re-annotations in plasmids",
    "Figure_1B_plasmids_NCBI_vs_Prodigal_CDS_mapping",
    "Plasmid sequence"
)

plot_mean_length(
    bacteria_df,
    "Mean CDS length change after Prodigal re-annotation in bacterial chromosomes",
    "Figure_2A_bacterial_chromosomes_mean_CDS_length_change",
    "Bacterial chromosome"
)

plot_mean_length(
    plasmid_df,
    "Mean CDS length change after Prodigal re-annotation in plasmids",
    "Figure_2B_plasmids_mean_CDS_length_change",
    "Plasmid sequence"
)

print("Done.")
print(f"Results saved in: {outdir}")
