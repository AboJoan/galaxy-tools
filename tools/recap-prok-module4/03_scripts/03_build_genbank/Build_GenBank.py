#!/usr/bin/env python3
import sys
import csv
import argparse
import re
from pathlib import Path
from collections import OrderedDict
from datetime import datetime

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


# -------------------------
# Utility
# -------------------------

def norm_empty(x) -> str:
    if x is None:
        return ""
    x = str(x).strip()
    if x.lower() in {"", "na", "none", "-", "nan"}:
        return ""
    return x

def strip_version(acc: str) -> str:
    return re.sub(r"\.\d+$", "", acc.strip())

def safe_first_token(x: str) -> str:
    return norm_empty(x).split()[0]

def parse_tsv(path: str) -> list[dict]:
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append({k: (v if v is not None else "") for k, v in row.items()})
    return rows

def find_metadata_row(rows: list[dict], replicon_accession: str) -> dict:
    """
    Accept:
      - real replicon accession (NZ_CP039126.1 or NZ_CP039126)
      - OR a file-key like Re-annotated_01 / Re-annotated_03 / PLRe-annotated_02 / PLRe-annotated_031
    Priority:
      1) match Replicon_accession / Replicon_accession_norm (with/without version)
      2) else parse index from key and match metadata by Genome_ID or Inf + replicon role/set
    """
    q = norm_empty(replicon_accession)
    q0 = strip_version(q)

    # ---- 1) direct accession match ----
    for r in rows:
        ra = strip_version(norm_empty(r.get("Replicon_accession")))
        rn = strip_version(norm_empty(r.get("Replicon_accession_norm")))
        if q0 and (q0 == ra or q0 == rn):
            return r

    # ---- 2) treat it as a key like Re-annotated_01 / Re-annotated_03 / PLRe-annotated_02 / PLRe-annotated_031 ----
    # Extract the first 2-digit index we see at the end: _01, _02, _03 ...
    m = re.search(r"_(\d{2})(?:\D.*)?$", q)
    if not m:
        raise ValueError(f"[ERROR] Replicon_accession not found in metadata: {replicon_accession} (normalized={q0})")

    idx = m.group(1)  # "01"
    is_plasmid_key = bool(re.match(r"^(PLRe-|Putative_PLRe-)", q))

    # Candidates matching by Genome_ID or Inf
    candidates = []
    for r in rows:
        genome_id = norm_empty(r.get("Genome_ID"))          # Re-annotated_01...02....03
        inf       = norm_empty(r.get("Inf"))                # our_genome_01
        rep_set   = norm_empty(r.get("Replicon_set"))       # bacterial_genome / plasmid ...
        rep_role  = norm_empty(r.get("Replicon_role"))      # chromosome_1 / plasmid ...
        # match index
        ok_idx = (
            genome_id.endswith(f"_{idx}") or
            inf.endswith(f"_{idx}")
        )
        if not ok_idx:
            continue

        if is_plasmid_key:
            # be permissive: anything that looks like plasmid
            if ("plasmid" in rep_set.lower()) or ("plasmid" in rep_role.lower()):
                candidates.append(r)
        else:
            # bacterial/chr
            if ("bacterial" in rep_set.lower()) or ("chromosome" in rep_role.lower()) or ("genome" in rep_set.lower()):
                candidates.append(r)

    # If strict filter found nothing, fall back to any idx match (still better than failing)
    if not candidates:
        for r in rows:
            genome_id = norm_empty(r.get("Genome_ID"))
            inf       = norm_empty(r.get("Inf"))
            if genome_id.endswith(f"{idx}") or inf.endswith(f"{idx}"):
                candidates.append(r)

    if not candidates:
        raise ValueError(
            f"[ERROR] Could not map key '{replicon_accession}' to a metadata row using Genome_ID/Inf index _{idx}."
        )

    # If multiple rows, choose the first (metadata should be consistent)
    return candidates[0]
    raise ValueError(f"[ERROR] Replicon_accession not found in metadata: {replicon_accession} (normalized={q0})")


# -------------------------
# FASTA / Prodigal parsing
# -------------------------

def read_fasta_as_ordered_records(path: str) -> OrderedDict:
    recs = OrderedDict()
    for rec in SeqIO.parse(path, "fasta"):
        recs[rec.id] = rec
    return recs

def parse_prodigal_faa(path: str):
    id2aa = {}
    gffid2protTok = {}
    for rec in SeqIO.parse(path, "fasta"):
        header_tok = rec.id
        desc = rec.description
        m = re.search(r'ID=([^; \t]+)', desc)
        if m:
            gff_id = m.group(1)
            gffid2protTok[gff_id] = header_tok
        id2aa[header_tok] = rec.seq
    return id2aa, gffid2protTok

def parse_gff_cds(path: str):
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype.lower() != "cds":
                continue

            a = {}
            for kv in attrs.split(";"):
                if not kv:
                    continue
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    a[k] = v
                else:
                    a[kv] = ""

            gid = a.get("ID") or a.get("locus_tag") or a.get("protein_id") or a.get("Name") or ""
            yield {
                "seqid": seqid,
                "start": int(start),
                "end":   int(end),
                "strand": strand,
                "phase": phase,
                "attrs": a,
                "gff_id": gid
            }


# -------------------------
# eggNOG parsing
# -------------------------

def load_eggnog_csv(path: str) -> dict:
    mapping = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        cols = [c.strip() for c in reader.fieldnames or []]
        qcol = "query" if "query" in cols else ("#query" if "#query" in cols else None)
        if qcol is None:
            raise ValueError("eggNOG CSV must have a 'query' or '#query' column.")

        for row in reader:
            key = safe_first_token(row.get(qcol, ""))
            if not key:
                continue
            mapping[key] = {
                "COG_category":  norm_empty(row.get("COG_category")),
                "Description":   norm_empty(row.get("Description")),
                "Preferred_name": norm_empty(row.get("Preferred_name")),
                "GOs":           norm_empty(row.get("GOs")) or norm_empty(row.get("GO_Term")),
                "EC":            norm_empty(row.get("EC")),
                "KEGG_ko":       norm_empty(row.get("KEGG_ko")) or norm_empty(row.get("KEGG_KO")),
                "PFAMs":         norm_empty(row.get("PFAMs")),
                "KEGG_Pathway":  norm_empty(row.get("KEGG_Pathway")),
                "KEGG_Module":   norm_empty(row.get("KEGG_Module")),
            }
    return mapping


# -------------------------
# PFAM DUF map (optional)
# -------------------------

UNKNOWN_WORDS_RE = re.compile(
    r"\b(unknown function|unknowen function|function unknown|uncharacteri[sz]ed|hypothetical)\b",
    flags=re.I
)

def load_pfam_duf_map(tsv_path: str | None) -> dict:
    """
    Flexible loader:
    returns map where keys can be:
      - PFxxxxx
      - DUF#### (if present)
    values: description string
    """
    if not tsv_path:
        return {}
    p = Path(tsv_path)
    if not p.exists() or p.stat().st_size == 0:
        return {}

    with open(p, newline="") as fh:
        # try tab-delimited with header
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            return {}

        # Find best columns
        cols = [c.strip() for c in reader.fieldnames]
        # candidate id columns
        id_cols = [c for c in cols if c.lower() in {"pfam", "pfam_id", "id", "name", "accession", "acc", "pfam_acc"}]
        desc_cols = [c for c in cols if c.lower() in {"desc", "description", "pfam_desc"}]

        mp = {}
        for row in reader:
            # pick id
            rid = ""
            for c in id_cols:
                rid = norm_empty(row.get(c))
                if rid:
                    break
            if not rid:
                continue

            # pick desc
            rdesc = ""
            for c in desc_cols:
                rdesc = norm_empty(row.get(c))
                if rdesc:
                    break

            rid = rid.strip()
            if rid:
                mp[rid] = rdesc

        return mp


# -------------------------
# Functional logic (YOUR FINAL RULES)
# -------------------------

BAD_PRODUCT_PATTERNS = [
    r"\bpsort\b",
    r"\blocation\b",
    r"\bsubcellular\b",
    r"\bcellular\s+location\b",
    r"\bextracellular\b",
    r"\bperiplasm(ic)?\b",
    r"\bmembrane\s+location\b",
    r"\bunknown\b",
    r"\bunknowen\b",
    r"\bfunction\s+unknown\b",
    r"\buncharacteri[sz]ed\b",
]
bad_product_re = re.compile("|".join(BAD_PRODUCT_PATTERNS), flags=re.I)

def split_tokens(x: str) -> list[str]:
    x = norm_empty(x)
    if not x:
        return []
    return [t for t in re.split(r"[,\s;|]+", x) if t.strip()]

def is_only_cog_s(note_txt: str) -> bool:
    """
    True if note contains ONLY COG=S and nothing else.
    Example note: "COG=S"
    or "COG=S; " (minor whitespace)
    """
    note_txt = (note_txt or "").strip()
    if not note_txt:
        return False
    # Remove spaces/newlines
    t = re.sub(r"\s+", "", note_txt)
    # allow "COG=S" alone or with trailing separators
    return t in {"COG=S", "COG=S;", ";COG=S", ";COG=S;", "COG=S;;"}

def build_note(ann: dict, strip_ko_prefix: bool = False) -> str | None:
    fields = []

    cog = norm_empty(ann.get("COG_category"))
    if cog:
        fields.append(f"COG={cog}")

    ko = norm_empty(ann.get("KEGG_ko"))
    if ko:
        if strip_ko_prefix:
            ko = ko.replace("ko:", "")
        fields.append(f"KO={ko}")

    go = norm_empty(ann.get("GOs"))
    if go:
        fields.append(f"GO={go}")

    ec = norm_empty(ann.get("EC"))
    if ec:
        fields.append(f"EC={ec}")

    pf = norm_empty(ann.get("PFAMs"))
    if pf:
        fields.append(f"PFAMs={pf}")

    pw = norm_empty(ann.get("KEGG_Pathway"))
    if pw:
        fields.append(f"KEGG_Pathway={pw}")

    km = norm_empty(ann.get("KEGG_Module"))
    if km:
        fields.append(f"KEGG_Module={km}")

    return "; ".join(fields) if fields else None

def note_has_any_evidence_except_cog_s(note_txt: str) -> bool:
    """
    Evidence defined as presence of any of:
      KO=, GO=, EC=, PFAMs=, KEGG_Pathway=, KEGG_Module=, or COG=<not S>
    """
    if not note_txt:
        return False
    if "KO=" in note_txt or "GO=" in note_txt or "EC=" in note_txt or "PFAMs=" in note_txt or "KEGG_Pathway=" in note_txt or "KEGG_Module=" in note_txt:
        return True
    # COG evidence: any letter except S
    return bool(re.search(r"(?:^|[;\s])COG=([A-Z])(?:$|[;\s])", note_txt, flags=re.I)) and not bool(re.search(r"(?:^|[;\s])COG=S(?:$|[;\s])", note_txt, flags=re.I))

def extract_cog_from_note(note_txt: str) -> str:
    m = re.search(r"(?:^|[;\s])COG=([A-Za-z])(?:$|[;\s])", note_txt or "")
    return (m.group(1).upper() if m else "")

def duf_only_and_unknown(pfam_field: str, pfam_map: dict) -> bool:
    """
    If PFAMs contains ONLY DUF tokens (or PFAMs that map to DUF),
    then consult mapping description:
      - if description has unknown/uncharacterized/hypothetical => treat as hypothetical
      - else treat as functional
    If mapping missing => conservative: unknown => hypothetical
    """
    toks = split_tokens(pfam_field)
    if not toks:
        return False

    # consider DUF-only if every token starts with DUF OR is PFxxxxx that maps to a DUF-like name in map
    def is_duf_token(t: str) -> bool:
        tu = t.upper()
        if tu.startswith("DUF"):
            return True
        # PFxxxxx case: if in map and its ID/name indicates DUF
        if tu.startswith("PF") and tu in pfam_map:
            desc = pfam_map.get(tu, "") or ""
            # if id in map is DUF? we only have desc; accept if desc mentions DUF
            return "DUF" in desc.upper()
        return False

    if not all(is_duf_token(t) for t in toks):
        return False  # not DUF-only

    # DUF-only: decide unknown vs functional using mapping
    # If ANY DUF desc indicates unknown/uncharacterized/hypothetical -> hypothetical
    # If all mapped and none unknown -> functional
    # If mapping missing -> hypothetical (unknown)
    any_missing = False
    for t in toks:
        key = t.upper()
        desc = ""
        if key in pfam_map:
            desc = pfam_map.get(key, "") or ""
        else:
            # sometimes map uses PFxxxxx while PFAMs contains DUF#### or vice versa
            # try to locate by scanning keys containing the token
            # (cheap fallback)
            hits = [k for k in pfam_map.keys() if k.upper() == key]
            if hits:
                desc = pfam_map.get(hits[0], "") or ""
            else:
                any_missing = True

        if desc and UNKNOWN_WORDS_RE.search(desc):
            return True

    if any_missing:
        return True

    return False

def choose_functional_product_name(ann: dict) -> str:
    """
    Must not output unknown/psort/etc.
    Prefer Preferred_name then Description; else 'putative protein'.
    """
    pref = norm_empty(ann.get("Preferred_name"))
    desc = norm_empty(ann.get("Description"))

    def ok_name(x: str) -> bool:
        if not x or x == "-":
            return False
        if bad_product_re.search(x):
            return False
        return True

    if ok_name(pref):
        return pref
    if ok_name(desc):
        return desc if len(desc) <= 200 else (desc[:197] + "...")
    return "putative protein"

def decide_product_and_note(ann: dict, pfam_map: dict, strip_ko_prefix: bool) -> tuple[str, str | None]:
    """
    FINAL DECISION:
      hypothetical ONLY if:
        - note is empty/None   OR
        - note only contains COG=S  OR
        - ONLY evidence is PFAMs and it's DUF-only AND mapping says unknown/uncharacterized/hypothetical
      otherwise functional (never hypothetical).
    """
    note = build_note(ann, strip_ko_prefix=strip_ko_prefix)
    note_txt = (note or "").replace("\n", " ").strip()

    # Case 1: note empty
    if not note_txt:
        return "hypothetical protein", None  # no note

    # Case 2: note only COG=S
    if is_only_cog_s(note_txt):
        return "hypothetical protein", note

    # Identify evidence types present
    cog = extract_cog_from_note(note_txt)
    has_ko = "KO=" in note_txt
    has_go = "GO=" in note_txt
    has_ec = "EC=" in note_txt
    has_pw = "KEGG_Pathway=" in note_txt
    has_km = "KEGG_Module=" in note_txt
    has_pfam = "PFAMs=" in note_txt

    # If any real evidence besides PFAM exists => functional
    if (cog and cog != "S") or has_ko or has_go or has_ec or has_pw or has_km:
        return choose_functional_product_name(ann), note

    # Now we are in a narrow case:
    # Note exists, but NO KO/GO/EC/Pathway/Module, and COG is either empty or S (but not only COG=S), leaving PFAMs as possible evidence.
    # If PFAM exists and is DUF-only -> consult map
    if has_pfam:
        pfam_field = norm_empty(ann.get("PFAMs"))
        if duf_only_and_unknown(pfam_field, pfam_map):
            return "hypothetical protein", note
        else:
            return choose_functional_product_name(ann), note

    # If we reach here: note exists but doesn't contain KO/GO/EC/PFAM/pathway/module and COG isn't non-S (already handled),
    # so treat as hypothetical (very rare).
    return "hypothetical protein", note


# -------------------------
# Metadata -> GenBank header
# -------------------------

def build_taxonomy_list(meta: dict) -> list[str]:
    ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    out = []
    for r in ranks:
        v = norm_empty(meta.get(r))
        if v:
            out.append(v)
    return out

def format_locus_date_today() -> str:
    return datetime.now().strftime("%d-%b-%Y").upper()


# -------------------------
# Main
# -------------------------

def main():
    print(f"[RUN] Using script: {Path(__file__).resolve()}", file=sys.stderr)

    ap = argparse.ArgumentParser(
        description="Build GenBank from assembly FASTA, Prodigal GFF/FAA, eggNOG CSV, and metadata TSV."
    )
    ap.add_argument("--metadata", required=True, help="Replicon metadata TSV")
    ap.add_argument("--replicon_accession", required=True, help="e.g. NZ_CP039126 or NZ_CP039126.1")

    ap.add_argument("--fna", required=True, help="ASSEMBLY FASTA (genome/contigs), NOT genes_cds.fna")
    ap.add_argument("--faa", required=True, help="Prodigal proteins .faa")
    ap.add_argument("--gff", required=True, help="Prodigal .gff")
    ap.add_argument("--eggnog_csv", required=True, help="eggNOG-mapper CSV")
    ap.add_argument("--out_gbk", required=True, help="Output GenBank file")

    ap.add_argument("--pfam_duf_map", default=None, help="Optional TSV mapping for DUF/PFAM descriptions")
    ap.add_argument("--topology", choices=["circular", "linear"], default="circular")
    ap.add_argument("--locus_prefix", default=None, help="Prefix for locus_tag, e.g. BAS01")
    ap.add_argument("--taxon", type=int, default=None, help="NCBI taxid (optional)")
    ap.add_argument("--host", default=None, help="Override host (optional)")
    ap.add_argument("--strip-ko-prefix", action="store_true", help="KO=Kxxxxx instead of KO=ko:Kxxxxx")
    ap.add_argument("--transl-table", type=int, default=11)

    args = ap.parse_args()

    # Load metadata row
    meta_rows = parse_tsv(args.metadata)
    meta = find_metadata_row(meta_rows, args.replicon_accession)

    organism = norm_empty(meta.get("Host")) or norm_empty(meta.get("Genus")) or "Unknown organism"
    strain = norm_empty(meta.get("Strain"))
    definition = strain if strain else f"{organism} genome"

    acc_norm = strip_version(
        norm_empty(meta.get("Replicon_accession_norm")) or
        norm_empty(meta.get("Replicon_accession")) or
        args.replicon_accession
    )
    record_id = acc_norm

    genome = read_fasta_as_ordered_records(args.fna)
    prot_seqs, gffid2protTok = parse_prodigal_faa(args.faa)
    egg = load_eggnog_csv(args.eggnog_csv)

    # PFAM DUF map
    pfam_map = load_pfam_duf_map(args.pfam_duf_map)

    # Build SeqRecords per contig
    contig_records = OrderedDict()
    locus_date = format_locus_date_today()
    taxonomy_list = build_taxonomy_list(meta)

    for contig_id, contig_rec in genome.items():
        if len(genome) == 1:
            rec_id = record_id
            rec_name = record_id
        else:
            rec_id = contig_id
            rec_name = contig_id

        r = SeqRecord(
            contig_rec.seq,
            id=rec_id,
            name=rec_name,
            description=definition
        )

        r.annotations["molecule_type"] = "DNA"
        r.annotations["topology"] = args.topology
        r.annotations["data_file_division"] = "BCT"
        r.annotations["date"] = locus_date
        if taxonomy_list:
            r.annotations["taxonomy"] = taxonomy_list

        r.annotations["comment"] = (
            "Re-annotated genome using Prodigal (CDS prediction) and eggNOG-mapper "
            "(functional annotation). GenBank file constructed manually for comparative annotation."
        )

        contig_records[contig_id] = r

    # SOURCE + rep_origin
    for contig_id, rec in contig_records.items():
        src_qual = {
            "organism": [organism],
            "mol_type": ["genomic DNA"],
            "topology": [args.topology],
        }
        if strain:
            src_qual["strain"] = [strain]
        if args.taxon:
            src_qual["db_xref"] = [f"taxon:{args.taxon}"]

        host_val = norm_empty(args.host) or norm_empty(meta.get("Source"))
        if host_val:
            src_qual["host"] = [host_val]

        rec.features.append(
            SeqFeature(
                FeatureLocation(0, len(rec.seq), strand=1),
                type="source",
                qualifiers=src_qual
            )
        )

        rep_txt = norm_empty(meta.get("Replication origins"))
        if rep_txt:
            nums = [int(re.sub(r"[^\d]", "", x)) for x in re.findall(r"[\d,]+", rep_txt)]
            if len(nums) >= 2:
                s, e = nums[0], nums[1]
                start0 = max(0, s - 1)
                end0 = min(len(rec.seq), e)
            else:
                start0, end0 = 0, 1

            rec.features.append(
                SeqFeature(
                    FeatureLocation(start0, end0, strand=1),
                    type="rep_origin",
                    qualifiers={
                        "note": [f"Predicted oriC; {rep_txt}"],
                        "inference": ["prediction:Ori-Finder 2022"]
                    }
                )
            )

    # Add CDS features
    n_cds = 0
    for feat in parse_gff_cds(args.gff):
        seqid = feat["seqid"]
        if seqid not in contig_records:
            sid = seqid.split()[0]
            if sid in contig_records:
                seqid = sid
            else:
                continue

        gff_id = feat["gff_id"]
        protTok = gffid2protTok.get(gff_id, None)
        aa_seq = prot_seqs.get(protTok, None)
        ann = egg.get(protTok, {}) if protTok else {}

        # ---- FINAL forced rules (your rule-set) ----
        product, note = decide_product_and_note(
            ann=ann,
            pfam_map=pfam_map,
            strip_ko_prefix=args.strip_ko_prefix
        )

        # location
        start0 = feat["start"] - 1
        end0 = feat["end"]
        strand = -1 if feat["strand"] == "-" else 1
        loc = FeatureLocation(start0, end0, strand=strand)

        # locus_tag
        if args.locus_prefix:
            locus_tag = f"{args.locus_prefix}_{n_cds + 1:05d}"
        else:
            locus_tag = gff_id if gff_id else f"cds_{n_cds + 1}"

        qualifiers = {
            "locus_tag": [locus_tag],
            "product": [product],
            "transl_table": [str(args.transl_table)],
            "protein_id": [locus_tag],
        }

        if note:
            qualifiers["note"] = [note]

        if aa_seq is not None and len(aa_seq) > 0:
            qualifiers["translation"] = [str(aa_seq)]

        contig_records[seqid].features.append(SeqFeature(loc, type="CDS", qualifiers=qualifiers))
        n_cds += 1

    if n_cds == 0:
        print(
            "WARNING: No CDS features were written. "
            "Check matching between GFF IDs and FAA headers (ID=...), and eggNOG query IDs.",
            file=sys.stderr
        )

    out_path = Path(args.out_gbk)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as out:
        SeqIO.write(list(contig_records.values()), out, "genbank")

    print(f"✔ Wrote: {out_path}")
    print(f"Contigs: {len(contig_records)}   CDS features: {n_cds}")


if __name__ == "__main__":
    main()
