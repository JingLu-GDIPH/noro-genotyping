#!/usr/bin/env python3
"""Download GenBank records for genome references and derive gene coordinates.

For each accession in ``gi_with_genotype.fasta`` / ``gii_with_genotype.fasta`` the
script fetches the corresponding GenBank record, parses its CDS features and
records the coordinates of every annotated gene (ORF1, VP1, VP2 ...).

Norovirus GenBank records almost never annotate RdRp as a separate CDS (it sits
inside the ORF1 polyprotein), so RdRp coordinates are derived by BLASTing the
pre-trimmed region reference library (``gi_rdrp`` / ``gii_rdrp``) against each
genome and taking the best subject span. VP1 coordinates are taken from GenBank
when present and otherwise recovered the same way from the VP1 region library.

The result is written to ``ref_seq/genome_gene_coordinates.tsv`` which the main
genotyping pipeline reads at run time to slice query sequences.
"""

import argparse
import csv
import io
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


REFERENCE_FILES = {
    ("GI", "genome"): ("gi_with_genotype.fasta", "gi_genome.fasta"),
    ("GI", "rdrp"): "gi_rdrp.fasta",
    ("GI", "vp1"): ("gi_VP1.fa", "gi_vp1.fasta"),
    ("GII", "genome"): ("gii_with_genotype.fasta", "gii_genome.fasta"),
    ("GII", "rdrp"): "gii_rdrp.fasta",
    ("GII", "vp1"): "gii_vp1.fasta",
}

# Canonical protein names used in the coordinate table. GenBank ``product`` /
# ``gene`` qualifiers are mapped onto these via ``match_capsid_product`` /
# ``match_polyprotein_product``.
CANONICAL_PROTEINS = ["ORF1", "RdRp", "VP1", "VP2"]

# VP1 / capsid product strings seen in norovirus GenBank records.
CAPSID_KEYWORDS = ("vp1", "major capsid", "capsid protein vp1", "viral capsid",
                   "58 kd capsid", "coat protein", "major structural")

# VP2 product strings.
MINOR_KEYWORDS = ("vp2", "minor structural", "minor capsid", "orf3", "minor protein")

# Norovirus ORF1 is translated as a single polyprotein that is proteolytically
# cleaved into six non-structural proteins (NS1/2 ... NS7). GenBank almost never
# annotates the individual cleavage products, so we carve them out of the ORF1
# span using the conserved relative cleavage-site positions of the Norwalk virus
# polyprotein (well established in the literature). Fractions are expressed as a
# proportion of the ORF1 nucleotide length, (start, end] half-open per protein.
# This lets off-region queries be annotated with their actual protein (e.g.
# NTPase, VPg, 3CLpro) rather than the generic "ORF1".
ORF1_PROTEIN_BOUNDARIES = [
    ("p48",      0.000, 0.260),   # NS1/2 (N-terminal)
    ("NTPase",   0.260, 0.410),   # NS3 (p41, NTPase/helicase-like)
    ("p22",      0.410, 0.490),   # NS4
    ("VPg",      0.490, 0.570),   # NS5 (genome-linked viral protein)
    ("3CLpro",   0.570, 0.680),   # NS6 (3C-like protease)
    ("RdRp",     0.680, 1.000),   # NS7 (overlaps the region-library RdRp)
]

# Genome regions outside ORF1/2/3 that GenBank sometimes annotates and that are
# useful labels for off-region queries.
UTR_FEATURE_TYPES = {"5'UTR", "3'UTR", "UTR"}

COORDINATE_COLUMNS = [
    "group",
    "accession",
    "genotype",
    "sequence_length",
    "protein",
    "start",
    "end",
    "strand",
    "source",
]

ACCESSION_RE = re.compile(r"(NC_\d+|[A-Z]{1,2}\d{5,8})", re.IGNORECASE)


def accession_from_header(header):
    """Return the NCBI accession embedded in a genome FASTA header.

    Genome headers look like ``GI.1_MH638228`` (genotype_accession); we fall
    back to the first whitespace-delimited token if no accession pattern matches.
    """
    match = ACCESSION_RE.search(header)
    return match.group(1).upper() if match else header.split()[0]


def genotype_from_header(header, group):
    """Parse the genotype (``GI.1`` / ``GII.NA2`` / ...) from a genome header."""
    token = header.split()[0]
    match = re.match(rf"({group})\.(NA\d+|\d+)_", token, re.IGNORECASE)
    if match:
        return f"{match.group(1).upper()}.{match.group(2).upper()}"
    match = re.search(rf"({group})\.(NA\d+|\d+)", token, re.IGNORECASE)
    if match:
        return f"{match.group(1).upper()}.{match.group(2).upper()}"
    return ""


def resolve_path(ref_dir, names):
    if isinstance(names, str):
        names = (names,)
    for name in names:
        path = ref_dir / name
        if path.is_file():
            return path
    expected = ", ".join(str(ref_dir / n) for n in names)
    raise FileNotFoundError(f"Reference file not found; expected one of: {expected}")


def match_capsid_product(text):
    """Return canonical name if ``text`` looks like a VP1/VP2 product string."""
    low = text.lower()
    if any(k in low for k in MINOR_KEYWORDS):
        return "VP2"
    if any(k in low for k in CAPSID_KEYWORDS):
        return "VP1"
    return None


def match_polyprotein_product(text):
    """Return ``ORF1`` for the non-structural polyprotein CDS."""
    low = text.lower()
    if "nonstructural polyprotein" in low or "non-structural polyprotein" in low \
            or "polyprotein" in low or "orf1" in low:
        return "ORF1"
    return None


def parse_genbank_cds(gb_text):
    """Parse a GenBank record and return ``{canonical_protein: (start, end, strand)}``.

    Coordinates are 1-based inclusive on the forward (genome) strand, matching
    BLAST subject coordinates. ORF1/VP1/VP2 are recovered from GenBank CDS
    annotations, the 5'/3' UTRs from their feature types; RdRp is derived later
    by region-library BLAST.
    """
    record = next(SeqIO.parse(io.StringIO(gb_text), "genbank"))
    coords = {}
    for feature in record.features:
        if feature.type == "CDS":
            product = feature.qualifiers.get("product", [""])[0]
            gene = feature.qualifiers.get("gene", [""])[0]
            haystack = f"{gene} {product}"
            canonical = match_capsid_product(product) or match_capsid_product(gene)
            if canonical is None:
                canonical = match_polyprotein_product(haystack)
            if canonical is None:
                continue
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = "+" if feature.location.strand in (1, None) else "-"
            # Keep the first occurrence; norovirus genomes have single copies.
            coords.setdefault(canonical, (start, end, strand))
        elif feature.type in UTR_FEATURE_TYPES:
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = "+" if feature.location.strand in (1, None) else "-"
            # Label 5'/3' ends by position relative to the genome.
            mid = (start + end) / 2.0
            label = "5'UTR" if mid <= len(record.seq) / 2 else "3'UTR"
            coords.setdefault(label, (start, end, strand))
    coords["_length"] = len(record.seq)
    coords["_accession"] = record.id.split(".")[0].upper()
    return coords


def fetch_genbank(accessions, cache_dir, email, api_key, retries, batch_size):
    """Download GenBank records for ``accessions`` with caching and retries.

    Records are written to ``cache_dir/<accession>.gb``. A dict
    ``{accession: genbank_text}`` is returned for the records that could be
    retrieved; missing accessions are silently omitted (the caller handles them
    via the region-library fallback).
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    cache_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    pending = []
    for accession in accessions:
        cache_path = cache_dir / f"{accession}.gb"
        if cache_path.is_file() and cache_path.stat().st_size > 0:
            results[accession] = cache_path.read_text()
        else:
            pending.append(accession)

    sleep_seconds = 0.11 if api_key else 0.34
    for batch_start in range(0, len(pending), batch_size):
        batch = pending[batch_start:batch_start + batch_size]
        text = None
        for attempt in range(1, retries + 1):
            try:
                handle = Entrez.efetch(
                    db="nucleotide",
                    id=",".join(batch),
                    rettype="gbwithparts",
                    retmode="text",
                )
                text = handle.read()
                break
            except Exception as exc:  # noqa: BLE001 - network errors are varied
                if attempt == retries:
                    sys.stderr.write(
                        f"  WARN: failed to fetch batch starting {batch[0]} "
                        f"after {retries} attempts: {exc}\n"
                    )
                    time.sleep(sleep_seconds)
                else:
                    time.sleep(2.0 * attempt)
        if not text:
            continue
        # NCBI returns a concatenation of LOCUS-separated records.
        for record_text in split_genbank_records(text):
            accession = parse_locus_accession(record_text)
            if accession is None:
                continue
            cache_path = cache_dir / f"{accession}.gb"
            cache_path.write_text(record_text)
            results[accession] = record_text
        time.sleep(sleep_seconds)
    return results


def split_genbank_records(text):
    """Yield individual LOCUS..// GenBank records from a concatenated string."""
    records = []
    current = []
    for line in text.splitlines():
        current.append(line)
        if line.strip() == "//":
            records.append("\n".join(current) + "\n")
            current = []
    if current and any(line.startswith("LOCUS") for line in current):
        records.append("\n".join(current) + "\n")
    return records


def parse_locus_accession(record_text):
    """Pull the accession out of a single GenBank record's LOCUS/VERSION line."""
    for line in record_text.splitlines():
        if line.startswith("ACCESSION"):
            parts = line.split()
            if len(parts) >= 2:
                return parts[1].upper()
    for line in record_text.splitlines():
        if line.startswith("VERSION"):
            parts = line.split()
            if len(parts) >= 2:
                return parts[1].split(".")[0].upper()
    return None


def run(command):
    subprocess.run(
        [str(item) for item in command],
        check=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        text=True,
    )


def region_spans_for_genomes(genome_records, region_records, work_dir, prefix):
    """BLAST every genome against one region library and pick the best hit each.

    Building the database once for the region library and querying all genomes in
    a single ``blastn`` invocation is dramatically faster than one makeblastdb +
    blastn pair per genome (the previous per-genome approach rebuilt the database
    ~3000 times). Genomes are the database (subject); region refs are the query,
    so the returned ``sstart``/``send`` fall on each genome.

    Returns ``{genome_accession: (start, end, identity_pct)}`` (1-based, inclusive,
    forward strand). Genomes with no hit above the e-value threshold are omitted.
    """
    genome_fasta = work_dir / f"{prefix}_genomes.fasta"
    region_fasta = work_dir / f"{prefix}_region.fasta"
    database = work_dir / f"{prefix}_db"
    output = work_dir / f"{prefix}_hits.tsv"

    SeqIO.write([rec for _, _, _, rec in genome_records], genome_fasta, "fasta")
    SeqIO.write(region_records, region_fasta, "fasta")
    run(["makeblastdb", "-in", genome_fasta, "-dbtype", "nucl", "-out", database])
    run([
        "blastn", "-task", "blastn",
        "-query", region_fasta, "-db", database,
        "-max_target_seqs", str(max(1, len(genome_records))),
        "-max_hsps", "1", "-evalue", "1e-20",
        "-outfmt", "6 qseqid sseqid pident length sstart send bitscore evalue",
        "-out", output,
    ])
    spans = {}
    with output.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 8:
                continue
            genome_header = row[1]
            bitscore = float(row[6])
            sstart, send = int(row[4]), int(row[5])
            span = (min(sstart, send), max(sstart, send), float(row[2]), bitscore)
            current = spans.get(genome_header)
            if current is None or bitscore > current[3]:
                spans[genome_header] = span
    # Normalise keys from the genome FASTA header back to accession.
    by_accession = {}
    for group, accession, _genotype, genome_record in genome_records:
        span = spans.get(genome_record.id)
        if span is not None:
            by_accession[accession] = span[:3]
    return by_accession


def load_genome_records(ref_dir):
    """Return ``[(group, accession, genotype, SeqRecord), ...]`` for the genomes."""
    records = []
    for group in ("GI", "GII"):
        path = resolve_path(ref_dir, REFERENCE_FILES[(group, "genome")])
        for record in SeqIO.parse(path, "fasta"):
            accession = accession_from_header(record.id)
            genotype = genotype_from_header(record.id, group)
            clean = Seq(str(record.seq).replace("-", "").replace(".", ""))
            records.append((group, accession, genotype,
                            SeqRecord(clean, id=record.id, description="")))
    return records


def load_region_records(ref_dir, group, region):
    """Return the region-library SeqRecords for ``(group, region)``."""
    path = resolve_path(ref_dir, REFERENCE_FILES[(group, region)])
    records = []
    for record in SeqIO.parse(path, "fasta"):
        clean = Seq(str(record.seq).replace("-", "").replace(".", ""))
        records.append(SeqRecord(clean, id=record.id, description=""))
    return records


def derive_coordinates(ref_dir, cache_dir, email, api_key, retries, batch_size, force):
    """Build the full coordinate table.

    The function returns a list of dict rows matching ``COORDINATE_COLUMNS``.
    """
    if force and cache_dir.exists():
        shutil.rmtree(cache_dir)

    genome_records = load_genome_records(ref_dir)
    accessions = sorted({accession for _, accession, _, _ in genome_records})
    print(f"Genome references: {len(genome_records)} sequences, "
          f"{len(accessions)} unique accessions")

    genbank_text = fetch_genbank(accessions, cache_dir, email, api_key,
                                 retries, batch_size)
    print(f"GenBank records fetched/cached: {len(genbank_text)} / {len(accessions)}")

    # Pre-load region libraries once.
    region_libs = {}
    for group in ("GI", "GII"):
        for region in ("rdrp", "vp1"):
            region_libs[(group, region)] = load_region_records(ref_dir, group, region)

    # Split genomes by genogroup so each BLAST stays within the right region
    # library (a GI genome should not be typed against GII references).
    by_group = {"GI": [], "GII": []}
    for entry in genome_records:
        by_group[entry[0]].append(entry)

    rows = []
    work_dir = Path(tempfile.mkdtemp(prefix="noro_coords_"))
    try:
        # Derive RdRp and VP1 coordinates with one batched BLAST per
        # (genogroup, region) pair instead of one pair per genome.
        rdrp_spans = {}
        vp1_spans = {}
        for group in ("GI", "GII"):
            group_genomes = by_group[group]
            if not group_genomes:
                continue
            print(f"  BLAST {group} genomes vs RdRp region library "
                  f"({len(group_genomes)} genomes)")
            rdrp_spans.update(region_spans_for_genomes(
                group_genomes, region_libs[(group, "rdrp")], work_dir,
                f"{group}_rdrp",
            ))
            print(f"  BLAST {group} genomes vs VP1 region library "
                  f"({len(group_genomes)} genomes)")
            vp1_spans.update(region_spans_for_genomes(
                group_genomes, region_libs[(group, "vp1")], work_dir,
                f"{group}_vp1",
            ))

        for index, (group, accession, genotype, genome_record) in enumerate(genome_records, 1):
            length = len(genome_record.seq)
            if index % 200 == 0 or index == len(genome_records):
                print(f"  assembling {index}/{len(genome_records)}: {accession}")
            coords = {}

            gb_text = genbank_text.get(accession)
            if gb_text:
                try:
                    coords = parse_genbank_cds(gb_text)
                except Exception as exc:  # noqa: BLE001
                    sys.stderr.write(
                        f"  WARN: could not parse GenBank for {accession}: {exc}\n"
                    )
                    coords = {}

            # VP1 from GenBank when available; otherwise use the batched span.
            if "VP1" in coords:
                rows.append(_coord_row(group, accession, genotype, length,
                                       "VP1", coords["VP1"], "genbank"))
            elif accession in vp1_spans:
                rows.append(_coord_row(group, accession, genotype, length, "VP1",
                                       vp1_spans[accession],
                                       "inferred_from_region_blast"))

            # VP2 from GenBank when present (minor structural protein).
            if "VP2" in coords:
                rows.append(_coord_row(group, accession, genotype, length,
                                       "VP2", coords["VP2"], "genbank"))

            # ORF1 polyprotein span from GenBank when present (helps label
            # off-region queries that fall elsewhere in ORF1).
            orf1_span = coords.get("ORF1")
            if orf1_span:
                rows.append(_coord_row(group, accession, genotype, length,
                                       "ORF1", orf1_span, "genbank"))
                # Carve the ORF1 polyprotein into its six non-structural cleavage
                # products using conserved relative boundaries. This lets
                # off-region queries be annotated as NTPase / VPg / 3CLpro etc.
                for protein, frac_lo, frac_hi in ORF1_PROTEIN_BOUNDARIES:
                    if protein == "RdRp":
                        # RdRp is reported from the region library below for
                        # higher accuracy; skip the rough ORF1-slice version.
                        continue
                    orf1_start, orf1_end, _strand = orf1_span
                    orf1_len = orf1_end - orf1_start + 1
                    p_start = orf1_start + int(round(frac_lo * orf1_len))
                    p_end = orf1_start + int(round(frac_hi * orf1_len)) - 1
                    rows.append(_coord_row(group, accession, genotype, length,
                                           protein, (p_start, p_end, "+"),
                                           "orf1_relative"))

            # UTRs from GenBank when annotated.
            for utr_label in ("5'UTR", "3'UTR"):
                if utr_label in coords:
                    rows.append(_coord_row(group, accession, genotype, length,
                                           utr_label, coords[utr_label], "genbank"))

            # RdRp is essentially never an independent CDS in norovirus GenBank
            # records, so always take the region-library span.
            if accession in rdrp_spans:
                rows.append(_coord_row(group, accession, genotype, length, "RdRp",
                                       rdrp_spans[accession],
                                       "inferred_from_region_blast"))
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)

    return rows


def _coord_row(group, accession, genotype, length, protein, coord_tuple, source):
    """Build a coordinate row from ``(start, end, identity_or_strand, ...)`` tuples.

    GenBank-derived tuples carry ``(start, end, strand)`` with strand as ``+/-``;
    region-BLAST tuples carry ``(start, end, identity_pct)``. We normalise both
    to ``+`` strand (BLAST spans are already reported on the forward strand).
    """
    start, end, extra = coord_tuple
    if isinstance(extra, str):
        strand = extra
    else:
        strand = "+"  # BLAST subject spans are forward-strand by construction
    return {
        "group": group,
        "accession": accession,
        "genotype": genotype,
        "sequence_length": length,
        "protein": protein,
        "start": start,
        "end": end,
        "strand": strand,
        "source": source,
    }


def write_coordinates(rows, output_path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COORDINATE_COLUMNS,
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Download GenBank records for genome references and derive "
                    "RdRp/VP1/VP2/ORF1 coordinates."
    )
    parser.add_argument("--ref-dir", type=Path, default=Path("ref_seq"))
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"),
                        help="NCBI contact email (or set NCBI_EMAIL)")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                        help="NCBI API key for higher rate limits (or set NCBI_API_KEY)")
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--cache-dir", type=Path, default=None,
                        help="Directory for cached .gb files "
                             "(default: <ref_dir>/genbank_cache)")
    parser.add_argument("--force", action="store_true",
                        help="Re-download even if cache entries exist")
    parser.add_argument("--output", type=Path, default=None,
                        help="Output TSV path "
                             "(default: <ref_dir>/genome_gene_coordinates.tsv)")
    args = parser.parse_args()

    if not args.email:
        parser.error("NCBI requires a contact email: pass --email or set NCBI_EMAIL")

    cache_dir = args.cache_dir or (args.ref_dir / "genbank_cache")
    output_path = args.output or (args.ref_dir / "genome_gene_coordinates.tsv")

    if shutil.which("makeblastdb") is None or shutil.which("blastn") is None:
        raise RuntimeError("makeblastdb/blastn not found on PATH")

    rows = derive_coordinates(args.ref_dir, cache_dir, args.email, args.api_key,
                              args.retries, args.batch_size, args.force)
    write_coordinates(rows, output_path)

    by_protein = {}
    for row in rows:
        by_protein.setdefault(row["protein"], 0)
        by_protein[row["protein"]] += 1
    by_source = {}
    for row in rows:
        by_source.setdefault(row["source"], 0)
        by_source[row["source"]] += 1
    print(f"Wrote {len(rows)} coordinate rows to {output_path}")
    print("By protein: " + ", ".join(f"{k}={v}" for k, v in sorted(by_protein.items())))
    print("By source:  " + ", ".join(f"{k}={v}" for k, v in sorted(by_source.items())))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
