#!/usr/bin/env python3
"""Single-pass norovirus GI/GII genotyping with RdRp/VP1 phylogenetic typing.

The pipeline runs BLAST once against the near-complete genome reference
library. The top-scoring genome hit fixes the genogroup and provides, via the
precomputed gene coordinate table (``genome_gene_coordinates.tsv``), the RdRp /
VP1 coordinates needed to slice the query. The slices are placed into the
per-genogroup RdRp and VP1 reference trees; a genotype is reported only when the
SH-like branch support at the placement clade clears the bootstrap threshold.

Sequences that do not cover RdRp or VP1 are annotated with their genomic
position (ORF1 / ORF2 / ORF3 and the overlapping protein) and a list of the
closest genome references, instead of being forced through the tree step.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


REFERENCE_FILES = {
    ("GI", "genome"): ("gi_with_genotype.fasta", "gi_genome.fasta"),
    ("GI", "rdrp"): "gi_rdrp.fasta",
    ("GI", "vp1"): ("gi_vp1.fasta", "gi_VP1.fa"),
    ("GII", "genome"): ("gii_with_genotype.fasta", "gii_genome.fasta"),
    ("GII", "rdrp"): "gii_rdrp.fasta",
    ("GII", "vp1"): "gii_vp1.fasta",
}

COORDINATES_FILENAME = "genome_gene_coordinates.tsv"

# SH-like branch support below which a genotype call is reported as "unsure".
BOOTSTRAP_THRESHOLD = 0.75

# Proteins used to label the genomic position of a query hit. RdRp/VP1 are the
# typed regions; the rest are reported for off-region annotation only.
TYPED_PROTEINS = {"RdRp", "VP1"}

# ORF1 is cleaved into non-structural proteins. When a query overlaps one of
# these fine-grained proteins we report the specific name (e.g. "NTPase") rather
# than the generic "ORF1". ORF1 itself is kept as a fallback label.
ORF1_GRANULAR = {"p48", "NTPase", "p22", "VPg", "3CLpro", "RdRp"}
ORF1_PROTEINS = ORF1_GRANULAR | {"ORF1"}
ORF2_PROTEINS = {"VP1"}
ORF3_PROTEINS = {"VP2"}
UTR_PROTEINS = {"5'UTR", "3'UTR"}

# Fine-grained proteins that subsume the generic "ORF1" label. When both are
# reported we drop the generic one to avoid noise in the annotation.
ORF1_GENERIC = {"ORF1"}

PROTEIN_TO_REGION = {"RdRp": "rdrp", "VP1": "vp1"}

OUTPUT_FIELDS = [
    "input_file",
    "sequence_name",
    "sequence_length",
    "status",
    "genogroup",
    # Stage 1: genome BLAST coarse classification
    "coarse_reference",
    "coarse_identity_pct",
    "coarse_alignment_length",
    "coarse_query_coverage_pct",
    "coarse_bitscore",
    # Genomic position annotation for every query (driven by gene coordinates)
    "query_genome_region",
    "query_genome_proteins",
    "query_genome_coords",
    # Stage 2: RdRp phylogenetic placement
    "rdrp_nearest_reference",
    "rdrp_genotype",
    "rdrp_genotype_confidence",
    "rdrp_bootstrap",
    "rdrp_tree_distance",
    "rdrp_blast_identity_pct",
    "rdrp_blast_alignment_length",
    "rdrp_reference_coverage_pct",
    # Stage 2: VP1 phylogenetic placement
    "vp1_nearest_reference",
    "vp1_genotype",
    "vp1_genotype_confidence",
    "vp1_bootstrap",
    "vp1_tree_distance",
    "vp1_blast_identity_pct",
    "vp1_blast_alignment_length",
    "vp1_reference_coverage_pct",
    # Off-region annotation: closest genome references
    "off_region_top_references",
]


def run(command, stdout=None, env=None):
    try:
        subprocess.run(
            [str(item) for item in command],
            check=True,
            stdout=stdout,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
    except subprocess.CalledProcessError as exc:
        detail = exc.stderr.strip() if exc.stderr else "no error message"
        raise RuntimeError(f"Command failed: {' '.join(map(str, command))}\n{detail}") from exc


def require_tools(aligner):
    missing = [
        tool
        for tool in ("makeblastdb", "blastn", aligner, "fasttree")
        if shutil.which(tool) is None
    ]
    if missing:
        raise RuntimeError("Missing required tools: " + ", ".join(missing))


def genotype_from_header(header, expected_group):
    token = header.split()[0]
    match = re.search(r"(GI|GII)\.(P?(?:NA\d+|\d+))", token, re.IGNORECASE)
    if match:
        return f"{match.group(1).upper()}.{match.group(2).upper()}"

    # GI VP1 references use headers such as I.1|M87661.
    match = re.match(r"I\.(NA\d+|\d+)\|", token, re.IGNORECASE)
    if match and expected_group == "GI":
        return f"GI.{match.group(1).upper()}"

    return "Unknown"


ACCESSION_RE = re.compile(r"(NC_\d+|[A-Z]{1,2}\d{5,8})", re.IGNORECASE)


def accession_from_header(header):
    """Return the NCBI accession embedded in a genome FASTA header.

    Genome headers look like ``GI.1_MH638228`` (genotype_accession); we fall
    back to the first whitespace-delimited token if no accession pattern matches.
    """
    match = ACCESSION_RE.search(header)
    return match.group(1).upper() if match else header.split()[0]


def resolve_reference_path(ref_dir, filenames):
    if isinstance(filenames, str):
        filenames = (filenames,)
    for filename in filenames:
        path = ref_dir / filename
        if path.is_file():
            return path
    expected = ", ".join(str(ref_dir / filename) for filename in filenames)
    raise FileNotFoundError(f"Reference file not found; expected one of: {expected}")


def load_references(ref_dir):
    references = {}
    metadata = {}
    counter = 0

    for (group, region), filenames in REFERENCE_FILES.items():
        path = resolve_reference_path(ref_dir, filenames)

        records = []
        for record in SeqIO.parse(path, "fasta"):
            counter += 1
            normalized_id = f"R{counter:05d}"
            clean_sequence = Seq(str(record.seq).replace("-", "").replace(".", ""))
            normalized = SeqRecord(clean_sequence, id=normalized_id, description="")
            records.append(normalized)
            metadata[normalized_id] = {
                "name": record.id,
                "group": group,
                "region": region,
                "genotype": (
                    genotype_from_header(record.id, group)
                    if region != "genome"
                    else ""
                ),
                "accession": (
                    accession_from_header(record.id)
                    if region == "genome"
                    else ""
                ),
                "length": len(clean_sequence),
            }
        if not records:
            raise ValueError(f"No reference sequences found in {path}")
        references[(group, region)] = records

    return references, metadata


def load_gene_coordinates(ref_dir):
    """Read ``genome_gene_coordinates.tsv`` into a nested lookup.

    Returns ``{accession: {protein: (start, end)}}`` with 1-based inclusive
    coordinates on the forward genome strand. Proteins without a table entry are
    simply absent from the inner dict. Returns ``{}`` if the file is missing so
    callers can degrade gracefully (and warn the user).
    """
    path = ref_dir / COORDINATES_FILENAME
    if not path.is_file():
        return {}
    coordinates = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            accession = row["accession"].strip().upper()
            protein = row["protein"].strip()
            try:
                start = int(row["start"])
                end = int(row["end"])
            except (ValueError, KeyError):
                continue
            coordinates.setdefault(accession, {})[protein] = (start, end)
    return coordinates


def load_queries(input_path):
    records = []
    metadata = {}
    for index, record in enumerate(SeqIO.parse(input_path, "fasta"), start=1):
        query_id = f"Q{index:06d}"
        clean_sequence = Seq(str(record.seq).replace("-", "").replace(".", ""))
        records.append(SeqRecord(clean_sequence, id=query_id, description=""))
        metadata[query_id] = {
            "name": record.id,
            "description": record.description,
            "length": len(clean_sequence),
            "sequence": clean_sequence,
        }
    if not records:
        raise ValueError(f"No sequences found in {input_path}")
    return records, metadata


def run_blast(queries, references, work_dir, threads, off_region_top_n):
    """BLAST queries once against the genome reference library.

    Only the genome references are searched: the top-scoring hit fixes the
    genogroup and (via the gene coordinate table) the RdRp/VP1 coordinates used
    to slice each query, which removes the previous second-pass region BLAST.

    ``-max_target_seqs`` is set high enough to also surface the next-best genome
    hits so off-region queries can be annotated with their closest references.
    Returns a list of hit dicts in the same shape as before.
    """
    query_fasta = work_dir / "queries.fasta"
    SeqIO.write(queries, query_fasta, "fasta")

    records = [
        record
        for (group, region), ref_records in references.items()
        if region == "genome"
        for record in ref_records
    ]
    reference_fasta = work_dir / "genome_references.fasta"
    database = work_dir / "genome_refs"
    output = work_dir / "genome_hits.tsv"
    SeqIO.write(records, reference_fasta, "fasta")
    run(
        [
            "makeblastdb",
            "-in",
            reference_fasta,
            "-dbtype",
            "nucl",
            "-out",
            database,
        ]
    )
    run(
        [
            "blastn",
            "-task",
            "blastn",
            "-query",
            query_fasta,
            "-db",
            database,
            "-num_threads",
            threads,
            "-max_target_seqs",
            str(max(len(records), off_region_top_n + 5)),
            "-max_hsps",
            "50",
            "-evalue",
            "1e-10",
            "-outfmt",
            "6 qseqid sseqid pident length qstart qend sstart send sstrand bitscore evalue",
            "-out",
            output,
        ]
    )
    hits = []
    with output.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            hits.append(
                {
                    "query": row[0],
                    "reference": row[1],
                    "identity": float(row[2]),
                    "alignment_length": int(row[3]),
                    "qstart": int(row[4]),
                    "qend": int(row[5]),
                    "sstart": int(row[6]),
                    "send": int(row[7]),
                    "sstrand": row[8],
                    "bitscore": float(row[9]),
                    "evalue": float(row[10]),
                }
            )
    return hits


def query_span(hit):
    return abs(hit["qend"] - hit["qstart"]) + 1


def reference_span(hit):
    return abs(hit["send"] - hit["sstart"]) + 1


def merged_interval_length(intervals):
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start + 1 for start, end in merged)


def classify_queries(
    query_metadata,
    hits,
    reference_metadata,
    min_identity,
    min_query_coverage,
    off_region_top_n,
):
    """Assign each query to a genogroup from its genome BLAST hits.

    Returns ``(assignments, top_genome_hits)`` where ``assignments`` maps
    ``query_id -> (group, best_hit)`` (best hit aggregated across HSPs) and
    ``top_genome_hits`` maps ``query_id -> [hit, ...]`` holding the next-best
    genome references for off-region annotation.
    """
    assignments = {}
    top_genome_hits = {}
    for query_id in query_metadata:
        grouped = {}
        for hit in hits:
            if hit["query"] != query_id:
                continue
            ref = reference_metadata[hit["reference"]]
            if ref["region"] != "genome":
                continue
            key = hit["reference"]
            grouped.setdefault(key, []).append(hit)

        candidates = []
        for reference_id, reference_hits in grouped.items():
            total_alignment = sum(hit["alignment_length"] for hit in reference_hits)
            if total_alignment == 0:
                continue
            covered = merged_interval_length(
                [
                    (min(hit["qstart"], hit["qend"]), max(hit["qstart"], hit["qend"]))
                    for hit in reference_hits
                ]
            )
            coverage = 100.0 * covered / query_metadata[query_id]["length"]
            identity = (
                sum(hit["identity"] * hit["alignment_length"] for hit in reference_hits)
                / total_alignment
            )
            bitscore = sum(hit["bitscore"] for hit in reference_hits)
            best_hsp = max(
                reference_hits,
                key=lambda item: (
                    item["bitscore"],
                    item["alignment_length"],
                    item["identity"],
                ),
            ).copy()
            best_hsp.update(
                {
                    "identity": identity,
                    "alignment_length": total_alignment,
                    "bitscore": bitscore,
                    "query_coverage_pct": coverage,
                }
            )
            if identity >= min_identity and coverage >= min_query_coverage:
                candidates.append(best_hsp)

        if not candidates:
            assignments[query_id] = (None, None)
        else:
            best = max(
                candidates,
                key=lambda item: (
                    item["bitscore"],
                    item["alignment_length"],
                    item["identity"],
                ),
            )
            assignments[query_id] = (
                reference_metadata[best["reference"]]["group"],
                best,
            )

        # Always collect the top-N genome references by aggregated bitscore for
        # off-region annotation, regardless of whether the query classifies.
        ranked = sorted(
            (
                _aggregate_reference(reference_id, reference_hits)
                for reference_id, reference_hits in grouped.items()
            ),
            key=lambda item: (item["bitscore"], item["alignment_length"], item["identity"]),
            reverse=True,
        )
        top_genome_hits[query_id] = ranked[:off_region_top_n]

    return assignments, top_genome_hits


def _aggregate_reference(reference_id, reference_hits):
    """Collapse the HSPs against one reference into a single summary hit."""
    total_alignment = sum(hit["alignment_length"] for hit in reference_hits)
    best_hsp = max(
        reference_hits,
        key=lambda item: (
            item["bitscore"],
            item["alignment_length"],
            item["identity"],
        ),
    ).copy()
    best_hsp.update(
        {
            "reference": reference_id,
            "alignment_length": total_alignment,
            "bitscore": sum(hit["bitscore"] for hit in reference_hits),
        }
    )
    return best_hsp


def slice_query_for_protein(query_seq, hit, ref_coords, protein):
    """Slice the query region corresponding to ``protein`` using the best hit.

    The best genome BLAST hit aligns query ``[qstart, qend]`` to the reference
    span ``[sstart, send]``. ``ref_coords`` gives the protein's 1-based inclusive
    coordinates ``(p_start, p_end)`` on the reference genome. Assuming approximate
    colinearity we project those reference coordinates back onto the query:

        plus  strand:  q = qstart + (p - sstart)
        minus strand:  q = qend   - (p - sstart)   (then reverse-complement)

    Returns a ``Seq`` in the protein's forward orientation, or ``None`` when the
    protein's coordinates fall outside the aligned span (nothing to slice) or no
    usable slice survives the boundary clamping.
    """
    if not ref_coords or protein not in ref_coords:
        return None
    p_start, p_end = ref_coords[protein]
    s_lo, s_hi = min(hit["sstart"], hit["send"]), max(hit["sstart"], hit["send"])

    # Clamp the requested protein window to the part that the HSP actually
    # covers; if there is no overlap at all, the query does not span this region.
    overlap_lo = max(p_start, s_lo)
    overlap_hi = min(p_end, s_hi)
    if overlap_hi < overlap_lo:
        return None

    if hit["sstrand"] == "minus":
        q_anchor = max(hit["qstart"], hit["qend"])
        q_lo = q_anchor - (overlap_hi - s_lo)
        q_hi = q_anchor - (overlap_lo - s_lo)
    else:
        q_anchor = min(hit["qstart"], hit["qend"])
        q_lo = q_anchor + (overlap_lo - s_lo)
        q_hi = q_anchor + (overlap_hi - s_lo)

    q_lo = max(1, q_lo)
    q_hi = min(len(query_seq), q_hi)
    if q_hi < q_lo:
        return None
    segment = query_seq[q_lo - 1:q_hi]
    if hit["sstrand"] == "minus":
        segment = segment.reverse_complement()
    return segment


def map_hit_to_genome_location(hit, ref_coords):
    """Project a query's reference span onto the gene coordinate table.

    Returns ``{"coords": (gstart, gend), "proteins": [...], "region": str}``.
    ``proteins`` lists the proteins overlapping the span, preferring the most
    specific label (e.g. ``NTPase`` over the generic ``ORF1``) and ordering the
    typed regions first. ``region`` summarises the ORF context
    (``ORF1``/``ORF2``/``ORF3``/``5'UTR``/``3'UTR``/``intergenic``/``multi``).
    """
    s_lo = min(hit["sstart"], hit["send"])
    s_hi = max(hit["sstart"], hit["send"])
    overlaps = []
    for protein, (p_start, p_end) in ref_coords.items():
        if p_end < s_lo or p_start > s_hi:
            continue
        overlap = min(p_end, s_hi) - max(p_start, s_lo) + 1
        overlaps.append((protein, overlap))
    # Drop the generic "ORF1" label whenever a fine-grained ORF1 protein also
    # overlaps, so the annotation reads "NTPase" instead of "NTPase, ORF1".
    has_granular = any(protein in ORF1_GRANULAR for protein, _ in overlaps)
    if has_granular:
        overlaps = [(p, o) for p, o in overlaps if p not in ORF1_GENERIC]
    # Typed proteins first, then by descending overlap.
    overlaps.sort(key=lambda item: (item[0] not in TYPED_PROTEINS, -item[1]))

    orfs = set()
    for protein, _ in overlaps:
        if protein in ORF1_PROTEINS:
            orfs.add("ORF1")
        elif protein in ORF2_PROTEINS:
            orfs.add("ORF2")
        elif protein in ORF3_PROTEINS:
            orfs.add("ORF3")
        elif protein in UTR_PROTEINS:
            orfs.add(protein)
    if len(orfs) > 1:
        region = "multi"
    elif orfs:
        region = next(iter(orfs))
    else:
        region = "intergenic"

    return {
        "coords": (s_lo, s_hi),
        "proteins": [protein for protein, _ in overlaps] or ["intergenic"],
        "region": region,
    }


def assign_one_query_by_tree(task):
    (
        query_id,
        query_sequence,
        region,
        reference_alignment,
        reference_ids,
        reference_metadata,
        work_dir,
        tree_dir,
    ) = task

    query_work_dir = work_dir / f"{region}_{query_id}"
    query_work_dir.mkdir(parents=True, exist_ok=True)
    query_fasta = query_work_dir / "query_segment.fasta"
    combined_alignment = query_work_dir / "combined_aligned.fasta"
    tree_path = tree_dir / f"{region}_{query_id}.nwk"

    SeqIO.write([SeqRecord(query_sequence, id=query_id, description="")], query_fasta, "fasta")

    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")

    with combined_alignment.open("w") as handle:
        run(
            [
                "mafft",
                "--quiet",
                "--addfragments",
                query_fasta,
                "--reorder",
                "--thread",
                "1",
                reference_alignment,
            ],
            stdout=handle,
            env=env,
        )
    # -nome -mllen makes FastTree emit SH-like branch support as the internal
    # node labels, which we read back as the placement confidence.
    with tree_path.open("w") as handle:
        run(
            ["fasttree", "-nt", "-nome", "-mllen", combined_alignment],
            stdout=handle,
            env=env,
        )

    tree = Phylo.read(tree_path, "newick")
    nearest_id = None
    nearest_distance = float("inf")
    for reference_id in reference_ids:
        distance = tree.distance(query_id, reference_id)
        if distance < nearest_distance:
            nearest_distance = distance
            nearest_id = reference_id

    bootstrap = clade_support(tree, query_id, nearest_id)

    nearest = reference_metadata[nearest_id]
    return (
        query_id,
        nearest_id,
        nearest_distance,
        nearest["genotype"],
        nearest["name"],
        bootstrap,
    )


def clade_support(tree, query_id, reference_id):
    """SH-like support (0-1) for the placement of ``query_id`` near ``reference_id``.

    FastTree writes SH-like support as internal node labels (read here as
    ``clade.confidence``). When a query is placed as a fragment it often forms a
    two-tip cherry with its nearest reference; that cherry node carries little
    information and FastTree reports 0 support for it. In that case we walk up
    to the first ancestor that contains at least one other reference and report
    *that* node's support, since it reflects how stably the query+neighbour
    group separates from the rest of the tree.

    Falls back to 0.0 when no support value can be recovered.
    """
    try:
        mrca = tree.common_ancestor(query_id, reference_id)
    except (ValueError, KeyError):
        return 0.0
    if mrca is None:
        return 0.0

    # Walk up from the placement clade until we reach a node that includes more
    # than just the query and its nearest neighbour. The first such ancestor is
    # the deepest clade that actually tests the placement.
    node = mrca
    tips_seen = {query_id, reference_id}
    while node is not None:
        tip_names = {t.name for t in node.get_terminals()}
        if len(tip_names - tips_seen) == 0:
            # Same two tips; this is the placement cherry itself — keep walking.
            node = _parent_of(tree, node)
            continue
        support = _read_support(node)
        if support is not None:
            return support
        node = _parent_of(tree, node)
    return 0.0


def _parent_of(tree, clade):
    """Return the direct parent of ``clade`` in ``tree`` or ``None`` at the root."""
    if clade is tree.root:
        return None
    for candidate in tree.get_nonterminals():
        if clade in candidate.clades:
            return candidate
    return None


def _read_support(clade):
    """Return the numeric support of ``clade`` (0-1) or ``None`` if unlabelled."""
    if clade.confidence is not None:
        try:
            return max(0.0, min(1.0, float(clade.confidence)))
        except (TypeError, ValueError):
            pass
    if clade.name:
        try:
            return max(0.0, min(1.0, float(clade.name)))
        except (TypeError, ValueError):
            return None
    return None


def build_tree_and_assign(
    group,
    region,
    query_ids,
    query_metadata,
    references,
    reference_metadata,
    gene_coordinates,
    assignments,
    work_dir,
    tree_dir,
    threads,
    min_slice_identity,
    min_slice_length,
):
    """Place each query's gene slice into the per-genogroup region tree.

    The query slice is taken from the genome coordinate table (via the best
    genome hit) rather than from a second region BLAST. ``min_slice_identity``
    and ``min_slice_length`` gate which queries are worth placing.
    """
    protein = next(
        protein for protein, reg in PROTEIN_TO_REGION.items() if reg == region
    )
    usable_queries = {}
    slice_meta = {}
    for query_id in query_ids:
        assigned_group, best_hit = assignments.get(query_id, (None, None))
        if best_hit is None:
            continue
        ref = reference_metadata[best_hit["reference"]]
        ref_coords = gene_coordinates.get(ref["accession"], {})
        if protein not in ref_coords:
            continue
        if best_hit["identity"] < min_slice_identity:
            continue
        segment = slice_query_for_protein(
            query_metadata[query_id]["sequence"], best_hit, ref_coords, protein
        )
        if segment is None or len(segment) < min_slice_length:
            continue
        usable_queries[query_id] = segment
        slice_meta[query_id] = {
            "best_hit": best_hit,
            "slice_length": len(segment),
            "protein_ref_length": abs(ref_coords[protein][1] - ref_coords[protein][0]) + 1,
        }

    if not usable_queries:
        return {}

    prefix = f"{group}_{region}"
    region_work_dir = work_dir / prefix
    region_work_dir.mkdir(parents=True, exist_ok=True)
    region_tree_dir = tree_dir / prefix
    region_tree_dir.mkdir(parents=True, exist_ok=True)
    reference_fasta = region_work_dir / f"{prefix}_references.fasta"
    reference_alignment = region_work_dir / f"{prefix}_references_aligned.fasta"

    SeqIO.write(references[(group, region)], reference_fasta, "fasta")
    with reference_alignment.open("w") as handle:
        run(
            [
                "mafft",
                "--quiet",
                "--auto",
                "--thread",
                threads,
                reference_fasta,
            ],
            stdout=handle,
        )

    reference_ids = {record.id for record in references[(group, region)]}
    results = {}
    tasks = [
        (
            query_id,
            segment,
            prefix,
            reference_alignment,
            reference_ids,
            reference_metadata,
            region_work_dir,
            region_tree_dir,
        )
        for query_id, segment in usable_queries.items()
    ]
    max_workers = max(1, int(threads))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_query = {
            executor.submit(assign_one_query_by_tree, task): task[0]
            for task in tasks
        }
        for future in as_completed(future_to_query):
            query_id = future_to_query[future]
            result = future.result()
            nearest_id = result[1]
            nearest_distance = result[2]
            bootstrap = result[5]
            results[query_id] = {
                "reference_id": nearest_id,
                "distance": nearest_distance,
                "bootstrap": bootstrap,
                "slice": slice_meta[query_id],
            }
    return results


def derive_status(group, coarse_hit, typed_regions, confident_regions, ref_coords):
    """Pick the final status label for a query.

    - ``unclassified``    : no genogroup (no usable genome hit)
    - ``off_region``      : genogroup assigned but the hit does not overlap RdRp
                            or VP1 (no typed region was sliced)
    - ``genotype_unsure`` : at least one region was placed but none cleared the
                            bootstrap threshold
    - ``partial_genotype``: exactly one of the two regions cleared bootstrap
    - ``typed``           : both regions cleared bootstrap
    - ``genogroup_only``  : genogroup assigned, region coordinates unavailable
                            or too short to place, yet the hit does overlap a
                            typed region (coordinate table missing, etc.)
    """
    if not group:
        return "unclassified"
    if not coarse_hit:
        return "genogroup_only"
    # Determine whether the genome hit overlapped a typed region at all.
    s_lo = min(coarse_hit["sstart"], coarse_hit["send"])
    s_hi = max(coarse_hit["sstart"], coarse_hit["send"])
    overlaps_typed = any(
        protein in ref_coords
        and not (coords[1] < s_lo or coords[0] > s_hi)
        for protein, coords in ref_coords.items()
        if protein in TYPED_PROTEINS
    )
    if not typed_regions:
        return "off_region" if not overlaps_typed else "genogroup_only"
    if not confident_regions:
        return "genotype_unsure"
    if len(confident_regions) == 1:
        return "partial_genotype"
    return "typed"


def build_fasta_annotation(row, typed_regions, confident_regions):
    """Compose the FASTA description line for an annotated query."""
    parts = [f"genogroup={row['genogroup']}"]
    for region in ("rdrp", "vp1"):
        genotype = row.get(f"{region}_genotype", "")
        if not genotype:
            continue
        confidence = row.get(f"{region}_genotype_confidence", "")
        label = f"{region.upper()}={genotype}"
        if confidence == "unsure":
            label += "(unsure)"
        parts.append(label)
    if row.get("query_genome_region"):
        parts.append(f"region={row['query_genome_region']}")
        proteins = row.get("query_genome_proteins", "")
        if proteins and proteins != "intergenic":
            parts.append(f"proteins={proteins}")
    return " ".join(parts)


def format_number(value, decimals=3):
    if value is None or value == "":
        return ""
    return f"{value:.{decimals}f}"


def write_html_report(rows, tree_dir, reference_metadata, report_path,
                      input_name, bootstrap_threshold):
    """Render the self-contained HTML report via ``generate_report``.

    Imported lazily so the main pipeline does not pay the import cost when the
    ``--report`` flag is not used, and so the report module can be developed or
    shipped independently.
    """
    import importlib

    try:
        report_module = importlib.import_module("generate_report")
    except ImportError as exc:
        print(f"WARNING: --report requested but generate_report module is "
              f"unavailable: {exc}", file=sys.stderr)
        return
    # Map normalised reference IDs (R00001..) to readable names for the tree
    # tip labels.
    reference_names = {rid: meta["name"] for rid, meta in reference_metadata.items()}
    # generate_report.render_report expects plain dicts; our rows already are.
    report_module.BOOTSTRAP_THRESHOLD = bootstrap_threshold
    report_module.render_report(
        rows, tree_dir, reference_names, report_path,
        title="Norovirus Genotyping Report",
        input_name=input_name,
    )
    print(f"Wrote HTML report to {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Classify norovirus GI/GII and assign RdRp and VP1 genotypes by phylogeny."
    )
    parser.add_argument("--input", required=True, type=Path, help="Input FASTA")
    parser.add_argument("--ref-dir", required=True, type=Path, help="Reference directory")
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--aligner",
        choices=("mafft",),
        default="mafft",
        help="Multiple sequence aligner used before FastTree; each query is aligned separately",
    )
    parser.add_argument("--min-identity", type=float, default=70.0)
    parser.add_argument(
        "--min-query-coverage",
        type=float,
        default=50.0,
        help="Minimum input-sequence coverage for GI/GII genome classification",
    )
    parser.add_argument(
        "--min-slice-identity",
        type=float,
        default=70.0,
        help="Minimum genome-hit identity for a query to enter the RdRp/VP1 tree step",
    )
    parser.add_argument(
        "--min-slice-length",
        type=int,
        default=300,
        help="Minimum sliced segment length (nt) for a query to enter the RdRp/VP1 tree step",
    )
    parser.add_argument(
        "--bootstrap-threshold",
        type=float,
        default=BOOTSTRAP_THRESHOLD,
        help="Minimum SH-like branch support for a confident genotype call; "
             "below this the genotype is reported as 'unsure'",
    )
    parser.add_argument(
        "--off-region-top-n",
        type=int,
        default=3,
        help="Number of closest genome references listed for off-region queries",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Also write a self-contained HTML report (with phylogenetic trees "
             "and SH support values) next to the TSV output",
    )
    args = parser.parse_args()

    require_tools(args.aligner)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tree_dir = Path(f"{args.output_prefix}_trees")
    tree_dir.mkdir(parents=True, exist_ok=True)

    references, reference_metadata = load_references(args.ref_dir)
    queries, query_metadata = load_queries(args.input)
    gene_coordinates = load_gene_coordinates(args.ref_dir)
    if not gene_coordinates:
        print(
            "WARNING: gene coordinate table not found "
            f"({args.ref_dir / COORDINATES_FILENAME}); RdRp/VP1 slicing and "
            "genomic position annotation will be unavailable. "
            "Run download_genome_annotations.py to generate it.",
            file=sys.stderr,
        )

    with tempfile.TemporaryDirectory(prefix="noro_genotyping_") as temp:
        work_dir = Path(temp)
        threads_arg = str(args.threads)
        hits = run_blast(
            queries, references, work_dir, threads_arg, args.off_region_top_n
        )
        assignments, top_genome_hits = classify_queries(
            query_metadata,
            hits,
            reference_metadata,
            args.min_identity,
            args.min_query_coverage,
            args.off_region_top_n,
        )

        tree_results = {}
        if gene_coordinates:
            for group in ("GI", "GII"):
                group_queries = [
                    query_id
                    for query_id, (assigned_group, _) in assignments.items()
                    if assigned_group == group
                ]
                for region in ("rdrp", "vp1"):
                    regional_results = build_tree_and_assign(
                        group,
                        region,
                        group_queries,
                        query_metadata,
                        references,
                        reference_metadata,
                        gene_coordinates,
                        assignments,
                        work_dir,
                        tree_dir,
                        threads_arg,
                        args.min_slice_identity,
                        args.min_slice_length,
                    )
                    for query_id, result in regional_results.items():
                        tree_results[(query_id, region)] = result

    rows = []
    annotated_records = []
    input_name = args.input.name
    for query in queries:
        query_id = query.id
        original = query_metadata[query_id]
        group, coarse_hit = assignments[query_id]
        row = {field: "" for field in OUTPUT_FIELDS}
        row.update(
            {
                "input_file": input_name,
                "sequence_name": original["name"],
                "sequence_length": original["length"],
                "status": "unclassified",
                "genogroup": group or "Unknown",
            }
        )

        ref_coords = {}
        if coarse_hit:
            ref = reference_metadata[coarse_hit["reference"]]
            row.update(
                {
                    "coarse_reference": ref["name"],
                    "coarse_identity_pct": format_number(coarse_hit["identity"]),
                    "coarse_alignment_length": coarse_hit["alignment_length"],
                    "coarse_query_coverage_pct": format_number(
                        coarse_hit.get(
                            "query_coverage_pct",
                            100.0 * query_span(coarse_hit) / original["length"],
                        )
                    ),
                    "coarse_bitscore": format_number(coarse_hit["bitscore"], 1),
                }
            )
            ref_coords = gene_coordinates.get(ref["accession"], {})
            location = map_hit_to_genome_location(coarse_hit, ref_coords)
            row.update(
                {
                    "query_genome_region": location["region"],
                    "query_genome_proteins": "+".join(location["proteins"]),
                    "query_genome_coords": f"{location['coords'][0]}-{location['coords'][1]}",
                }
            )

        typed_regions = []
        confident_regions = []
        for region in ("rdrp", "vp1"):
            result = tree_results.get((query_id, region))
            if result is None:
                continue
            nearest = reference_metadata[result["reference_id"]]
            slice_info = result["slice"]
            best_hit = slice_info["best_hit"]
            slice_length = slice_info["slice_length"]
            protein_ref_length = slice_info["protein_ref_length"]
            confident = result["bootstrap"] >= args.bootstrap_threshold
            confidence_label = "typed" if confident else "unsure"
            row.update(
                {
                    f"{region}_nearest_reference": nearest["name"],
                    f"{region}_genotype": nearest["genotype"],
                    f"{region}_genotype_confidence": confidence_label,
                    f"{region}_bootstrap": format_number(result["bootstrap"], 3),
                    f"{region}_tree_distance": format_number(result["distance"], 6),
                    # Slice statistics: identity comes from the genome hit, but the
                    # alignment length and coverage describe the extracted gene
                    # segment rather than the whole-genome span.
                    f"{region}_blast_identity_pct": format_number(best_hit["identity"]),
                    f"{region}_blast_alignment_length": slice_length,
                    f"{region}_reference_coverage_pct": format_number(
                        100.0 * slice_length / protein_ref_length
                        if protein_ref_length else 0.0
                    ),
                }
            )
            typed_regions.append(region)
            if confident:
                confident_regions.append(region)

        # Off-region closest references (need at least the coarse hit).
        proteins = row.get("query_genome_proteins", "")
        typed_overlap = bool(proteins) and any(
            protein in TYPED_PROTEINS for protein in proteins.split("+")
        )
        if coarse_hit and not typed_overlap:
            top_refs = top_genome_hits.get(query_id, [])
            row["off_region_top_references"] = ",".join(
                f"{reference_metadata[hit['reference']]['name']}({hit['identity']:.1f}%)"
                for hit in top_refs
            )

        row["status"] = derive_status(
            group, coarse_hit, typed_regions, confident_regions, ref_coords
        )

        annotation = build_fasta_annotation(row, typed_regions, confident_regions)
        annotated_records.append(
            SeqRecord(
                original["sequence"],
                id=original["name"],
                description=annotation,
            )
        )
        rows.append(row)

    table_path = Path(f"{args.output_prefix}_genotyping.tsv")
    fasta_path = Path(f"{args.output_prefix}_genotyped.fasta")
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    SeqIO.write(annotated_records, fasta_path, "fasta")

    print(f"Wrote {len(rows)} sequence records to {table_path}")
    print(f"Wrote annotated FASTA to {fasta_path}")
    print(f"Wrote phylogenetic trees to {tree_dir}")

    if args.report:
        report_path = Path(f"{args.output_prefix}_report.html")
        write_html_report(
            rows, tree_dir, reference_metadata, report_path,
            args.input.name, args.bootstrap_threshold,
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
