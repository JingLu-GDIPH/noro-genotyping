#!/usr/bin/env python3
"""Supplement genotype-labelled genome FASTA files with missing VP1 types."""

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from genotype_norovirus import genotype_from_header


FILES = {
    "GI": {
        "original": "gi_with_genotype-original.fasta",
        "genome": "gi_genome.fasta",
        "vp1": ("gi_VP1.fa", "gi_vp1.fasta"),
        "output": "gi_with_genotype.fasta",
    },
    "GII": {
        "original": "gii_with_genotype-original.fasta",
        "genome": "gii_genome.fasta",
        "vp1": ("gii_vp1.fasta",),
        "output": "gii_with_genotype.fasta",
    },
}


def clean_sequence(sequence):
    return Seq(str(sequence).upper().replace("-", "").replace(".", ""))


def original_genotype(identifier, group):
    match = re.match(
        rf"({group})\.(NA\d+|\d+)_",
        identifier,
        re.IGNORECASE,
    )
    if not match:
        raise ValueError(f"Cannot parse genotype from original header: {identifier}")
    return f"{match.group(1).upper()}.{match.group(2).upper()}"


def accession_from_identifier(identifier):
    match = re.search(r"(NC_\d+|[A-Z]{1,2}\d{5,8})", identifier, re.IGNORECASE)
    return match.group(1).upper() if match else identifier.split()[0]


def resolve_path(ref_dir, names):
    for name in names:
        path = ref_dir / name
        if path.is_file():
            return path
    raise FileNotFoundError(
        "Reference file not found; expected one of: "
        + ", ".join(str(ref_dir / name) for name in names)
    )


def run_blast(vp1_path, genome_path, group, work_dir):
    vp1_records = []
    vp1_metadata = {}
    for index, record in enumerate(SeqIO.parse(vp1_path, "fasta"), start=1):
        reference_id = f"R{index:04d}"
        genotype = genotype_from_header(record.id, group)
        if genotype == "Unknown":
            continue
        vp1_records.append(
            SeqRecord(clean_sequence(record.seq), id=reference_id, description="")
        )
        vp1_metadata[reference_id] = {
            "name": record.id,
            "genotype": genotype,
        }

    genome_records = []
    genome_metadata = {}
    for index, record in enumerate(SeqIO.parse(genome_path, "fasta"), start=1):
        query_id = f"Q{index:04d}"
        genome_records.append(
            SeqRecord(clean_sequence(record.seq), id=query_id, description="")
        )
        genome_metadata[query_id] = record

    vp1_fasta = work_dir / f"{group}_vp1.fasta"
    genome_fasta = work_dir / f"{group}_genomes.fasta"
    database = work_dir / f"{group}_vp1_db"
    hit_path = work_dir / f"{group}_hits.tsv"
    SeqIO.write(vp1_records, vp1_fasta, "fasta")
    SeqIO.write(genome_records, genome_fasta, "fasta")

    subprocess.run(
        ["makeblastdb", "-in", vp1_fasta, "-dbtype", "nucl", "-out", database],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    subprocess.run(
        [
            "blastn",
            "-query",
            genome_fasta,
            "-db",
            database,
            "-max_target_seqs",
            "20",
            "-max_hsps",
            "1",
            "-evalue",
            "1e-20",
            "-outfmt",
            "6 qseqid sseqid pident length qstart qend sstart send bitscore",
            "-out",
            hit_path,
        ],
        check=True,
        stderr=subprocess.PIPE,
        text=True,
    )

    best = {}
    with hit_path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            query_id, reference_id = row[:2]
            hit = {
                "query_id": query_id,
                "reference_id": reference_id,
                "identity": float(row[2]),
                "alignment_length": int(row[3]),
                "qstart": int(row[4]),
                "qend": int(row[5]),
                "sstart": int(row[6]),
                "send": int(row[7]),
                "bitscore": float(row[8]),
            }
            rank = (hit["bitscore"], hit["alignment_length"], hit["identity"])
            if query_id not in best or rank > best[query_id]["rank"]:
                hit["rank"] = rank
                best[query_id] = hit

    typed = []
    for query_id, hit in best.items():
        vp1 = vp1_metadata[hit["reference_id"]]
        typed.append(
            {
                "record": genome_metadata[query_id],
                "genotype": vp1["genotype"],
                "nearest_vp1_reference": vp1["name"],
                "identity": hit["identity"],
                "alignment_length": hit["alignment_length"],
                "bitscore": hit["bitscore"],
            }
        )
    return typed, {item["genotype"] for item in vp1_metadata.values()}


def main():
    parser = argparse.ArgumentParser(
        description="Create final GI/GII genome references from original and supplemental FASTA files."
    )
    parser.add_argument("--ref-dir", type=Path, default=Path("ref_seq"))
    args = parser.parse_args()

    for tool in ("makeblastdb", "blastn"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"Required tool not found: {tool}")

    audit_rows = []
    with tempfile.TemporaryDirectory(prefix="merge_noro_genomes_") as temp:
        work_dir = Path(temp)

        for group, filenames in FILES.items():
            original_path = args.ref_dir / filenames["original"]
            genome_path = args.ref_dir / filenames["genome"]
            vp1_path = resolve_path(args.ref_dir, filenames["vp1"])
            output_path = args.ref_dir / filenames["output"]

            originals = list(SeqIO.parse(original_path, "fasta"))
            original_types = {
                original_genotype(record.id, group) for record in originals
            }
            original_accessions = {
                accession_from_identifier(record.id) for record in originals
            }
            original_index_by_accession = {
                accession_from_identifier(record.id): index
                for index, record in enumerate(originals)
            }

            candidates, target_types = run_blast(
                vp1_path,
                genome_path,
                group,
                work_dir,
            )
            missing_types = target_types - original_types
            best_by_type = {}
            for candidate in candidates:
                genotype = candidate["genotype"]
                rank = (
                    candidate["bitscore"],
                    candidate["alignment_length"],
                    candidate["identity"],
                )
                if genotype not in best_by_type or rank > best_by_type[genotype]["rank"]:
                    candidate["rank"] = rank
                    best_by_type[genotype] = candidate

            output_records = [
                SeqRecord(
                    clean_sequence(record.seq),
                    id=record.id,
                    description="",
                )
                for record in originals
            ]

            for genotype in sorted(missing_types):
                candidate = best_by_type.get(genotype)
                if candidate is None:
                    audit_rows.append(
                        {
                            "group": group,
                            "vp1_genotype": genotype,
                            "status": "not_available",
                            "accession": "",
                            "nearest_vp1_reference": "",
                            "vp1_identity_pct": "",
                            "vp1_alignment_length": "",
                            "sequence_length": "",
                        }
                    )
                    continue

                record = candidate["record"]
                accession = accession_from_identifier(record.id)
                if accession in original_accessions:
                    index = original_index_by_accession[accession]
                    previous_id = output_records[index].id
                    output_records[index] = SeqRecord(
                        clean_sequence(output_records[index].seq),
                        id=f"{genotype}_{accession}",
                        description=(
                            f"relabelled_from={previous_id} nearest_vp1="
                            f"{candidate['nearest_vp1_reference']}"
                        ),
                    )
                    status = "relabelled_existing"
                else:
                    output_records.append(
                        SeqRecord(
                            clean_sequence(record.seq),
                            id=f"{genotype}_{accession}",
                            description=(
                                f"supplemental_genome nearest_vp1="
                                f"{candidate['nearest_vp1_reference']}"
                            ),
                        )
                    )
                    original_accessions.add(accession)
                    original_index_by_accession[accession] = len(output_records) - 1
                    status = "added"

                audit_rows.append(
                    {
                        "group": group,
                        "vp1_genotype": genotype,
                        "status": status,
                        "accession": accession,
                        "nearest_vp1_reference": candidate[
                            "nearest_vp1_reference"
                        ],
                        "vp1_identity_pct": f"{candidate['identity']:.3f}",
                        "vp1_alignment_length": candidate["alignment_length"],
                        "sequence_length": len(clean_sequence(record.seq)),
                    }
                )

            SeqIO.write(output_records, output_path, "fasta")
            print(
                f"{group}: {len(originals)} original + "
                f"{len(output_records) - len(originals)} supplemental = "
                f"{len(output_records)} sequences"
            )

    audit_path = args.ref_dir / "genome_merge_audit.tsv"
    fields = [
        "group",
        "vp1_genotype",
        "status",
        "accession",
        "nearest_vp1_reference",
        "vp1_identity_pct",
        "vp1_alignment_length",
        "sequence_length",
    ]
    with audit_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)
    print(f"Audit table: {audit_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
