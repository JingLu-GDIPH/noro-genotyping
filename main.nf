#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.sequences = null
params.ref_dir = "${projectDir}/ref_seq"
params.outdir = "./results"
params.threads = 4
params.aligner = "mafft"
params.min_identity = 70.0
params.min_query_coverage = 50.0
params.min_slice_identity = 70.0
params.min_slice_length = 300
params.bootstrap_threshold = 0.75
params.off_region_top_n = 3
params.generate_report = true

if (!params.sequences) {
    error "Please provide input FASTA files with --sequences"
}

process NOROVIRUS_GENOTYPING {
    tag "${sample_name}"
    publishDir "${params.outdir}/samples", mode: 'copy'

    conda "${projectDir}/environment.yml"

    input:
    tuple val(sample_name), path(sequence_file)
    tuple path(gi_genome), path(gi_rdrp), path(gi_vp1), path(gii_genome), path(gii_rdrp), path(gii_vp1), path(gene_coordinates)

    output:
    path "${sample_name}_genotyping.tsv", emit: tables
    path "${sample_name}_genotyped.fasta", emit: fastas
    path "${sample_name}_trees", emit: trees
    path "${sample_name}_report.html", emit: reports, optional: true

    script:
    """
    mkdir ref_seq
    cp ${gi_genome} ref_seq/gi_genome.fasta
    cp ${gi_rdrp} ref_seq/gi_rdrp.fasta
    cp ${gi_vp1} ref_seq/gi_vp1.fasta
    cp ${gii_genome} ref_seq/gii_genome.fasta
    cp ${gii_rdrp} ref_seq/gii_rdrp.fasta
    cp ${gii_vp1} ref_seq/gii_vp1.fasta
    cp ${gene_coordinates} ref_seq/genome_gene_coordinates.tsv

    python3 ${projectDir}/genotype_norovirus.py \
        --input ${sequence_file} \
        --ref-dir ref_seq \
        --output-prefix ${sample_name} \
        --threads ${params.threads} \
        --aligner ${params.aligner} \
        --min-identity ${params.min_identity} \
        --min-query-coverage ${params.min_query_coverage} \
        --min-slice-identity ${params.min_slice_identity} \
        --min-slice-length ${params.min_slice_length} \
        --bootstrap-threshold ${params.bootstrap_threshold} \
        --off-region-top-n ${params.off_region_top_n} \
        ${params.generate_report ? '--report' : ''}
    """
}

process MERGE_RESULTS {
    tag "merge genotyping tables"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tables
    path fastas

    output:
    path "norovirus_genotyping.tsv"
    path "norovirus_genotyped.fasta"

    script:
    """
    python3 - <<'PY'
import csv
from pathlib import Path

tables = sorted(Path(".").glob("*_genotyping.tsv"))
with open("norovirus_genotyping.tsv", "w", newline="") as output:
    writer = None
    for table in tables:
        with table.open() as handle:
            reader = csv.DictReader(handle, delimiter="\\t")
            if writer is None:
                writer = csv.DictWriter(output, fieldnames=reader.fieldnames, delimiter="\\t")
                writer.writeheader()
            writer.writerows(reader)

with open("norovirus_genotyped.fasta", "w") as output:
    for fasta in sorted(Path(".").glob("*_genotyped.fasta")):
        output.write(fasta.read_text())
PY
    """
}

workflow {
    sequence_channel = Channel
        .fromPath(params.sequences, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }

    reference_channel = Channel.value(
        tuple(
            file("${params.ref_dir}/gi_with_genotype.fasta", checkIfExists: true),
            file("${params.ref_dir}/gi_rdrp.fasta", checkIfExists: true),
            file("${params.ref_dir}/gi_VP1.fa", checkIfExists: true),
            file("${params.ref_dir}/gii_with_genotype.fasta", checkIfExists: true),
            file("${params.ref_dir}/gii_rdrp.fasta", checkIfExists: true),
            file("${params.ref_dir}/gii_vp1.fasta", checkIfExists: true),
            file("${params.ref_dir}/genome_gene_coordinates.tsv", checkIfExists: true)
        )
    )

    NOROVIRUS_GENOTYPING(sequence_channel, reference_channel)
    MERGE_RESULTS(
        NOROVIRUS_GENOTYPING.out.tables.collect(),
        NOROVIRUS_GENOTYPING.out.fastas.collect()
    )
}

workflow.onComplete {
    log.info "Norovirus dual-region genotyping completed."
    log.info "Results: ${params.outdir}/norovirus_genotyping.tsv"
}
