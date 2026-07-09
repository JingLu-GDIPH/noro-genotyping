#
# Norovirus GI/GII genotyping pipeline — self-contained Docker image.
#
# The image bundles every runtime dependency (Biopython, BLAST+, MAFFT,
# FastTree, Nextflow) in a conda environment so the pipeline runs anywhere
# Docker is available, with no host-side installs.
#
# Build:   docker build -t noro-genotyping:2.2.0 .
# Run:     docker run --rm -v "$PWD/input:/input" -v "$PWD/ref_seq:/ref_seq:ro" \
#                   -v "$PWD/results:/results" noro-genotyping:2.2.0 \
#                   nextflow run /pipeline/main.nf --sequences "/input/*.fasta" \
#                   --ref_dir /ref_seq --outdir /results --threads 8

FROM mambaorg/micromamba:1.5.10-jammy

LABEL org.opencontainers.image.title="noro-genotyping"
LABEL org.opencontainers.image.description="Norovirus GI/GII single-pass BLAST genotyping with RdRp/VP1 phylogenetic typing"
LABEL org.opencontainers.image.version="2.2.0"
LABEL org.opencontainers.image.source="https://github.com/JingLu-GDIPH/noro-genotyping"
LABEL org.opencontainers.image.licenses="MIT"

# --- 1. Install the conda environment from environment.yml -------------------
# Copy only the environment file first so dependency installation is cached
# across code-only changes.
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# mambaorg/micromamba uses $MAMBA_USER (non-root) and roots the env in
# /opt/conda. The container's default env name comes from environment.yml
# ("norovirus-typing"); we install it into the base env for simplicity.
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all -y \
    && rm -rf /tmp/environment.yml

# --- 2. Copy the pipeline code and reference data ----------------------------
WORKDIR /pipeline
COPY --chown=$MAMBA_USER:$MAMBA_USER \
    main.nf nextflow.config \
    genotype_norovirus.py \
    download_genome_annotations.py \
    merge_genome_references.py \
    generate_report.py \
    VERSION.md README.md \
    /pipeline/

# Reference FASTAs + the precomputed gene coordinate table. The GenBank cache
# is intentionally NOT baked in (it is regenerable and large); users who want
# to refresh coordinates mount their own ref_seq or run the download script.
COPY --chown=$MAMBA_USER:$MAMBA_USER ref_seq/ /pipeline/ref_seq/

# --- 3. Runtime entrypoint ---------------------------------------------------
# The conda env is installed into base (/opt/conda), so putting /opt/conda/bin
# on PATH is enough for python3/blastn/mafft/fasttree/nextflow to resolve. We do
# NOT wrap the entrypoint in `micromamba run` because that helper writes to a
# cache dir that is not writable when Nextflow mounts work dirs as root.
ENV PATH=/opt/conda/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV CONDA_PREFIX=/opt/conda

# Working directory for the user's data; the pipeline writes outputs wherever
# --outdir/--output-prefix point (typically a mounted volume).
WORKDIR /data

CMD ["nextflow", "run", "/pipeline/main.nf", "--help"]
