# Norovirus GI/GII Genotyping Pipeline

![version](https://img.shields.io/badge/version-2.2.0-blue)
![license](https://img.shields.io/badge/license-MIT-green)
![nextflow](https://img.shields.io/badge/nextflow-%3E%3D20.0.0-34a853?logo=nextflow&logoColor=white)

诺如病毒（Norovirus）GI/GII 核酸序列的单轮 BLAST 分型流程：根据 bitscore 最高的基因组参考判定 genogroup，用基因坐标表切割 RdRp/VP1 区段并建树确认型别，结合 SH-like 支持值判定置信度，并对区域外序列注释其所在的基因（NTPase、VPg、3CLpro、5'UTR 等）。

版本说明见 [`VERSION.md`](VERSION.md)。

## 目录

- [工作流程](#工作流程)
- [Docker 快速开始](#docker-快速开始)
- [本地安装](#本地安装)
- [参考库与坐标表准备](#参考库与坐标表准备)
- [运行](#运行)
- [输出](#输出)
- [HTML 报告](#html-报告)
- [目录结构](#目录结构)
- [可调参数](#可调参数)

## 工作流程

该流程对输入的 Norovirus 核酸序列进行单轮 BLAST 分型：

1. 将输入序列与 GI/GII 近完整基因组库进行一次 BLAST 比对，根据 bitscore 最高的命中参考判断为 GI 或 GII，因此位于非 RdRp/VP1 区域的片段也可以完成 genogroup 判定。
2. 用最佳基因组命中的 accession 在基因坐标表（`genome_gene_coordinates.tsv`）中查询 RdRp/VP1/ORF1/VP2 的基因组坐标，再结合 BLAST 比对坐标把 query 中对应的 RdRp/VP1 区段切出来，无需再做第二轮 region BLAST。
3. 同时把 query 的命中区间映射到 ORF1/ORF2/ORF3 和具体的蛋白（RdRp/VP1/VP2/NTPase/VPg/3CLpro/p48/5'UTR 等），注释 query 在基因组上的位置。
4. 对切出的 RdRp/VP1 片段分别用 MAFFT 加入对应区域参考序列比对，再用 FastTree（`-nome -mllen`，SH-like 支持值）建树。
5. 按最小树枝距离确定最接近的参考序列；当 query 与最近参考所在 clade 的 SH-like 支持值 ≥ `--bootstrap-threshold`（默认 0.75）时输出型别为 `typed`，否则输出型别但仍标记 `unsure`。
6. 对于完全落在 RdRp/VP1 之外的序列，输出 `status=off_region` 并列出最相似的几株参考基因组，供人工参考。

注意：建树不会使用完整基因组。流程用基因坐标表从 query 中切取 RdRp 或 VP1 区段，然后每次只取一条 query 片段，用 MAFFT 加入对应区域参考序列，再单独建树确认该条序列的最近参考。`--threads` 控制同时处理多少条 query。完整基因组只用于 GI/GII 粗分类和坐标推算。

## Docker 快速开始

最简单的运行方式是使用预构建的 Docker 镜像，镜像内已包含所有依赖（Biopython、BLAST+、MAFFT、FastTree、Nextflow）。

```bash
# 拉取镜像（或本地构建，见下文）
docker pull ghcr.io/jinglu-gdiph/noro-genotyping:2.2.0

# 运行（挂载输入、参考库、输出目录）
docker run --rm \
  -v "$(pwd)/input:/input" \
  -v "$(pwd)/ref_seq:/ref_seq:ro" \
  -v "$(pwd)/results:/results" \
  ghcr.io/jinglu-gdiph/noro-genotyping:2.2.0 \
  nextflow run /pipeline/main.nf \
    --sequences "/input/*.fasta" \
    --ref_dir /ref_seq \
    --outdir /results \
    --threads 8
```

也可以直接调用 Python 入口（跳过 Nextflow）：

```bash
docker run --rm \
  -v "$(pwd)/input:/input" \
  -v "$(pwd)/ref_seq:/ref_seq:ro" \
  -v "$(pwd)/results:/results" \
  ghcr.io/jinglu-gdiph/noro-genotyping:2.2.0 \
  python3 /pipeline/genotype_norovirus.py \
    --input /input/sample.fasta \
    --ref-dir /ref_seq \
    --output-prefix /results/sample \
    --threads 8 --report
```

### 本地构建镜像

```bash
docker build -t noro-genotyping:2.2.0 .
```

## 本地安装

需要 Python ≥ 3.9 和以下命令行工具：

| 工具 | 用途 | conda 安装 |
|---|---|---|
| BLAST+ (`makeblastdb`, `blastn`) | 序列比对 | `conda install -c bioconda blast` |
| MAFFT | 多序列比对 | `conda install -c bioconda mafft` |
| FastTree | 建树（含 SH-like 支持值） | `conda install -c bioconda fasttree` |
| Nextflow | 流程编排（可选） | `conda install -c bioconda nextflow` |
| Biopython + numpy | 序列/树解析 | `conda install -c conda-forge biopython numpy` |

一键安装所有依赖：

```bash
mamba env create -f environment.yml   # 或 conda env create -f environment.yml
conda activate norovirus-typing
```

## 参考库与坐标表准备

### 基因坐标表（一次性生成）

Norovirus 的 GenBank 记录通常只标注 ORF1（多聚蛋白）、ORF2（VP1）、ORF3（VP2）三个 CDS，**不单独标注 RdRp**（RdRp 位于 ORF1 的 3' 端）。因此坐标表的获取方式是：

- VP1/VP2/ORF1/UTR：直接从 GenBank 的 CDS/feature 注释解析。
- ORF1 细分蛋白（p48/NTPase/p22/VPg/3CLpro）：用文献已知的保守相对切割位点切分 ORF1。
- RdRp：用预裁剪好的 region 参考库（`gi_rdrp`/`gii_rdrp`）BLAST 回比到每条基因组，取最佳命中的 subject 坐标。
- 部分 accession 的 GenBank 拿不到或缺少 VP1 注释时，同样用 region 库回比补齐。

一次性生成坐标表：

```bash
python3 download_genome_annotations.py \
  --ref-dir ref_seq \
  --email your.email@example.com
```

可通过 `NCBI_EMAIL` / `NCBI_API_KEY` 环境变量提供联系邮箱和 API key。脚本会把 GenBank 缓存到 `ref_seq/genbank_cache/`（**不应提交到版本库**），再次运行只增量下载缺失项；坐标表写到 `ref_seq/genome_gene_coordinates.tsv`，列：

```text
group  accession  genotype  sequence_length  protein  start  end  strand  source
```

`source` ∈ `genbank` / `inferred_from_region_blast` / `orf1_relative`。

### 参考序列

默认读取：

```text
ref_seq/
  gi_with_genotype.fasta      # GI 近完整基因组（含型别标签）
  gi_rdrp.fasta               # GI RdRp 参考片段
  gi_VP1.fa                   # GI VP1 参考片段
  gii_with_genotype.fasta     # GII 近完整基因组
  gii_rdrp.fasta              # GII RdRp 参考片段
  gii_vp1.fasta               # GII VP1 参考片段
  genome_gene_coordinates.tsv # download_genome_annotations.py 生成
```

区域参考序列标题必须包含基因型，当前支持：

```text
>U07611_GII.P1
>AJ277606_GII.1
>I.1|M87661
```

基因组库由带型别的 original 文件作为主体，并从 `gi_genome.fasta`、`gii_genome.fasta` 中补充缺失 VP1 型别：

```bash
python3 merge_genome_references.py --ref-dir ref_seq
```

该命令生成 `gi_with_genotype.fasta`、`gii_with_genotype.fasta` 和 `genome_merge_audit.tsv`。original 文件不会被修改。

## 运行

### Nextflow

```bash
nextflow run main.nf \
  --sequences "input/*.fasta" \
  --ref_dir ref_seq \
  --outdir results \
  --threads 8 \
  --aligner mafft
```

默认使用 MAFFT 的片段加入模式，并启用 `--quiet` 减少 stderr 进度输出。每条 query 序列会分别与 RdRp/VP1 参考序列比对并建树，避免多条 query 互相影响树拓扑。

其他流程调用时推荐固定以下入口：

```bash
nextflow run /path/to/noro_genotyping/main.nf \
  --sequences "/path/to/query/*.fasta" \
  --ref_dir /path/to/noro_genotyping/ref_seq \
  --outdir /path/to/output \
  --threads 8
```

### Python 直接调用

```bash
python3 genotype_norovirus.py \
  --input query.fasta \
  --ref-dir ref_seq \
  --output-prefix output/sample \
  --threads 8 \
  --report
```

## 输出

`results/norovirus_genotyping.tsv` 为所有输入文件的合并信息表，主要字段包括：

- `genogroup`：相似性粗分类得到的 GI/GII。
- `coarse_reference`：基因组粗分类阶段的最佳参考序列。
- `query_genome_region`、`query_genome_proteins`、`query_genome_coords`：query 在最佳基因组参考上的 ORF 区/蛋白/坐标注释。
- `rdrp_nearest_reference`、`rdrp_genotype`、`rdrp_genotype_confidence`、`rdrp_bootstrap`：RdRp 树上最近参考、P-type、置信度（typed/unsure）和 SH-like 支持值。
- `vp1_nearest_reference`、`vp1_genotype`、`vp1_genotype_confidence`、`vp1_bootstrap`：VP1 树上最近参考、capsid genotype、置信度和 SH-like 支持值。
- `*_tree_distance`：查询序列与最近参考之间的树枝距离。
- `*_blast_identity_pct`、`*_reference_coverage_pct`：切片阶段的相似度和参考区域覆盖度。
- `off_region_top_references`：落在 RdRp/VP1 之外的 query 列出的最相似参考基因组。
- `status`：`typed`、`partial_genotype`、`genotype_unsure`、`off_region`、`genogroup_only` 或 `unclassified`。

其他输出：

```text
results/norovirus_genotyped.fasta
results/samples/<sample>_genotyping.tsv
results/samples/<sample>_genotyped.fasta
results/samples/<sample>_trees/
results/samples/<sample>_report.html
```

## HTML 报告

`--report`（Nextflow 默认开启，`--generate_report false` 关闭）会在每个样本旁生成一个
**自包含 HTML 表格报告** `<sample>_report.html`，包括：

- 顶部统计卡片（typed / partial_genotype / genotype_unsure / off_region / unclassified 数量）。
- RdRp（P-type）基因型表：每条 query 的型别、置信度（typed/unsure）、SH 支持值、最近参考、树枝距离、BLAST 一致性、参考覆盖度。
- VP1（capsid）基因型表：同上。
- 基因组位置表：每条 query 的 ORF 区（ORF1/ORF2/ORF3/multi）、具体蛋白（RdRp/VP1/NTPase/VPg/3CLpro/p48/5'UTR 等，完整基因组覆盖多蛋白时只显示前 5 个并标注 "more"）、参考坐标、粗分类参考、以及 off-region 序列的最相似参考列表。

也可单独对已有 TSV 重新生成报告（无需 trees 目录）：

```bash
python3 generate_report.py \
  --tsv results/samples/<sample>_genotyping.tsv \
  --output results/samples/<sample>_report.html
```

## 目录结构

```text
noro_genotyping/
├── main.nf                      # Nextflow 工作流入口
├── nextflow.config              # Nextflow 配置（资源、profile）
├── environment.yml              # conda 环境定义
├── Dockerfile                   # 容器镜像构建
├── .gitignore
├── genotype_norovirus.py        # 核心分型脚本（单轮 BLAST + 坐标切割 + 建树）
├── download_genome_annotations.py  # 一次性下载 GenBank 并生成坐标表
├── merge_genome_references.py   # 基因组参考库合并（补充缺失 VP1 型别）
├── generate_report.py           # HTML 表格报告生成
├── VERSION.md                   # 版本与变更日志
└── ref_seq/                     # 参考序列与坐标表
    ├── gi_with_genotype.fasta
    ├── gii_with_genotype.fasta
    ├── gi_rdrp.fasta
    ├── gii_rdrp.fasta
    ├── gi_VP1.fa
    ├── gii_vp1.fasta
    ├── genome_gene_coordinates.tsv
    └── genbank_cache/           # （git-ignored）GenBank 缓存
```

## 可调参数

```text
--min_identity 70            # 粗分类最低一致性
--min_query_coverage 50      # 粗分类最低 query 覆盖度
--min_slice_identity 70      # 进入 RdRp/VP1 建树的最低基因组命中一致性
--min_slice_length 300       # 进入建树的最小切片长度（nt）
--bootstrap_threshold 0.75   # SH 支持值阈值，低于则标 unsure
--off_region_top_n 3         # off-region 序列列出的最相似参考数量
```

GI/GII 粗分类覆盖度以输入序列长度为分母；RdRp/VP1 切片覆盖度以对应参考区域长度为分母。低于 `--min_slice_identity` 或切片短于 `--min_slice_length` 的序列不会进入建树。bootstrap 支持值低于 `--bootstrap_threshold` 的型别仍会输出但标记 `unsure`。
