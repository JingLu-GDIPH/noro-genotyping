# Norovirus GI/GII Genotyping Pipeline

![version](https://img.shields.io/badge/version-2.2.0-blue)
![license](https://img.shields.io/badge/license-MIT-green)

诺如病毒（Norovirus）GI/GII 核酸序列的单轮 BLAST 分型流程：根据 bitscore 最高的基因组参考判定 genogroup，用基因坐标表切割 RdRp/VP1 区段并建树确认型别，结合 SH-like 支持值判定置信度，并对区域外序列注释其所在的基因（NTPase、VPg、3CLpro、5'UTR 等）。

## 工作流程

1. 将输入序列与 GI/GII 近完整基因组库进行一次 BLAST 比对，根据 bitscore 最高的命中参考判断为 GI 或 GII。
2. 用最佳基因组命中的 accession 在基因坐标表（`genome_gene_coordinates.tsv`）中查询 RdRp/VP1 等坐标，结合 BLAST 比对坐标把 query 中对应的 RdRp/VP1 区段切出来，无需再做第二轮 region BLAST。
3. 同时把 query 的命中区间映射到具体的蛋白（RdRp/VP1/NTPase/VPg/3CLpro/p48/5'UTR 等），注释 query 在基因组上的位置。
4. 对切出的 RdRp/VP1 片段分别用 MAFFT 加入对应区域参考序列比对，再用 FastTree（SH-like 支持值）建树。
5. 按最小树枝距离确定最接近的参考序列；当 SH-like 支持值 ≥ `--bootstrap-threshold`（默认 0.75）时输出型别为 `typed`，否则输出型别但标记 `unsure`。
6. 对于完全落在 RdRp/VP1 之外的序列，输出 `status=off_region` 并列出最相似的几株参考基因组。

注意：建树不会使用完整基因组。每次只取一条 query 的 RdRp 或 VP1 切片，用 MAFFT 加入对应区域参考序列，再单独建树确认该条序列的最近参考。

## 快速开始

只需 Nextflow 和 Docker，无需安装任何其他依赖。Nextflow 会自动拉取已包含全部依赖（Biopython、BLAST+、MAFFT、FastTree）的 Docker 镜像。

```bash
nextflow run JingLu-GDIPH/noro-genotyping \
  -profile docker \
  --sequences "input/*.fasta" \
  --outdir results \
  --threads 8
```

首次运行会自动拉取镜像（约 2 GB）。参考库（GI/GII 基因组、RdRp/VP1 片段、基因坐标表）已随仓库提供，无需额外下载。

<details>
<summary>不使用 Docker（本地 conda 运行）</summary>

```bash
# 安装依赖
mamba env create -f environment.yml && conda activate norovirus-typing

# 克隆并运行
git clone https://github.com/JingLu-GDIPH/noro-genotyping.git
cd noro-genotyping
nextflow run main.nf \
  --sequences "input/*.fasta" \
  --ref_dir ref_seq \
  --outdir results \
  --threads 8
```

也可直接调用 Python 入口：

```bash
python3 genotype_norovirus.py \
  --input query.fasta \
  --ref-dir ref_seq \
  --output-prefix output/sample \
  --threads 8 \
  --report
```

</details>

## 输出

`results/norovirus_genotyping.tsv` 为合并信息表，主要字段：

- `genogroup`：GI/GII 粗分类。
- `coarse_reference`：粗分类最佳参考序列。
- `query_genome_region`、`query_genome_proteins`、`query_genome_coords`：query 在基因组上的 ORF 区/蛋白/坐标注释。
- `rdrp_genotype`、`rdrp_genotype_confidence`、`rdrp_bootstrap`：RdRp P-type、置信度（typed/unsure）、SH 支持值。
- `vp1_genotype`、`vp1_genotype_confidence`、`vp1_bootstrap`：VP1 capsid 型别、置信度、SH 支持值。
- `off_region_top_references`：落在 RdRp/VP1 之外的 query 列出的最相似参考。
- `status`：`typed`、`partial_genotype`、`genotype_unsure`、`off_region`、`genogroup_only` 或 `unclassified`。

`--report`（Nextflow 默认开启）额外生成自包含的 **HTML 表格报告** `<sample>_report.html`，汇总状态统计、RdRp/VP1 基因型表和基因组位置表。

其他输出：

```text
results/norovirus_genotyped.fasta
results/samples/<sample>_genotyping.tsv
results/samples/<sample>_genotyped.fasta
results/samples/<sample>_trees/
results/samples/<sample>_report.html
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

## 参考库与坐标表更新

参考库已随仓库提供。如需重新生成基因坐标表（例如补充了新的基因组参考）：

```bash
python3 download_genome_annotations.py \
  --ref-dir ref_seq \
  --email your.email@example.com
```

该脚本从 GenBank 解析 VP1/VP2/ORF1/UTR 坐标，并用 RdRp 区域库 BLAST 回比推算 RdRp 坐标，同时把 ORF1 细分为 p48/NTPase/p22/VPg/3CLpro。可通过 `NCBI_EMAIL` / `NCBI_API_KEY` 环境变量提供邮箱和 API key。

## 许可证

MIT，详见 [LICENSE](LICENSE)。版本说明见 [VERSION.md](VERSION.md)。
