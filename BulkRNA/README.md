# Bulk RNA-seq

Gene-level differential expression, over-representation analysis (ORA), gene set enrichment analysis (GSEA) and expression heatmaps.

## Requirements

R with data.table, DESeq2, clusterProfiler, fgsea, BiocParallel, ggplot2, ggtext, patchwork, pheatmap, gtable and ragg. The statistical implementation uses R 4.2.1, DESeq2 1.38.3, clusterProfiler 4.6.2 and fgsea 1.24.0. Plots use Arial; select an installed font if needed.

## Run

Set the working directory to BulkRNA:

```r
source('00_prepare_shared_data.R', encoding = 'UTF-8')
source('Extended_Data_Figure_8.R', encoding = 'UTF-8')
source('Extended_Data_Figure_10.R', encoding = 'UTF-8')
```

Preparation writes shared tables to prepared/tables/. Plotting scripts read those tables and write PNG/PDF figures and plotting data to results/<figure>/plots/ and tables/. To change plot appearance only, rerun the relevant figure script.

| Script | Contents |
|---|---|
| Extended_Data_Figure_8.R | Pregnancy volcano, public/EPO GSEA and Pregnancy/EPO heatmap |
| Extended_Data_Figure_10.R | EPO/EMP9 volcanoes, ORA, GSEA and treatment heatmap |
| Figure_3b.R | GO barplots from the two specified CSV inputs |
| 00_process_reads.sh | Optional FASTQ processing |
| 01_import_featurecounts.R | Assemble reprocessed count, length and TPM tables |

## Included inputs

- data/counts/: raw gene counts; one gene per row and sample per column.
- data/samples.tsv: unique source/sample identifiers, condition and tissue.
- data/comparisons.tsv and comparison_samples.tsv: comparison definitions, model and sample membership.
- data/exclusions.tsv: recorded exclusions and balancing selections. A balancing selection is not a QC failure.
- data/annotation/: gene lengths, frozen GO/Hallmark memberships, gene exclusions and plot labels.
- data/tpm/PregEMP90807.tsv.gz: the supplied full 0807 TPM matrix; Paper/0804 TPM is calculated from counts and gene lengths.

The 35 source-labelled samples span Paper, 0804, 0807 and public Pregnancy datasets. Reusing a sample across comparisons does not create another biological replicate. Models are unpaired; no animal pairing is inferred from sample names.

| Comparison | Samples | Design |
|---|---|---|
| Pregnancy SP/BM | 6/6, Paper + 0807 | ~ cohort + group |
| EPO SP/BM | 3/3 | ~ group |
| Pregnancy + EMP9 SP/BM | 3/3 | ~ group |
| EPO/saline SP | 3/3 | ~ group |
| Pregnancy + EMP9/Pregnancy SP | 3/3 | ~ group |
| Public Pregnancy SP/BM | 4/4 | ~ group |

## Statistical definitions

DESeq2 uses raw counts and unshrunken log2 fold changes. DEG selection requires raw P < 0.01, BH-adjusted P < 0.20 and abs(log2FC) >= log2(1.5), not abs(log2FC) >= 1.5. Result extraction uses alpha = 0.25 for independent filtering before the final 0.20 selection. All-zero rows are removed for Pregnancy and EPO SP/BM as specified in comparisons.tsv.

Pregnancy SP-high Top200 genes are positive DEGs sorted by decreasing log2FC, then adjusted P, raw P and gene. The EPO-responsive signature contains all EPO-up SP DEGs. Heatmap-specific exclusions do not alter the GSEA signatures.

ORA uses the frozen TERM2GENE table, a comparison-specific background of finite-P genes mapped to that annotation, set sizes 5–500 and BH correction across the tested terms. Display terms are selected after testing; full results, query and background are exported.

GSEA uses all finite log2FC values with weighted exponent 1, standard score, sizes 10–1000, eps = 0, serial execution and seed 20261003. Deterministic gene ordering, the 1e-12 tie offset and reversed rank/score direction are explicit in R/prepare_GSEA.R. One predefined signature is tested per comparison. NES/FDR labels come from the calculated table.

## Figure values

Volcanoes show log2FC versus −log10(raw P). Display coordinates are capped at x ±10 and y 20/14/40 for Pregnancy/EPO/EMP9; classification uses uncapped statistics. Missing values are not replaced by zero.

Heatmaps show gene-wise Z-scores of log2(TPM + 1), standardized separately within Paper, 0807 and 0804 sample sets. Display colors are capped at ±2; undefined Z-scores remain missing. The Pregnancy/EPO matrix contains 35 genes and 18 samples. The treatment matrix contains 24 genes and 12 samples: 22 strict overlap genes plus Klf1/Epor as display genes, not additional strict dual-cutoff DEGs.

Figure_3b.R separately reads data/legacy_Figure3b/preg sp vc preg bm.csv and preg sp vs non-preg sp.csv. These CSVs are not bundled. Its read sequence and same-name output reuse follow that script; it does not use the count-based comparisons above.

## Optional FASTQ processing

Place the required files at the relative locations in data/fastq_files.tsv and set the matching mouse reference paths:

```bash
STAR_INDEX=/path/to/mm10_index GTF=/path/to/genes.gtf bash 00_process_reads.sh
Rscript 01_import_featurecounts.R
```

The recipe uses fastp, STAR and paired-fragment featureCounts (-p --countReadPairs -g gene_name). It is separate from reproduction using the included counts and does not overwrite them. FASTQs and reference files are not bundled. Some listed files already contain merged sequencing lanes; do not add the corresponding unmerged files again. The complete Paper/Public raw-read processing history is not represented by this August processing recipe.
