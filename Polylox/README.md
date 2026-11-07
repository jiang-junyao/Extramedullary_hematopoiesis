# Polylox clonal analysis

Barcode tissue distribution, tissue sharing, sampling controls, lineage-output and split-sample analyses.

## Requirements and run

R 4.4.1 with dplyr, tidyr, ggplot2 and scales. PDF export requires Cairo support.

Set the working directory to Polylox and prepare the shared data:

```r
source('00_prepare_shared_data.R', encoding = 'UTF-8')
```

Run the plotting script required:

```r
source('Figure_1.R', encoding = 'UTF-8')
source('Figure_2.R', encoding = 'UTF-8')
source('Extended_Data_Figure_1.R', encoding = 'UTF-8')
source('Extended_Data_Figure_2.R', encoding = 'UTF-8')
source('Supplementary_Figure_3.R', encoding = 'UTF-8')
```

Preparation produces tissue, simulation and lineage-output tables under prepared/. Scripts write figures and underlying values to their own results/ subdirectories. Inputs are reused rather than copied for each plot.

Additional plotting scripts are provided in optional/. Run them from the Polylox working directory after shared preparation; do not change into optional/. Lineage-output scripts use prepared/output/tables/output_data.rds.

For split-sample and barcode-richness analyses:

```r
source('optional/In_silico_split_sample.R', encoding = 'UTF-8')
source('optional/Barcode_richness_sensitivity.R', encoding = 'UTF-8')
```

## Inputs

| Location | Contents |
|---|---|
| data/tissue/raw/ | Barcode count/Pgen tables; Code is a character identifier |
| data/tissue/metadata/ | Library metadata and progenitor-population strategy |
| data/tissue/formal/ | Supplied source-data percentages, with workbook/sheet/cell provenance |
| data/simulation/raw/ | Count/Pgen inputs for sampling and mixing analyses |
| data/simulation/library_metadata.tsv | Library-to-record and tissue/population mapping |
| data/simulation/shuffle_data_all.csv | Supplied uniform-mixing iterations |
| data/simulation/reference/ | Reference matrices used to check count preparation |
| data/output/raw/ | Lineage-output count/Pgen tables |
| data/output/libraries.tsv and mice.csv | Technical libraries and animal conditions |

Required columns include Code, pgen and sample/library counts. Metadata identifies mouse, sample, tissue, cell type and replicate; sample labels must match the count columns. Record-specific population strategies distinguish HSC/MPP, LSK and LSM inputs.

## Definitions

- A clone key is mouse/record + Code. Identical barcode strings from different animals are not merged into one clone.
- Read counts are nonnegative integers. Single-read observations are removed. Technical repeats are not biological replicates.
- Normalization and merging are specified separately in R/prepare_tissue.R and R/prepare_output.R; their profile-combination rules are not interchangeable.
- Pgen is used as supplied, not recalculated. Clonal selection uses Pgen<1e-5; zero is valid and missing values are not imputed.
- Tissue-specific means observed in one measured tissue; shared means observed in multiple measured tissues within the same animal. Co-detection is not direct evidence of migration.
- Tissue percentage denominators and selected populations are exported with the prepared tables. Heatmap row percentages describe the distribution of an individual barcode across the selected tissues.
- Lineage-output preparation distinguishes spleen-local output from combined Leg/Spine/Spleen output. Its pooled barcode-level Wilcoxon tests are exploratory, not tests with independent mice as the observations.

## Plot scope

Some percentage plots use the supplied source-data values and do not recompute inferential statistics from barcode counts. Heatmap reconstruction scripts use the included records identified in each script; exact manuscript panel-to-record correspondence is unverified for some examples. The scripts retain that distinction.

Sampling scripts separate observed values, supplied mixing iterations and regenerated simulations. Seeds, downsampling fractions, iterations and error/display summaries are specified in the relevant R/prepare_simulation*.R and plotting scripts.
