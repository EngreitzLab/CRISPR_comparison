# Benchmark enhancer-gene prediction models against CRISPR data

This workflow is designed to evaluate the performance of an enhancer-gene linking model against
experimental data from CRISPR enhancer screes. It supports the evaluation of multiple predictors
against a single experimental data file. Comparing performance of a single predictor against
multiple experiments currently not supported.

Fundamentally, the idea is to overlap each experimentally tested element-gene-celltype tuple with a
predicted element-gene-celltype tuple. Diagnostic plots such as PR curves are then produced on this
overlapped dataset. Care must be taken when an experimentally tested element overlaps multiple
predicted elements, or when an experimentally tested element does not overlap any predicted
elements. See the configuration section below for how to handle these cases. 

Other notes:

 * There must be both experimental positives and negatives in the experimental data file in order to
 produce PR curves.
 * The code currently overlaps based on gene symbols (not on gene TSS coordinates)

## Requirements
 * Inputs (see below for formats):
 	* One experimental data file 
 	* At least one predictions file
 	* (optional) Configuration file describing metadata of predictors, including how to aggreagate
  multiple predicted enhancers overlapping one experimental enhancer
 	* (optional) Cell type mappings between cell types in predictions and experiment
 	
### Dependencies
Dependencies are automatically handled via conda when executing the workflow via snakemake (see
below). Conda environment files can be found in `workflow/envs`. If comparisons are performed
without conda, following R dependencies are required:

```sh
# base R:
R (>=4.1.1)

# R packages:
R.utils (>=2.11.0)
data.table (>=1.14.2)
tidyverse (>=1.3.1)
cowplot (>=1.1.1)
ggpubr (>=0.4.0)
ggcorrplot (>=0.1.3)
upsetr (>=1.4.0)
plotly (>=4.10.0)
rocr (>=1.0_11)
catools (>=1.18.2)
boot (>=1.3_28.1)
dt (>=0.22)
rmarkdown (>=2.13)
bookdown (>=0.25)
optparse (>=1.7.1)
GenomicRanges (>=1.46.1)
rtracklayer (>=1.54.0)
BiocParallel (>=1.28.3)
```

## File Formats

* Experimental Data
  * <https://docs.google.com/spreadsheets/d/1xold84upBFigZPFQGUMeNC5t2id-Bub_Y1o5Xid_MHw/edit?usp=sharing>. 
  * See `resources/example/experimental_data_chrX.txt.gz` as an example.

* Predictions
  * <https://docs.google.com/spreadsheets/d/1xold84upBFigZPFQGUMeNC5t2id-Bub_Y1o5Xid_MHw/edit?usp=sharing>
  * See `resources/example/K562_ABC_K562HiC_chrx.txt.gz` as an example. 

## Configuring the snakemake workflow
The `config/config.yml` file is used to specify comparisons that should be performed. See the
comparison `"example"` in this file as an example. In addition to the predictions and experiment
input files, each comparison can take a prediction config file (pred_config) in .txt format as input.
This file specifies how the different predictors should be handled and the behavior of the comparison
code for this predictor depends on its content:

 * pred_id: Short name for each predictor. Same as the names of 'pred' in the `config.yml` file.
 * pred_col: Column name in prediction file containing the predictor values. Following the
 aforementioned file format this is typlically 'Score'.
 * boolean: (TRUE/FALSE) is this predictor binary, or os this a quantitative prediction score?
 * alpha: Predictor score cutoff for plots. Can be NA if not applicable or unknown.
 * aggregate_function: In the case that an experimentally tested element overlaps multiple predicted
 elements, how should the predicted elements be aggregated to the level of the tested element. 
 * fill_value: In the case that an experimentally tested element does not overlap any predicted
 elements, what value should be filled in.
 * inverse_predictor: Set to TRUE if lower values of the predictor signify more confidence in the
 prediction. This is appropriate for predictors such as linear distance or pvalue.
 * pred_name_long: A pretty name (2-3 words) for the predictor to make plots and tables look nicer.
 
See `resources/example/pred_config.txt` for an example. If this file is left out (`NULL` in 
`config.txt`), a file with default values will be generated, however they might not be appropriate
for the provided predictors.

## Sample command
Following command can be used to perform all specified comparisons using snakemake with conda to
automatically handle dependencies. This requires that snakemake (>=5.10.0) and conda
(e.g. [miniconda](https://docs.conda.io/en/latest/miniconda.html)) are installed.

```sh
# perform all comparisons specified in config.yml (-n = dryrun, remove for execution)
snakemake --use-conda -j1 -n
```

All generated output including the main .html document are saved to `results/example/`. For other 
comparison, the subdirectory in results will be named after the comparison name as specified in the
`config.yml` file.
