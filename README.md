# Evaluate Enhancer-Gene prediction models against experimental data

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
 produce PR curves. If only positives are available, the plotting module will fail but useful
 intermediate files will still be produced
 * Each prediction file may have multiple score columns. The names of these columns may be the same
 across prediction files - they will be differentiated based on dataset name (as defined by the
 prediction table). The combination of dataset name and prediction column should be unique. TO DO:
 Check for this.
 * The code currently overlaps on GeneSymbol (not on gene TSS coordinates)

## Requirements
 * Required Inputs (see below for formats)
 	* One experimental data file 
 	* At least one predictions file
 	* Configuration file describing how predictor should be aggregated
 	
### Dependencies
Dependencies are automatically handled via conda when executing the workflow via snakemake (see
below). Conda environment files can be found in `workflow/envs`. If comparisons are performed
without conda, following R dependencies are required:

```
# base R:
R (>=4.0.0)

# R packages:
tidyverse (>=1.3.0)
here (>=1.0.1)
data.table (>=1.13.6)
GenomicRanges (>=1.42.0)
ROCR (>=1.0_11)
caTools (>=1.18.1)
rmarkdown (>=2.7)
bookdown (>=0.21)
DT (>=0.18)
plotly (>=4.9.3)
```

## File Formats

* Experimental Data
  * <https://docs.google.com/spreadsheets/d/1Tl_fdPeAeiVkZnettWxeMJW3zL5zfFRcC2TV-vmQRLQ/edit#gid=257851280>. 
  * Not all of these columns are required. See `resources/example/input/K562.ExperimentalData.slim.txt`
  for list of required columns

* Predictions
  * <https://docs.google.com/spreadsheets/d/1BQBFC4PzPA8v3tA_OkSp2lU1YpO74uwmkEa2TJTt-ic/edit#gid=0>
  * See `resources/example/input/K562.ABC.Predictions.AvgHiC.chrX.ENCODE.format.txt.gz` for example. 

## Configuring the snakemake workflow
The config/config.yml file is used to specify comparisons that should be performed. See the
comparison "example" in this file as an example. In addition to the predictions and experiment input
files, each comparison requires a prediction config file (pred_config) in .txt format. This file
specifies how the predictions should be handled and the behavior of the comparison code for this
predictor depends on it's content:

 * pred.col: must match the name of the column in the predictions file
 * agg.func: In the case that an experimentally tested element overlaps multiple predicted elements,
 how should the predicted elements be aggregated to the level of the tested element. 
 * fill.val: In the case that an experimentally tested element does not overlap any predicted
 elements, what value should be filled in.
 * lowerIsMoreConfident: Set to TRUE if lower values of the predictor signify more confidence in the
 prediction. This is appropraite for predictors such as linear distance or pvalue. It is generally
 preferred for this column to be FALSE.
 
See `resources/example/input/pred_confg.txt` for an example.

## Sample command
Following command can be used to perform all specified comparisons using snakemake with conda to
automatically handle dependencies. This requires that snakemake (>=5.10.0) and conda
(e.g. miniconda) are installed.

```
# perform all comparisons specified in config.yml (-n = dryrun, remove for execution)
snakemake --use-conda -n
```
