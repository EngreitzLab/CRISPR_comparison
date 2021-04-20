## Merge predictions with experimental data for downstream comparisons of predictions with CRISPR
## CRE perturbation data

# attach required packages and functions
suppressPackageStartupMessages(library(data.table))
source(file.path(snakemake@scriptdir, "crisprComparisonFunctions.R"))

# load experimental data
message("Reading experimental data in: ", snakemake@input$experiment)
expt <- loadFileString(snakemake@input$experiment)
message("\tLoaded experimental data with ", nrow(expt), " rows.")

# table containing prediction files and parameters
pred_table <- data.frame(
  name = snakemake@params$pred_names,
  path = snakemake@input$predictions,
  threshold = unlist(snakemake@params$pred_threshold[snakemake@params$pred_names]),
   stringsAsFactors = FALSE)

# load all predictions
pred_list <- loadPredictions(pred_table)

# load pred_config file
pred_config <- fread(snakemake@input$pred_config)

# load cell mapping file if provided
if (!is.null(snakemake@params$cell_name_mapping)) {
  cellMapping <- fread(snakemake@params$cell_name_mapping)
} else {
  cellMapping <- ""
}

# QC input data
qcExpt(expt, snakemake@params$expt_positive_column)
qcPrediction(pred_list, pred_config)
outdir <- dirname(snakemake@output$merged)
checkExistenceOfExperimentalGenesInPredictions(expt, pred_list, outdir)

# merge experimental data with predictions and write to output file
message("Merging experimentals data and predictions")
merged <- combineAllExptPred(expt = expt, 
                             pred_list = pred_list,
                             #threshold = threshold.table,
                             config = pred_config,
                             cellMapping = cellMapping, 
                             outdir = outdir,
                             fill_missing = snakemake@params$ignore_expt_missing_predictions)

# compute linear for all E-G pairs. this is useful when evaluating the distance for non-ABC
# predictors where distance is probably not part of the predictions
merged$comp.distance <- with(
  merged, abs((startPerturbationTarget + endPerturbationTarget)/2 - (startTSS + endTSS)/2)
)

# write merged data to output file
write.table(merged, snakemake@output$merged, sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)

# generate and write experimental data summary
writeExptSummary(merged, outdir)

message("All done!")
