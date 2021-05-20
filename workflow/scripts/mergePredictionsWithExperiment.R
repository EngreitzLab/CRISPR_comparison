## Merge predictions with experimental data for downstream comparisons of predictions with CRISPR
## CRE perturbation data

# open log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# attach required packages and functions
code_file <- file.path(snakemake@scriptdir, "crisprComparisonMergeFunctions.R")
suppressPackageStartupMessages(source(code_file))


## load data ---------------------------------------------------------------------------------------

# load experimental data
message("Reading experimental data in: ", snakemake@input$experiment)
expt <- fread(file = snakemake@input$experiment, showProgress = FALSE)
message("\tLoaded experimental data with ", nrow(expt), " rows.\n")

# load prediction files
pred_files <- snakemake@input$predictions
names(pred_files) <- snakemake@params$pred_names
pred_list <- lapply(pred_files, FUN = loadPredictions, show_progress = FALSE)

# load gene universe file
genes <- fread(snakemake@input$gene_universe)

# load pred_config file
pred_config <- fread(snakemake@input$pred_config, colClasses = c("alpha" = "numeric"))

# load cell mapping file if provided
if (!is.null(snakemake@params$cell_name_mapping)) {
  cellMapping <- fread(snakemake@params$cell_name_mapping)
} else {
  cellMapping <- NULL
}

# QC input data
qcPredConfig(pred_config)
qcExpt(expt, snakemake@params$expt_positive_column)
qcPrediction(pred_list, pred_config)

## process input data ------------------------------------------------------------------------------

# base output directory for any output
outdir <- dirname(snakemake@output$merged)

# filter experimental data for genes in gene universe
missing_file <- file.path(outdir, "expt_missing_from_gene_universe.txt")
expt <- filterExptGeneUniverse(expt, genes = genes, missing_file = missing_file)

# check if genes in experimental data are also found in predictions and write to file
genes_summary_file <- file.path(outdir, "experimental_genes_in_predictions.txt")
checkExistenceOfExperimentalGenesInPredictions(expt, pred_list, summary_file = genes_summary_file)

# merge experimental data with predictions
message("\nMerging experimentals data and predictions:")
merged <- combineAllExptPred(expt = expt, 
                             pred_list = pred_list,
                             config = pred_config,
                             cellMapping = cellMapping, 
                             outdir = outdir)

# add baseline predictor of distance to TSS for all E-G pairs
message("Adding baseline predictors:\n\tdistance to TSS")
dist_to_tss <- computeDistToTSS(expt)
merged <- rbind(merged, dist_to_tss)

# generate and write summary for merged data (TODO: replace by more informative output)
merged_summary_file <- file.path(outdir, "expt_pred_merged_summary.txt")
writeExptSummary(unique(merged[, seq_len(ncol(merged) - 4), with = FALSE]),
                 summary_file = merged_summary_file)

# write merged data to main output file
outfile <- snakemake@output$merged
if (tools::file_ext(outfile) == ".gz") {
  
  # write output to gzip file
  gz <- gzfile(outfile, open = "w")
  write.table(merged, file = gz, sep = "\t", quote = FALSE, col.names = TRUE,
              row.names = FALSE)
  close(gz)
  
} else {
  
  # write output to non-compressed text file
  write.table(merged, file = outfile, sep = "\t", quote = FALSE, col.names = TRUE,
              row.names = FALSE)
  
}

message("\nAll done!")

# close log file connection
sink()
sink(type = "message")
