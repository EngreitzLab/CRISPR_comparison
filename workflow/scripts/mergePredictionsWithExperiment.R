## Merge predictions with experimental data for downstream comparisons of predictions with CRISPR
## CRE perturbation data

# open log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# attach required packages and functions
merge_functions_file <- file.path(snakemake@scriptdir, "crisprComparisonMergeFunctions.R")
simple_predictors_file <- file.path(snakemake@scriptdir, "crisprComparisonSimplePredictors.R")
suppressPackageStartupMessages(source(merge_functions_file))
suppressPackageStartupMessages(source(simple_predictors_file))


## load data ---------------------------------------------------------------------------------------

# load experimental data
message("Reading experimental data in: ", snakemake@input$experiment)
expt <- fread(file = snakemake@input$experiment, showProgress = FALSE,
              colClasses = c("ValidConnection" = "character"))
message("\tLoaded experimental data with ", nrow(expt), " rows.\n")

# load prediction files
pred_files <- snakemake@input$predictions
names(pred_files) <- snakemake@params$pred_names
pred_list <- lapply(pred_files, FUN = loadPredictions, show_progress = FALSE)

# load tss and gene universe files
tss_annot <- fread(snakemake@input$tss_universe)
colnames(tss_annot) <- c("chrTSS", "startTSS", "endTSS", "gene", "score", "strandTSS")
tss_annot <- tss_annot[, -c("score", "strandTSS")]
gene_annot <- fread(snakemake@input$gene_universe)
colnames(gene_annot) <- c("chr", "start", "end", "gene", "score", "strand")

# load pred_config file
pred_config <- fread(snakemake@input$pred_config, colClasses = c("alpha" = "numeric"))

# load cell mapping file if provided
if (!is.null(snakemake@params$cell_type_mapping)) {
  cell_mappings <- lapply(snakemake@params$cell_type_mapping, FUN = fread)
  qcCellMapping(cell_mappings)
} else {
  cell_mappings <- list()
}

# QC input data
qcPredConfig(pred_config)
qcExperiment(expt, experimentalPositiveColumn = "Significant")
qcPredictions(pred_list, pred_config)

## process input data ------------------------------------------------------------------------------

# base output directory for any output
outdir <- dirname(snakemake@output$merged)

# filter experimental data for genes in gene universe
missing_file <- file.path(outdir, "expt_missing_from_gene_universe.txt")
expt <- filterExptGeneUniverse(expt, genes = tss_annot, missing_file = missing_file)

# cell type matching and filter predictions for cell types in experimental data
message("Mapping cell types in predictions to cell types in experimental data")
pred_list <- mapCellTypes(pred_list, cell_mappings = cell_mappings)
pred_list <- lapply(pred_list, FUN = function(p) p[p$ExperimentCellType %in% expt$CellType, ] )

# check if genes in experimental data are also found in predictions and write to file
# TODO: make this per cell type
genes_summary_file <- file.path(outdir, "experimental_genes_in_predictions.txt")
checkExistenceOfExperimentalGenesInPredictions(expt, pred_list, summary_file = genes_summary_file)

# merge experimental data with predictions
message("\nMerging experimentals data and predictions:")
merged <- combineAllExptPred(expt = expt, 
                             pred_list = pred_list,
                             config = pred_config,
                             outdir = outdir)

# add simple default baseline predictors
message("Adding baseline predictors:\n\tdistance to TSS\n\tnearest TSS")
dist_to_tss <- computeDistToTSS(expt)
nearest_tss <- nearestFeaturePred(expt, features = tss_annot, name = "nearestTSS")
nearest_gene <- nearestFeaturePred(expt, features = gene_annot, name = "nearestGene")
merged <- rbind(merged, dist_to_tss, nearest_tss, nearest_gene)

# rename 'CellType' column from experimental data
colnames(merged)[colnames(merged) == "CellType"] <- "ExperimentCellType"

# generate and write summary for merged data (TODO: replace by more informative output)
#merged_summary_file <- file.path(outdir, "expt_pred_merged_summary.txt")
#writeExptSummary(unique(merged[, seq_len(ncol(merged) - 4), with = FALSE]),
#                 summary_file = merged_summary_file)

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
