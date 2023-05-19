## Merge predictions with experimental data for downstream comparisons of predictions with CRISPR
## CRE perturbation data

# save.image("merge.rda")
# stop()

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

# load pred_config file
pred_config <- fread(snakemake@input$pred_config,
                     colClasses = c("alpha" = "numeric", "color" = "character"))

# config entry for this comparison is used to load named list of input files
config <- snakemake@config$comparisons[[snakemake@wildcards$comparison]]

# load experimental data
message("Reading experimental data in: ", snakemake@input$experiment)
expt <- fread(file = snakemake@input$experiment, showProgress = FALSE,
              colClasses = c("ValidConnection" = "character"))
message("\tLoaded experimental data with ", nrow(expt), " rows.\n")

# load prediction files
pred_files <- config$pred
pred_list <- lapply(pred_files, FUN = loadPredictions, show_progress = FALSE)

# load tss and gene universe files
tss_annot <- fread(snakemake@input$tss_universe)
colnames(tss_annot) <- c("chrTSS", "startTSS", "endTSS", "gene", "score", "strandTSS")
tss_annot <- tss_annot[, -c("score", "strandTSS")]
gene_annot <- fread(snakemake@input$gene_universe)
colnames(gene_annot) <- c("chr", "start", "end", "gene", "score", "strand")

# load optional cell mapping files if provided
ct_map_files <- config$cell_type_mapping
if (!is.null(ct_map_files)) {
  cell_mappings <- lapply(ct_map_files, FUN = fread)
  qcCellMapping(cell_mappings)
} else {
  cell_mappings <- list()
}

# load expressed genes files if provided
if (!is.null(snakemake@input$expressed_genes)) {
  expressed_genes <- fread(snakemake@input$expressed_genes)
} else {
  expressed_genes <- NULL
}

# QC pred_config file
qcPredConfig(pred_config, pred_list = pred_list)

# QC predictions and experimental data
qcPredictions(pred_list, pred_config = pred_config, one_tss = FALSE)
expt <- qcExperiment(expt, pos_col = snakemake@params$pos_col, remove_na_pos = TRUE)

## process input data ------------------------------------------------------------------------------

# base output directory for any output
outdir <- dirname(snakemake@output$merged)

# filter experimental data for genes in gene universe
missing_file <- file.path(outdir, "expt_missing_from_gene_universe.txt")
expt <- filterExptGeneUniverse(expt, genes = tss_annot, missing_file = missing_file)

# add expression information to experimental data if specified
if (!is.null(snakemake@input$expressed_genes)) {
  expt <- addGeneExpression(expt, expressed_genes = expressed_genes)
}

# cell type matching and filter predictions for cell types in experimental data
message("Mapping cell types in predictions to cell types in experimental data")
pred_list <- mapCellTypes(pred_list, cell_mappings = cell_mappings)
pred_list <- lapply(pred_list, FUN = function(p) p[p$ExperimentCellType %in% expt$CellType, ] )

# verify if bad cell matching resulted in no data for some predictions after matching
pred_rows <- vapply(pred_list, FUN = nrow, FUN.VALUE = integer(1))
if (any(pred_rows == 0)) {
  stop("No predictions left for ", paste(names(pred_rows[pred_rows == 0]), collapse = ", "),
       " after cell type matching. Check that cell type mapping files are correct.", call. = FALSE)
}

## overlap experimental data with predictions ------------------------------------------------------

# check if genes in experimental data are also found in predictions and write to file
# TODO: make this per cell type
genes_summary_file <- file.path(outdir, "experimental_genes_in_predictions.txt")
checkExistenceOfExperimentalGenesInPredictions(expt, pred_list, summary_file = genes_summary_file)

# merge experimental data with predictions
message("\nMerging experimentals data and predictions:")
merged <- combineAllExptPred(expt = expt, 
                             pred_list = pred_list,
                             config = pred_config,
                             outdir = outdir,
                             fill_pred_na = TRUE)

## compute baseline predictors ---------------------------------------------------------------------

# get simple baseline predictors to compute
baseline_pred_ids <- c("distToTSS", "nearestTSS", "nearestGene", "within100kbTSS")
if (!is.null(expressed_genes)) {
  baseline_pred_ids <- c(baseline_pred_ids,
                         c("nearestExprTSS", "nearestExprGene", "within100kbExprTSS"))
}

# compute and add baseline predictors
message("Computing baseline predictors:\n\t", paste(baseline_pred_ids, collapse = "\n\t"))
baseline_preds <- computeBaselinePreds(expt, preds = baseline_pred_ids, tss_annot = tss_annot,
                                       gene_annot = gene_annot, expressed_genes = expressed_genes)
merged <- rbind(merged, baseline_preds)

## format output and write to file -----------------------------------------------------------------

# rename 'CellType' column from experimental data
colnames(merged)[colnames(merged) == "CellType"] <- "ExperimentCellType"

# write merged data to main output file
readr::write_tsv(merged, file = snakemake@output$merged)

message("\nAll done!")

# close log file connection
sink()
sink(type = "message")
