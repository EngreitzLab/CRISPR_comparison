## Merge predictions with experimental data for downstream comparisons of predictions with CRISPR
## CRE perturbation data

# save.image("merge.rda")
# stop()

# open log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# attach required packages and functions
suppressPackageStartupMessages({
  source(file.path(snakemake@scriptdir, "crisprComparisonLoadInputData.R"))
  source(file.path(snakemake@scriptdir, "crisprComparisonMergeFunctions.R"))
  source(file.path(snakemake@scriptdir, "crisprComparisonSimplePredictors.R"))
})

## load data ---------------------------------------------------------------------------------------

# directory for all output
outdir <- dirname(snakemake@output$merged)

# config entry for this comparison is used to load named list of input files
config <- snakemake@config$comparisons[[snakemake@wildcards$comparison]]

# load pred_config file
include_col <- ifelse(is.null(snakemake@params$include_col), "include", snakemake@params$include_col)
pred_config <- importPredConfig(snakemake@input$pred_config,
                                expr = !is.null(snakemake@input$expressed_genes),
                                include_col = include_col,
                                filter = snakemake@params$filter_include_col)

# load experimental data
message("Reading CRISPR data in: ", snakemake@input$experiment)
expt <- fread(file = snakemake@input$experiment, showProgress = FALSE,
              colClasses = c("ValidConnection" = "character"))
message("\tLoaded CRISPR data for ", nrow(expt), " E-G pairs\n")

# load tss and gene universe files
tss_annot <- fread(snakemake@input$tss_universe, select = 1:6,
                   col.names = c("chrTSS", "startTSS", "endTSS", "gene", "score", "strandTSS"))
gene_annot <- fread(snakemake@input$gene_universe, select = 1:6,
                    col.names = c("chr", "start", "end", "gene", "score", "strand"))

# load all prediction files
pred_list <- loadPredictions(config$pred, show_progress = FALSE)

# if specified, filter out any predictions where elements overlap annotated gene TSS
if (snakemake@params$filter_tss == TRUE) {
  tss_filt_file <- file.path(outdir, "filter_predictions_tss_results.txt")
  pred_list <- filterPredictionsTSS(pred_list, tss_annot = tss_annot, summary_file = tss_filt_file)
}

# combined files per predictor, if files for multiple cell types were provided
pred_list <- lapply(pred_list, FUN = rbindlist)

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
  expressed_genes <- loadGeneExpr(snakemake@input$expressed_genes)
} else {
  expressed_genes <- NULL
}

# QC pred_config file
qcPredConfig(pred_config, pred_list = pred_list)

# QC predictions and experimental data
pred_list <- qcPredictions(pred_list, pred_config = pred_config, one_tss = FALSE)
expt <- qcExperiment(expt, pos_col = snakemake@params$pos_col, remove_na_pos = TRUE)

## process input data ------------------------------------------------------------------------------

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

# get all simple baseline predictors to compute
baseline_pred_ids <- c("distToTSS", "distToGene", "nearestTSS", "nearestGene", "within100kbTSS")
if (!is.null(expressed_genes)) {
  baseline_pred_ids <- c(baseline_pred_ids,
                         c("nearestExprTSS", "nearestExprGene", "within100kbExprTSS"))
}

# only retain baseline predictors to include in benchmark
baseline_pred_ids <- intersect(baseline_pred_ids, pred_config$pred_col)

# compute and add baseline predictors
message("Computing baseline predictors:\n\t", paste(baseline_pred_ids, collapse = "\n\t"))
baseline_preds <- computeBaselinePreds(expt, preds = baseline_pred_ids, tss_annot = tss_annot,
                                       gene_annot = gene_annot, expressed_genes = expressed_genes)
merged <- rbind(merged, baseline_preds)

## write to file -----------------------------------------------------------------------------------

# write merged data to main output file
fwrite(merged, file = snakemake@output$merged, sep = "\t", quote = FALSE, na = "NA")

message("\nAll done!")

# close log file connection
sink()
sink(type = "message")
