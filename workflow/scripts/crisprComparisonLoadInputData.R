## Functions to load, process and QC benchmarking pipeline input data files

library(data.table)
library(dplyr)


#' Import predictor config
#' 
#' Import predictor config file specifying parameters for each predictor in the benchmark and add
#' all required columns and baseline predictors.
#' 
#' @param pred_config_file Path to the predictor config file.
#' @param expr (logical) Should baseline predictors use information on whether genes are expressed?
#'   Can later be switched using `toggleBaselinePredictors`.
#' @param include_col Column containing information on which predictors to include in benchmark
#'   (TRUE/FALSE). Default column name is 'include', which is added if not in predictor config file
#'   with setting all predictors to be included in benchmark.
#' @param filter (logical) Specifying whether predictor config file should be filtered for
#'   predictors to include in benchmark only (default: TRUE).
importPredConfig <- function(pred_config_file, expr = FALSE, include_col = "include",
                             filter = TRUE) {
  
  # classes of standard columns in pred_config file
  cols <- c("pred_id" = "character", "pred_col" = "character", "boolean" = "logical",
            "alpha" = "numeric",  "aggregate_function" = "character", "fill_value" = "numeric",
            "inverse_predictor" = "logical", "pred_name_long" = "character", "color" = "character")
  
  # load pred_config file and subset to required columns
  pred_config <- fread(pred_config_file, colClasses = cols)
  select_cols <- c(names(cols), intersect(include_col, colnames(pred_config)))
  pred_config <- pred_config[, ..select_cols]
  
  # rename include_col or add default if missing from pred_config 
  if (include_col %in% colnames(pred_config)) {
    colnames(pred_config)[colnames(pred_config) == include_col] <- "include"
  } else {
    pred_config$include <- TRUE
  }
  
  # add unique identifier (pred_uid) for each predictor
  pred_config$pred_uid <- paste(pred_config$pred_id, pred_config$pred_col, sep = ".")
  
  # create default baseline predictors config if needed
  baseline_config <- create_baseline_pred_config(expr)
  
  # if no colors are specified for any of the predictors to include in the benchmark, set colors
  # of all default baseline predictors to NA as well
  include_preds <- pred_config[pred_config$include == TRUE, ]
  if (all(is.na(include_preds$color)) == TRUE) {
    baseline_config$color <- NA_character_
  } 
  
  # config for baseline predictors can be set in pred_config, so only add default values for
  # baseline predictors not specified in pred config
  add_baseline_preds <- setdiff(baseline_config$pred_uid, pred_config$pred_uid)
  pred_config <- rbind(pred_config, subset(baseline_config, pred_uid %in% add_baseline_preds))
  
  # check that generated pred_uid and long predictor names are unique TODO: MOVE TO QC FUNCTION?
  check_unique_identifier(pred_config, col = "pred_uid")
  check_unique_identifier(pred_config, col = "pred_name_long")
  
  # if specified filter out any predictors not included in benchmark
  if (filter == TRUE) {
    pred_config <- pred_config[pred_config$include == TRUE, ]
  }
  
  return(pred_config)
  
}


#' Load predictions
#' 
#' Load input prediction files and create a list of tables, containing data for all predictors.
#' 
#' @param pred_files Named list containing all input prediction files.
#' @param load_function Function used to load prediction files (default = 'fread'). Used to provide
#'   functions to load custom file formats.
#' @param show_progress (logical) Should detailed loading progress be showed?
loadPredictions <- function(pred_files, load_function = fread, show_progress = FALSE) {
  
  message("Loading prediction files:")
  
  pred_ids <- structure(names(pred_files), names = names(pred_files))
  preds <- lapply(pred_ids, FUN = load_pred_files_predictor, pred_files = pred_files,
                  load_function = load_function, show_progress = show_progress)

  return(preds)
  
}


#' Toggle baseline predictors on or off for benchmarks
#' 
#' Set internal baseline predictors on predictor config table on or off to include in benchmarks.
#' 
#' @param pred_config Table containing predictor config information.
#' @param on,off Vectors specifying which baseline predictors should be turned on or off
toggleBaselinePredictors <- function(pred_config, on = NULL, off = NULL) {
  
  # parse input arguments
  preds <- c("distToTSS", "distToGene", "nearestTSS", "nearestGene", "within100kbTSS",
             "nearestExprTSS", "nearestExprGene", "within100kbExprTSS")
  if(!is.null(on))  on  <- match.arg(on,  choices = preds, several.ok = TRUE)
  if(!is.null(off)) off <- match.arg(off, choices = preds, several.ok = TRUE)
  
  # toggle on of off specific baseline predictors
  pred_config[pred_config$pred_col %in% on,  "include"] <- TRUE
  pred_config[pred_config$pred_col %in% off, "include"] <- FALSE
  
  return(pred_config)
  
}


#' Load information on expressed genes
#' 
#' Load files containing information on which genes are expressed in different cell types.
#' 
#' @param ... Paths to one or more input files containing information on which genes are expressed.
#'   If more than one file are passed, these are simply concatenated.
loadGeneExpr <- function(...) {
  expr <- rbindlist(lapply(..., FUN = fread))
  return(expr)
}

## Helper functions --------------------------------------------------------------------------------

# add baseline predictors to pred config table
create_baseline_pred_config <- function(expr) {
  
  # all internal baseline identifiers
  baseline_preds <- c("distToTSS", "distToGene", "nearestTSS", "nearestGene", "within100kbTSS",
                      "nearestExprTSS", "nearestExprGene", "within100kbExprTSS")
  baseline_preds_uid <- paste("baseline", baseline_preds, sep = ".")
  
  # create pred_config table 
  baseline_config <- data.table(
    pred_id = "baseline",
    pred_col = baseline_preds,
    boolean = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    alpha = c(NA_real_, 1, 1, 1, 1, 1, 1, 1),
    aggregate_function = c("mean", "max", "max", "max", "max", "max", "max", "max"),
    fill_value = c(Inf, 0, 0, 0, 0, 0, 0, 0),
    inverse_predictor = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    pred_name_long = c("Distance to TSS", "Distance to Gene", "Nearest TSS", "Nearest gene",
                       "Within 100kb of TSS", "Nearest expr. TSS", "Nearest expr. gene",
                       "Within 100kb of expr. TSS"),
    color = c("#ffa600", "#ffa600", "#595959", "#bebebe", "#000000", "#595959", "#bebebe",
              "#000000"),
    include = c(TRUE, FALSE, rep(NA, times = 6)),
    pred_uid = baseline_preds_uid
  )
  
  # set default include column values based on expr argument, which specifies whether baseline
  # predictors based on expressed genes should be included or not
  expr_baselines <- c("nearestExprTSS", "nearestExprGene", "within100kbExprTSS")
  non_expr_baselines <- c("nearestTSS", "nearestGene", "within100kbTSS")
  if (expr == TRUE) {
    baseline_config <- toggleBaselinePredictors(baseline_config, on = expr_baselines,
                                                off = non_expr_baselines)
  } else {
    baseline_config <- toggleBaselinePredictors(baseline_config, on = non_expr_baselines,
                                                off = expr_baselines)
  }
  
  return(baseline_config)
  
}

# check that a given column is a unique identifier in a predictor config file
check_unique_identifier <- function(pred_config, col) {
  
  # get specified column for all predictors to include in benchmark
  ids <- pred_config[pred_config$include == TRUE, ][[col]]
  
  # raise error if ids in specified column are not unique identifiers
  if (any(table(ids) > 1)) {
    stop("'", col, "' in pred_config_file is not a unique identifier.", call. = FALSE)
  }
  
}

# load and concatenate all prediction files for one predictor provided via pred_id
load_pred_files_predictor <- function(pred_id, pred_files, load_function, show_progress) {
  
  message("\tLoading predictions for: ", pred_id)
  
  # load all files for the given predictor
  preds <- lapply(pred_files[[pred_id]], FUN = function(file, show_progress) {
    message("\t\tReading: ", file)
    load_function(file, showProgress = show_progress)
  }, show_progress = show_progress)
  
  # combine into one table
  preds <- rbindlist(preds)
  
  message("\tLoaded predictions with a total of ", nrow(preds), " rows.")
  
  return(preds)
}

# load one predictor in igvf format and convert to generic format
load_igvf_pred_file <- function(file, showProgress) {
  
  # get cell type from header
  header <- grep(readLines(file, n = 1000), pattern = "^#", value = TRUE)
  cell_type <- grep(header, pattern = "BiosampleString", value = 1)
  cell_type <- sub(".*BiosampleString:[ ]*", "", cell_type)
  
  # load predictions table and add cell type
  pred <- fread(file, skip = length(header), showProgress = showProgress)
  pred$CellType <- cell_type
  
  # select relevant columns for benchmarking pipeline
  base_cols <- c("ElementChr", "ElementStart", "ElementEnd", "ElementName", "ElementClass",
                 "ElementStrand", "GeneSymbol", "GeneEnsemblID", "GeneTSS", "CellType")
  score_cols <- setdiff(colnames(pred), base_cols)
  select_cols <- c(base_cols[c(1:3, 7, 10)], score_cols)
  pred <- pred[, ..select_cols]
  
  # rename column to CRISPR benchmark internal columns
  colnames(pred)[1:4] <- c("chr", "start", "end", "TargetGene")
  
  return(pred)
  
}
