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
  pred_config <- unique(pred_config[, ..select_cols])
  
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

  # if specified filter out any predictors not included in benchmark
  if (filter == TRUE) {
    pred_config <- pred_config[pred_config$include == TRUE, ]
  }
  
  # check that generated pred_uid names are unique and all colors are valid color names/codes
  check_unique_identifier(pred_config, col = "pred_uid")
  check_colors(pred_config)
  
  return(pred_config)
  
}


#' Load predictions
#' 
#' Load input prediction files and create a list of tables, containing data for all predictors.
#' 
#' @param pred_files Named list containing all input prediction files.
#' @param show_progress (logical) Should detailed loading progress be showed?
loadPredictions <- function(pred_files, show_progress = FALSE) {
  
  message("Loading prediction files:")
  
  pred_ids <- structure(names(pred_files), names = names(pred_files))
  preds <- lapply(pred_ids, FUN = load_pred_files_predictor, pred_files = pred_files,
                  show_progress = show_progress)

  return(preds)
  
}

#' Filter predictions for TSS overlaps
#' 
#' Filter out any predictions, where elements overlap annotated gene TSS
#' 
#' @param pred_list List containing predictions for all predictors.
#' @param tss_annot Table containing TSS annotations from tss_universe file.
#' @param summary_file Path to file where number of filtered out elements and E-G pairs will be
#'  reported.
filterPredictionsTSS <- function(pred_list, tss_annot, summary_file) {
  
  # create a GRanges object containing TSS coordinates for filtering
  tss <- makeGRangesFromDataFrame(tss_annot, seqnames.field = "chrTSS", start.field = "startTSS",
                                  end.field = "endTSS", starts.in.df.are.0based = TRUE)
  
  # filter all predictions based on TSS overlaps
  message("Filtering out predictions overlapping gene TSS for:")
  pred_ids <- structure(names(pred_list), names = names(pred_list))
  pred_list <- lapply(pred_ids, FUN = function(x, pred_list, tss) {
    message("\t", x)
    lapply(pred_list[[x]], FUN = filter_pred_tss, tss = tss)
  }, pred_list = pred_list, tss = tss)
  
  # extract TSS filtering results and combined into one table
  tss_filter <- lapply(pred_list, FUN = lapply, `[[`, 2)
  tss_filter <- rbindlist(lapply(tss_filter, FUN = rbindlist, idcol = "file"), idcol = "pred_id")
  
  # write TSS filtering stats to summary file
  fwrite(tss_filter, file = summary_file, sep = "\t", quote = FALSE, na = "NA")
  
  # remove filter from 'preds' list
  pred_list <- lapply(pred_list, FUN = lapply, `[[`, 1)
  
  return(pred_list)
  
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
  if (any(table(pred_config[[col]]) > 1)) {
    stop("'", col, "' in pred_config_file is not a unique identifier.", call. = FALSE)
  }
}

# check if colors in pred_config file are valid color ids
check_colors <- function(pred_config) {
 
   # check that colors are valid options
  valid_colors <- vapply(pred_config$color, FUN = function(col) {
    tryCatch(expr = is.matrix(col2rgb(col)),
             error = function(err) return(FALSE))
  }, FUN.VALUE = logical(1))
  
  # raise error if invalid color specification was found
  invalid_colors <- names(valid_colors[valid_colors == FALSE])
  if (length(invalid_colors > 0)) {
    stop("Invalid color specification(s): ", paste(invalid_colors, collapse = ", "), call. = FALSE)
  }
  
}

# load and concatenate all prediction files for one predictor provided via pred_id
load_pred_files_predictor <- function(pred_id, pred_files, show_progress) {
  
  message("\tLoading predictions for: ", pred_id)
  
  # load all files for the given predictor
  files <- structure(pred_files[[pred_id]], names = pred_files[[pred_id]])
  preds <- lapply(files, FUN = function(file, show_progress) {
    message("\t\tReading: ", file)
    fread(file, showProgress = show_progress)
  }, show_progress = show_progress)
  
  # report the total number of loaded E-G pairs
  n_pairs <- sum(vapply(preds, FUN = nrow, FUN.VALUE = integer(1)))
  message("\t\tLoaded predictions for ", n_pairs, " E-G pairs")
  
  return(preds)
}

# filter one prediction file for elements not overlapping any provided TSS 
filter_pred_tss <- function(pred, tss) {
  
  # get unique element coordinates and names for predictions
  element_cols <- c(colnames(pred)[[1]], c("start", "end", "name"))  # 1st col can be 'chr' / '#chr'
  elements <- unique(pred[, ..element_cols])
  
  # create GRanges object for both elements and TSSs
  elements <- makeGRangesFromDataFrame(elements, seqnames.field = element_cols[[1]],
                                       keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  
  # filter out elements that are overlapping any TSS
  elements_filt <- subsetByOverlaps(elements, tss, invert = TRUE)
  pred_filt <- pred[pred$name %in% elements_filt$name, ]
  
  # get the number of filtered out elements and predictions
  filter <- data.table(removed_elements = length(elements) - length(elements_filt),
                       removed_predictions = nrow(pred) - nrow(pred_filt))
  
  return(list(pred_filt, filter))
  
}
