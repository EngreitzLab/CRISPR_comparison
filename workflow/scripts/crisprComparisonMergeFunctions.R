
# requried packages for functions in this file
library(data.table)
library(GenomicRanges)
library(dplyr)

# load a file containing a prediction set and report number of E-G pairs
loadPredictions <- function(pred_file, show_progress = FALSE) {
  message("Reading predictions in: ", pred_file)
  pred <- fread(pred_file, showProgress = show_progress)
  message("\tLoaded predictions with ", nrow(pred), " rows.\n")
  return(pred)
}

# check format of an experiment dataset and report any issues
qcExperiment <- function(expt, pos_col, remove_na_pos = FALSE) {
  
  message("Running QC on experimental data:")
  
  # make sure that all column names are valid
  illegal_cols <- c("ExperimentCellType")
  illegal_expt_cols <- intersect(colnames(expt), illegal_cols)
  if (length(illegal_expt_cols) > 0) {
    stop("Illegal columns in experiment: ", paste(illegal_expt_cols, collapse = ", "),
         call. = FALSE)
  }
  
  # check for duplicate perturbation-gene pairs
  expt_filt <- subset(expt, ValidConnection == "TRUE")
  eg_id_cols <- c("CellType", "measuredGeneSymbol", "chrom", "chromStart", "chromEnd")
  duplicates <- duplicated(expt_filt[, ..eg_id_cols])
  if (any(duplicates)) {
    stop("The experimental data file contains duplicate perturbation - gene pairs", call. = FALSE)
  }
  
  # check to make sure regulated column contains TRUE/FALSE
  pos_vals <- unique(expt[[pos_col]])
  if (!all(pos_vals %in% c(FALSE, TRUE, NA))) {
    stop("The experimental data column must contain TRUE/FALSE or 1/0", call. = FALSE)
  }
  
  # filter out any perturbation-gene pairs where pos_col is NA if specified
  if (any(is.na(pos_vals)) & remove_na_pos == TRUE) {
    na_pos_col <- is.na(expt[[pos_col]])
    message("Filtering out ", sum(na_pos_col) , " case(s) with ", pos_col, " == NA")
    expt <- expt[!na_pos_col, ]
  }
    
  # raise warning if all perturbation-gene pairs have the same pos_col value
  if (sum(!is.na(pos_vals)) == 1) {
    warning("All values in ", pos_col, " are either positives or negative", call. = FALSE)
  }
  
  message("Done")
  return(expt)
  
}

# check format of a list of predictions and report any issues
qcPredictions <- function(pred_list, pred_config, one_tss = TRUE)  {
  
  message("Running QC on predictions:")
  
  # check that there is no 'baseline' prediction set as this is used internally
  if ("baseline" %in% names(pred_list)) {
    stop("Prediction set called 'baseline' not allowed. Please rename.", call. = FALSE)
  }
  
  # make sure that minimum required columns are present
  base_cols <- c("chr", "start", "end", "TargetGene", "CellType")
  invisible(lapply(names(pred_list), FUN = check_min_cols, pred_list = pred_list,
                   pred_config = pred_config, base_cols = base_cols))
  
  # make sure all column names are valid
  illegal_cols <- c("experiment", "MappedCellType", "PredCellType")
  invisible(lapply(names(pred_list), FUN = check_illegal_cols, pred_list = pred_list,
                   illegal_cols = illegal_cols))
  
  # check that pred_col formats are ok
  invisible(lapply(names(pred_list), FUN = function(pred_name) { 
    pred <- pred_list[[pred_name]]
    conf <- pred_config[pred_config$pred_id == pred_name, ]
    invisible(mapply(FUN = check_pred_col, pred_col = conf$pred_col, boolean = conf$boolean,
                     MoreArgs = list(df = pred)))
  }))
  
  # check that predictions contain only one TSS per gene
  if (one_tss == TRUE) {
    invisible(lapply(names(pred_list), FUN = check_one_tss, pred_list = pred_list))
  }
  
  message("Done")
  
}

# check for valid pred_config file
qcPredConfig <- function(pred_config, pred_list) {
  
  message("Running QC on pred_config file")
  
  # check that all predictors in pred_list are found in the pred_config file
  missing_preds <- setdiff(names(pred_list), pred_config$pred_id)
  if (length(missing_preds) > 0) {
    stop("Not all predictors in prediction files are listed in pred_config!", call. = FALSE)
  }
  
  # check that pred_id and pred_col create a unique identifier
  config_rows <- nrow(pred_config)
  unique_ids <- nrow(unique(pred_config[, c("pred_id", "pred_col"), with = FALSE]))
  if (config_rows != unique_ids) {
    stop("pred_id and pred_col in pred_config do not provide unique identifiers for predictions.",
         call. = FALSE)
  }
  
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
  
  message("Done")
  
}

qcCellMapping <- function(cell_mappings) {
  
  dummy <- lapply(cell_mappings, FUN = function(cm) {
    dup_expt <- any(duplicated(cm$experiment))
    dup_pred <- any(duplicated(cm$predictions))
    if (any(c(dup_expt, dup_pred))) {
      stop("Currently only unique cell type mappings allowed", call. = FALSE)
    } 
  })
  
}

# map cell types from predictions to cell types in experimental data
mapCellTypes <- function(pred_list, cell_mappings) {
  
  # create empty data.table for predictions not having cell_mappings
  missing_cell_mappings <- setdiff(names(pred_list), names(cell_mappings))
  names(missing_cell_mappings) <- missing_cell_mappings
  missing_cell_mappings <- lapply(missing_cell_mappings, FUN = function(i) data.table() )
  
  # combine with provided cell mappings and sort
  cell_mappings <- c(cell_mappings, missing_cell_mappings)[names(pred_list)]
  
  # map cell types
  mapped_preds <- mapply(FUN = map_cell_type, pred = pred_list, ct_map = cell_mappings,
                         SIMPLIFY = FALSE)
  
  return(mapped_preds)
  
}

# function to map cell types for one prediction set
map_cell_type <- function(pred, ct_map) {
  
  # map cell types if cell type mapping is provided
  if (nrow(ct_map) == 0) {
    
    # simply assume that experimental cell types are identical to prediction cell types
    pred_ct <- cbind(pred, ExperimentCellType = pred$CellType)
    
  } else {
    
    # rename columns in ct_map
    colnames(ct_map)[colnames(ct_map) == "experiment"] <- "ExperimentCellType"
    
    # merge pred and cell type mapping
    pred_ct <- merge(x = pred, y = ct_map, by.x = "CellType", by.y = "predictions", all.x = TRUE,
                     all.y = FALSE)
    
    # assign new CellType if an experimental cell type was specified
    to_map <- !is.na(pred_ct[["ExperimentCellType"]])
    map_values <- pred_ct[to_map, ][["ExperimentCellType"]]
    pred_ct[to_map, "CellType"] <- map_values
    
  }
  
  # rename original CellType column in output
  pred_ct <- pred_ct[, c(colnames(pred), "ExperimentCellType"), with = FALSE]
  colnames(pred_ct)[colnames(pred_ct) == "CellType"] <- "PredictionCellType"
  
  return(pred_ct)
  
}

# check if which genes in experimental also occur in predictions
checkExistenceOfExperimentalGenesInPredictions <- function(expt, pred_list, summary_file) {
  
  # all genes in experimental data
  expt_genes <- sort(unique(expt$measuredGeneSymbol))
  
  # check which genes occur in predictions
  expt_genes_in_pred <- lapply(pred_list, function(pred) {expt_genes %in% unique(pred$TargetGene) })
  
  # create summary containing all experimental genes and their occurence in each prediction set
  summary <- data.table(experimental_genes = expt_genes, as.data.table(expt_genes_in_pred))
  
  # write to output file
  write.table(summary, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
}

# combine a list of predictions with experimental data
combineAllExptPred <- function(expt, pred_list, config, outdir, fill_pred_na = TRUE) {
  
  # merge each set of predictions with experimental data
  prediction_sets <- structure(names(pred_list), names = names(pred_list))
  merged_list <- lapply(
    prediction_sets,
    function(p) {
      combineSingleExptPred(expt = expt, pred = pred_list[[p]], pred_name = p, config = config,
                            outdir = outdir, fill_pred_na = fill_pred_na)
    })
  
  # combine merged data into one table
  output <- rbindlist(merged_list, idcol = "pred_id")
  
  # rearrange columns for output
  left_cols <- colnames(output)[seq(2, ncol(output) - 3)]
  output_col_order <- c(left_cols, "pred_id", "pred_col", "pred_value", "Prediction")
  output <- output[, output_col_order, with = FALSE]
  
  return(output)
  
}

# merge predictions with experimental data
combineSingleExptPred <- function(expt, pred, pred_name, config, outdir, fill_pred_na) {
  
  message("Overlapping predictions with experimental data for: ", pred_name)
  
  # replace any NA in predictions with fill value to avoid NAs in output
  if (fill_pred_na == TRUE) {
    pred <- fill_pred_na(pred, pred_name = pred_name, config = config)
  }
  
  # Step 1: merging overlapping enhancer - gene pairs ----------------------------------------------
  
  # subset config to specified predictions
  config_pred <- subset(config, pred_id == pred_name)
  if (nrow(config_pred) == 0) {
    stop("Predictions ", pred_name, " missing from pred_config!", call. = FALSE)
  }
  
  # subset config to prediction columns appearing in data
  config_filt <- subset(config_pred, pred_col %in% colnames(pred))
  missing_pred <- setdiff(config_pred$pred_col, config_filt$pred_col)
  if (length(missing_pred) > 0) {
    warning("Following predictor(s) specified in config file not found for ", pred_name, ": ",
            paste(missing_pred, collapse = ", "), call. = FALSE)
  }
  
  # create GenomicRanges for CRE-G links for both experimental data and predictions. this applies a
  # trick with using the seqnames to restrict overlaps to E-G pairs involving the same genes and in
  # the same cell type
  expt_gr <- with(expt, GRanges(seqnames = paste0(CellType,":", chrom,":", measuredGeneSymbol),
                                ranges = IRanges(chromStart, chromEnd)))
  pred_gr <- with(pred, GRanges(seqnames = paste0(ExperimentCellType,":", chr,":", TargetGene),
                                ranges = IRanges(start, end)))
  
  # make sure that seqnames are the same, else GRanges will report unnecessary warnings. this could
  # be removed and replaces with suppressWarnings() when calling findOverlaps() for a small gain
  # in computation time
  seqlevels_all_pairs <- as.character(unique(c(seqnames(expt_gr), seqnames(pred_gr))))
  seqlevels(expt_gr) <- seqlevels_all_pairs
  seqlevels(pred_gr) <- seqlevels_all_pairs
  
  # find overlaps between predictions and experimental data
  ovl <- findOverlaps(expt_gr, pred_gr)
  
  # merge predictions with experimental data
  pred_merge_cols <- c("PredictionCellType", config_filt$pred_col)  # columns in predictions to add
  merged <- cbind(expt[queryHits(ovl)], pred[subjectHits(ovl), pred_merge_cols, with = FALSE])
  
  # Step 2: aggregating pairs with multiple overlaps -----------------------------------------------
  
  # sometimes perturbed elements will overlap multiple predicted elements (eg in the case of a large
  # deletion). in these cases need to summarize, e.g., sum ABC.Score across model elements
  # overlapping the deletion. this requires a config file describing how each prediction column
  # should be aggregated
  agg_cols <- setdiff(colnames(merged), config_filt$pred_col)
  merged <- collapseEnhancersOverlappingMultiplePredictions(merged, config = config_filt,
                                                            agg_cols = agg_cols)
  
  # Step 3: experimental data missing predictions --------------------------------------------------
  
  # A tested enhancer element may not have a prediction. For ABC this is typically the case if the
  # tested element does not overlap a DHS peak. In this case we need to fill the predictions table
  
  # add 'Prediction' column to merged with value 1, indicating pairs that overlap with predictions
  merged$Prediction <- 1
  
  # get pairs from experimental data that are missing from predictions 
  expt_missing_preds <- expt[setdiff(seq_len(nrow(expt)), queryHits(ovl)), ]
  
  # write these to a text file in the output sub-directory
  expt_missing_pred_dir <- file.path(outdir, "experimentalDataMissingPredictions")
  dir.create(expt_missing_pred_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    x = expt_missing_preds[, -c("chrTSS", "startTSS", "endTSS")],
    file = file.path(expt_missing_pred_dir, paste0(pred_name, "_expt_missing_predictions.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  # fill in missing values
  expt_missing_preds$PredictionCellType <- NA_character_
  expt_missing_preds <- fillMissingPredictions(expt_missing_preds, config = config_filt,
                                               agg_cols = agg_cols)
  
  # add 'Prediction' column with value 0, indicating pairs that were not found in predictions
  expt_missing_preds$Prediction <- 0
  
  # merge filled data with merged data
  merged <- rbind(merged, expt_missing_preds[, colnames(merged), with = FALSE])
  
  # Step 4: create output --------------------------------------------------------------------------
  
  ### WHY IS THIS NECESSARY??
  #score_colname <- colnames(merged)[grep("Score", colnames(merged))]
  #merged[[paste0(score_colname,'.Prediction')]] <- merged[[score_colname]] * merged$Prediction
  
  # get distance for predictions that overlap experiments
  #dist_colname <- colnames(merged)[grep("distance", colnames(merged))]
  #merged[[paste0(dist_colname,'.Prediction')]] <- merged[[dist_colname]] * merged$Prediction
  
  # rename new columns based on prediction dataset name
  #pred_cols <- colnames(merged) %in% c(config_filt$pred_col, "Prediction")
  #colnames(merged)[pred_cols] <- paste0(pred_name, ".", colnames(merged)[pred_cols])
  
  # convert to long format to generate output
  output <- melt(merged, measure.vars =  config_filt$pred_col, variable.name = "pred_col",
                 value.name = "pred_value")
  
  # sort output according to genomic coordinates of enhancers and target gene
  sortcols <- c("chrom", "chromStart", "chromEnd", "measuredGeneSymbol")
  setorderv(output, sortcols)
  
  return(output)
  
}

# aggregate experimental enhancers overlapping multiple predicted enhancers
collapseEnhancersOverlappingMultiplePredictions <- function(df, config, agg_cols) {
  
  # summarize columns as defined in config
  all_list <- mapply(FUN = function(pred_col, agg_func) {
    agg_func <- get(agg_func)  # get function from string
    df[, setNames(.(agg_func(get(pred_col))), pred_col), by = agg_cols]
  }, pred_col = config$pred_col, agg_func = config$aggregate_function, SIMPLIFY = FALSE)
  
  # special handling for aggregating the class column
  class_agg <- function(x) {
    if ("promoter" %in% x) {
      return("promoter")
    } else if ("tss" %in% x) {
      return("tss")
    } else if ("genic" %in% x) {
      return("genic")
    } else if ("distal" %in% x) {
      return("distal")
    } else if ("intergenic" %in% x) {
      return("intergenic")
    } else {
      return("UNKNOWN")
    }
  }
  
  # TODO: implement if class from predictions is added during merging
  #if ("class" %in% colnames(df)) {
  #  class_temp <- aggregate(df$class, by = list_for_agg, FUN = class_agg)
  #  colnames(class_temp)[colnames(class_temp) == "x"] <- "class"
  #  all_list$class <- class_temp
  #}
  
  # merge all the aggregates together to make collapsed data.frame
  # TODO: AVOID INTERMEDIATE DATA.FRAME BY Reduce
  output <- Reduce(function(df1, df2) merge(df1, df2, by = agg_cols), all_list)
  output <- as.data.table(output)
  
  return(output)
  
}

# fill in prediction values for experimental pairs missing from prediction data
fillMissingPredictions <- function(df, config, agg_cols) {
  
  # fill in missing predictions as described in the config file
  for (i in seq_len(nrow(config))) {
    df[, config$pred_col[[i]]] <- config$fill_value[i]
  }
  
  # fill in unknown columns (columns appearing in agg_cols, but not experiment or pred_cols)
  # unk_cols <- setdiff(c("class", agg_cols), unique(c(colnames(df), config$pred_col)))
  # df[, unk_cols] <- "Merge:UNKNOWN"
  
  return(df)
  
}

# filter experimental data for genes in gene universe and add TSS coordinates to experimental data
filterExptGeneUniverse <- function(expt, genes, missing_file = NULL) {
  
  # remove any existing TSS annotations from expt data
  expt_cols <- colnames(expt)
  tss_cols <- expt_cols %in% c("chrTSS", "startTSS", "endTSS")
  expt <- expt[, !tss_cols, with = FALSE]
  
  # add gene TSSs from specified gene universe to experimental data
  expt <- merge(expt, genes, by.x = "measuredGeneSymbol", by.y = "gene", all.x = TRUE, sort = FALSE)
  expt <- expt[, ..expt_cols]
  
  # get experimental data from genes that do not appear in the gene universe
  expt_missing <- expt[is.na(expt$chrTSS), ]
  
  # if there are any, write these to output file (if specified)
  if (nrow(expt_missing) > 0) {
    message("cre-gene pairs in experimental data missing from gene universe: ", nrow(expt_missing))
    if (!is.null(missing_file)) {
      write.table(expt_missing, file = missing_file, quote = FALSE, row.names = FALSE, sep = "\t")
    }
  }
  
  # filter out missing expt data, arrange columns and return as output
  expt_filt <- expt[!is.na(expt$chrTSS), ]
  
  return(expt_filt)
  
}

# create summary of experimental data
writeExptSummary <- function(df, summary_file) {
  
  # create summary table
  df_summary <- as.data.frame(list(
    numConnections = nrow(df),
    numIncludeInModel = sum(df$IncludeInModel),
    numIncludeInModelRegulated = sum(df$IncludeInModel & df$Regulated)
  ))
  
  # write to specified file
  write.table(df_summary, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
}

## HELPER FUNCTIONS ================================================================================

# check input --------------------------------------------------------------------------------------

# make sure that minimum required columns are present in a prediction set
check_min_cols <- function(pred_list, pred_config, pred, base_cols) {
  score_col <- pred_config[pred_config$pred_id == pred, ][["pred_col"]]  # get score cols for pred
  missing_cols <- setdiff(c(base_cols, score_col), colnames(pred_list[[pred]]))
  if (length(missing_cols) > 0) {
    stop("Missing columns in predictions '", pred, "': ", paste(missing_cols, collapse = ", "), ".",
         call. = FALSE)
  }
}

# check if a prediction set contains illegal column names
check_illegal_cols <- function(pred_list, pred, illegal_cols) {
  df <- pred_list[[pred]]
  wrong_pred_cols <- intersect(colnames(df), illegal_cols)
  if (length(wrong_pred_cols) > 0) {
    stop("Illegal columns in predictions '", pred, "': ", paste(wrong_pred_cols, collapse = ", "),
         ".", call. = FALSE)
  }
}

# check if a given predictor column has the correct format
check_pred_col <- function(df, pred_col, boolean) {
  pred_values <- df[[pred_col]]
  if (boolean == TRUE) {
    if (any(!unique(as.numeric(pred_values)) %in% c(0, 1))) {
      stop("Incorrect format for boolean predictor. Must be either 1/0 or TRUE/FALSE.",
           call. = FALSE)
    }
  } else {
    if (!is.numeric(pred_values)) {
      stop("Continuous predictors must have numeric values.", call. = FALSE)
    }
  }
}

# check that a prediction set only uses one TSS per gene, i.e. each enhancer-gene pair occurs once
check_one_tss <- function(pred_list, pred) {
  df <- pred_list[[pred]]
  total_pairs <- nrow(df)
  unique_pairs <- nrow(unique(df[, c("chr", "start", "end", "TargetGene")]))
  if (unique_pairs < total_pairs) {
    stop("Predictions '", pred, "' uses more than one TSS per gene.", call. = FALSE)
  }
}

# replace NAs in prediction scores with appropriate fill values
fill_pred_na <- function(pred, pred_name, config) {
  
  # get fill values for columns containing prediction scores in pred
  fill_values <- config[config$pred_id == pred_name, c("pred_col", "fill_value")]
  fill_values <- as.list(structure(fill_values$fill_value, names = fill_values$pred_col))
  
  # convert any integer scores to numeric to avoid errors when replacing NAs
  pred <- mutate(pred, across(where(is.integer) & names(fill_values), as.numeric)) 
  
  # replace NAs with fill values
  pred <- tidyr::replace_na(pred, replace = fill_values)
  
  return(pred)
  
}

## Deprecated import functions =====================================================================

# fread_ignore_empty <- function(f, ...) {
#   tryCatch({
#     return(fread(f, fill = TRUE, ...))
#   }, error = function(e){
#     print(paste0("Could not open file: ", f))
#     return()
#   })
# }
# 
# fread_gz_ignore_empty <- function(f, ...) {
#   tryCatch({
#     return(fread(paste0("gunzip -c ", f), ...))
#   }, error = function(e){
#     print(paste0("Could not open file: ", f))
#     return()
#   })
# }
# 
# smart_fread <- function(f, ...) {
#   if (summary(file(f))$class == "gzfile") {
#     out <- fread_gz_ignore_empty(f, ...)
#   } else {
#     out <- fread_ignore_empty(f, ...)
#   }
# 
#   #con <- file(f)
#   #on.exit(close(con), add = TRUE)
#   tryCatch({
#     closeAllConnections()
#   }, error = function(e) {
#     print(e)
#   }
#   )
# 
#   return(out)
# }
# 
# loadDelimited <- function(file.list) {
#   data.list <- lapply(file.list, smart_fread)
#   return(rbindlist(data.list, fill = TRUE))
# }
# 
# loadFileString <- function(file.str, delim = ",") {
#   file.list <- strsplit(file.str, split = delim)[[1]]
#   return(loadDelimited(file.list))
# }
# 
# loadPredictions <- function(pred.table) {
#   pred.list <- lapply(pred.table$path, function(s) {
#     message("Loading dataset: ", s)
#     df <- loadFileString(s)
#     message("\tDataset loaded with ", nrow(df), " rows")
#     return(df)
#     })
#   names(pred.list) <- pred.table$name
#   return(pred.list)
# }
# 
# # Deprecated TODO: re-implement with new prediction format
# qcPrediction <- function(pred_list, pred_config)  {
#   # Ensure that the fill value for each prediction column is at the extreme end of its range
#   message("Running QC on predictions")
#   
#   doOnePred <- function(pred, config) {
#     this_cols <- intersect(colnames(pred), config$pred_col)
#     lapply(this_cols, function(s) {
#       qcCol(col_name = s, 
#             colData = pred[, ..s],
#             fill_val = subset(config, pred_col == s)$fill_val, 
#             isInverted = subset(config, pred_col == s)$inverse_predictor)
#     })
#   }
#   
#   qcCol <- function(col_name, colData, fill_val, isInverted) {
#     # For each prediction column check that its missing fill val is at the extreme end of its range
#     isBad <- (isInverted & fill_val < pmin(colData)) | (!isInverted & fill_val > pmin(colData))
#     suppressWarnings(if (isBad) stop(paste0("Fill val for column ", col_name, " is not at the extreme of its range!", fill_val, " ", pmin(colData))))
#   }
#   
#   dummy <- lapply(pred_list, function(s) doOnePred(s, config = pred_config))
#   message("Done")
# }
