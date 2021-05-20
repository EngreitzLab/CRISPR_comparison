library(data.table)
library(GenomicRanges)

loadPredictions <- function(pred_file, show_progress = FALSE) {
  message("Reading predictions in: ", pred_file)
  pred <- fread(pred_file, showProgress = show_progress)
  message("\tLoaded predictions with ", nrow(pred), " rows.\n")
  return(pred)
}

qcExpt <- function(expt, experimentalPositiveColumn) {
  message("Running QC on experimental data")
  expt <- subset(expt, IncludeInModel)
  #Check for duplicate experiments
  dupe <- any(duplicated(expt[, c("CellType","GeneSymbol","chrPerturbationTarget","startPerturbationTarget","endPerturbationTarget")] ))
  if (dupe) {
    stop("The experimental data file contains duplicate experiments!", call. = FALSE)
  }
  
  #check to make sure regulated column contains TRUE/FALSE
  # reg.vals <- sort(unique(expt[, get(experimentalPositiveColumn)]))
  # if (!(identical(reg.vals, c(FALSE, TRUE)) | identical(reg.vals, c(0, 1)))) {
  #   print("Error: The experimental data column must contain exactly two distinct values: TRUE and FALSE")
  #   stop()
  # }

  #check to make sure regulated column contains TRUE/FALSE
  reg.vals <- sort(unique(expt[, get(experimentalPositiveColumn)]))
  if (!(all(reg.vals %in% c(FALSE, TRUE)) | all(reg.vals %in% c(0, 1)))) {
    stop("The experimental data column must contain TRUE/FALSE", call. = FALSE)
  }
  if (length(reg.vals) == 1) {
    warning("All values are either positives or negatives. Plotting code will fail, but merged prediction/experimental table will be output.",
            call. = FALSE)
  }
  message("Done")
}

qcPrediction <- function(pred.list, pred.config)  {
  # Ensure that the fill value for each prediction column is at the extreme end of its range
  message("Running QC on predictions")

  doOnePred <- function(pred, config) {
    pred <- as.data.table(pred)
    this.cols <- intersect(colnames(pred), config$pred.col)
    lapply(this.cols, function(s) {
      qcCol(s, 
            pred[, ..s],
            subset(config, pred.col == s)$fill.val, 
            subset(config, pred.col == s)$lowerIsMoreConfident)
      })
  }
  
  qcCol <- function(col.name, colData, fill.val, isInverted) {
    # For each prediction column check that its missing fill val is at the extreme end of its range
    isBad <- (isInverted & fill.val < pmin(colData)) | (!isInverted & fill.val > pmin(colData))
    suppressWarnings(if (isBad) stop(paste0("Fill val for column ", col.name, " is not at the extreme of its range!", fill.val, " ", pmin(colData))))
  }

  dummy <- lapply(pred.list, function(s) doOnePred(s, config = pred.config))
  message("Done")
}

# check for valid pred_config.txt input
qcPredConfig <- function(pred_config) {
  
  # check that there is no baseline prediction set. baseline is used internally
  if ("baseline" %in% pred_config$pred_id) {
    stop("Prediction set called 'baseline' not allowed. Please rename.", call. = FALSE)
  }
  
  # check that pred_id and pred_col create a unique identifier
  config_rows <- nrow(pred_config)
  unique_ids <- nrow(unique(pred_config[, c("pred_id", "pred_col"), with = FALSE]))
  if (config_rows != unique_ids) {
    stop("pred_id and pred_col in pred_config do not provide unique identifiers for predictions.",
         call. = FALSE)
  }
  
}

# check if which genes in experimental also occur in predictions
checkExistenceOfExperimentalGenesInPredictions <- function(expt, pred_list, summary_file) {
  
  # all genes in experimental data
  expt_genes <- sort(unique(expt$GeneSymbol))
  
  # check which genes occur in predictions
  expt_genes_in_pred <- lapply(pred_list, function(pred) {expt_genes %in% unique(pred$GeneSymbol) })
  
  # create summary containing all experimental genes and their occurence in each prediction set
  summary <- data.table(experimental_genes = expt_genes, as.data.table(expt_genes_in_pred))
  
  # write to output file
  write.table(summary, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
}

# combine a list of predictions with experimental data
combineAllExptPred <- function(expt, pred_list, config, cellMapping, outdir) {
  
  # merge each set of predictions with experimental data
  prediction_sets <- structure(names(pred_list), names = names(pred_list))
  merged_list <- lapply(
    prediction_sets,
    function(p) {
      combineSingleExptPred(expt = expt, pred = pred_list[[p]], pred_name = p,  config = config,
                            cellMapping = cellMapping, outdir = outdir)
                          })
  
  # columns to merge each set of merged predictions
  #merge_by_cols <- colnames(expt)
  
  # add class to merging columns if provided in input
  # if ("class" %in% colnames(expt)) merge_by_cols <- c(merge_by_cols, "class")
  
  # merge all data sets in merged_list
  #output <- Reduce(function(x, y) merge(x, y, by = merge_by_cols, all = TRUE), merged_list)
  
  # combine merged data into one table
  output <- rbindlist(merged_list, idcol = "pred_id")
  
  # rearrange columns for output
  expt_cols <- colnames(output)[seq(2, ncol(output) - 3)]
  output_col_order <- c(expt_cols, "pred_id", "pred_col", "pred_value", "Prediction")
  output <- output[, output_col_order, with = FALSE]
  
  return(output)

}

# merge predictions with experimental data
combineSingleExptPred <- function(expt, pred, pred_name, config, cellMapping, outdir) {
  
  # Step 1: merging overlapping enhancer - gene pairs ----------------------------------------------
  
  message("Overlapping predictions for predictions: ", pred_name)
  
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
  
  # map cell types
  if (!is.null(cellMapping)) {
    message("Mapping cell types of predictions and experimental data")
    pred <- applyCellTypeNameMapping(pred, cellMapping)
  }
  
  # create GenomicRanges for CRE-G links for both experimental data and predictions. this applies a
  # trick with using the seqnames to restrict overlaps to E-G pairs involving the same genes and in
  # the same cell type
  expt_gr <- with(expt, GRanges(seqnames = paste0(CellType,":",chrPerturbationTarget,":",GeneSymbol),
                                ranges = IRanges(startPerturbationTarget, endPerturbationTarget)))
  pred_gr <- with(pred, GRanges(seqnames = paste0(CellType,":",chrElement,":",GeneSymbol),
                                ranges = IRanges(startElement, endElement)))
  
  # make sure that seqnames are the same, else GRanges will report unnecessary warnings. this could
  # be removed and replaces with suppressWarnings() when calling findOverlaps() for a small gain
  # in computation time
  seqlevels_all_pairs <- as.character(unique(c(seqnames(expt_gr), seqnames(pred_gr))))
  seqlevels(expt_gr) <- seqlevels_all_pairs
  seqlevels(pred_gr) <- seqlevels_all_pairs
 
  # find overlaps between predictions and experimental data
  ovl <- findOverlaps(expt_gr, pred_gr)
  
  # merge predictions with experimental data
  pred_merge_cols <- config_filt$pred_col  # columns in predictions to add to experimental data
  merged <- cbind(expt[queryHits(ovl)], pred[subjectHits(ovl), pred_merge_cols, with = FALSE])

  # Step 2: aggregating pairs with multiple overlaps -----------------------------------------------
  
  # sometimes perturbed elements will overlap multiple predicted elements (eg in the case of a large
  # deletion). in these cases need to summarize, E.g., sum ABC.Score across model elements
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
  
  # write these to a text file in the output subdirectory
  expt_missing_pred_dir <- file.path(outdir, "experimentalDataMissingPredictions")
  dir.create(expt_missing_pred_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(
    x = expt_missing_preds[, -c("chrTSS", "startTSS", "endTSS")],
    file = file.path(expt_missing_pred_dir, paste0(pred_name, "_expt_missing_predictions.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
    )

  # fill in missing values
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
  output <- melt(merged, measure.vars = pred_merge_cols, variable.name = "pred_col",
                 value.name = "pred_value")
  
  # sort output according to genomic coordinates of enhancers and target gene
  sortcols <- c("chrPerturbationTarget", "startPerturbationTarget", "GeneSymbol")
  setorderv(output, sortcols)
  
  return(output)
  
}

# aggregate experimental enhancers overlapping multiple predicted enhancers
collapseEnhancersOverlappingMultiplePredictions <- function(df, config, agg_cols) {
  
  # summarize columns as defined in config
  list_for_agg <- as.list(df[, ..agg_cols])
  all_list <- mapply(FUN = function(pred_col, agg_func) {
      aggregate(df[, ..pred_col], by = list_for_agg, FUN = agg_func)
    }, pred_col = config$pred_col, agg_fun = config$aggregate_function, SIMPLIFY = FALSE)
  
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

  if ("class" %in% colnames(df)) {
    class_temp <- aggregate(df$class, by = list_for_agg, FUN = class_agg)
    colnames(class_temp)[colnames(class_temp) == "x"] <- "class"
    all_list$class <- class_temp
  }
  
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

applyCellTypeNameMapping <- function(df, cellMapping) {
  #Map CellType in predictions file to match experimental data file
  for (ii in seq(nrow(cellMapping))) {
    this.from <- strsplit(cellMapping$from[ii], split=",")[[1]]
    this.to <- cellMapping$to[ii]
    df$CellType[df$CellType %in% this.from] <- this.to
  }

  return(df)
  
}

# filter experimental data for genes in gene universe and add TSS coordinates to experimental data
filterExptGeneUniverse <- function(expt, genes, missing_file = NULL) {
  
  # add gene TSSs from specified gene universe to experimental data
  expt_cols <- colnames(expt)  # original column order used for output later
  expt <- merge(expt, genes, by = "GeneSymbol", all.x = TRUE, sort = FALSE)
  
  # get experimental data from genes that do not appear in the gene universe
  expt_missing <- expt[is.na(expt$chrTSS), ]
  
  # if there are any, write these to output file (if specified)
  if (nrow(expt_missing) > 0) {
    message("cre-gene pairs in experimental data missing from gene universe: ", nrow(expt_missing))
    if (!is.null(missing_file)) {
      write.table(expt_missing[, -c("chrTSS", "startTSS", "endTSS")], file = missing_file,
                  quote = FALSE, row.names = FALSE, sep = "\t")
    }
  }
  
  # filter out missing expt data, arrange columns and return as output
  expt_filt <- expt[!is.na(expt$chrTSS), c(expt_cols, "chrTSS", "startTSS", "endTSS"), with = FALSE]
  
  return(expt_filt)
  
}

# compute baseline 'distance to TSS' predictor
computeDistToTSS <- function(expt) {
  
  # compute distance to TSS
  dist_to_tss <- with(
    expt,
    abs((startPerturbationTarget + endPerturbationTarget) / 2 - (startTSS + endTSS) / 2)
  )
  
  # create output table
  output <- data.table(
    expt,
    pred_id = "baseline",
    pred_col = "distToTSS",
    pred_value = dist_to_tss,
    Prediction = 1
  )
  
  return(output)
  
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
