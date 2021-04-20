library(GenomicRanges)

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

checkExistenceOfExperimentalGenesInPredictions <- function(expt, pred.list, outdir) {
  experimentalGenes <- unique(expt$GeneSymbol)
  res <- sapply(pred.list, function(df) {experimentalGenes %in% unique(df$GeneSymbol)})
  df <- cbind(experimentalGenes, as.data.table(res))
  write.table(df, file.path(outdir, "ExperimentalGenesAppearingInPredictions.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
}

# combine a list of predictions with experimental data
combineAllExptPred <- function(expt, pred_list, config, cellMapping, outdir, fill_missing = TRUE) {
  
  # merge each set of predictions with experimental data
  merged_list <- lapply(
    names(pred_list),
    function(s) {
      combineSingleExptPred(expt = expt, pred = pred_list[[s]], #threshold = threshold[[s]],
                            pred_name = s, config = config, cellMapping = cellMapping, 
                            outdir = outdir, fill_missing = fill_missing)
                          })
  
  # columns to merge each set of merged predictions
  merge_by_cols <- c("chrPerturbationTarget", "startPerturbationTarget", "endPerturbationTarget",
                     "GeneSymbol", "CellType", "Significant", "Regulated", "EffectSize",
                     "IncludeInModel", "chrTSS", "startTSS", "endTSS")
  
  # add class to merging columns if provided in input
  if ("class" %in% colnames(expt)) merge_by_cols <- c(merge_by_cols, "class")
  
  # merge all data sets in merged_list
  output <- Reduce(function(x, y) merge(x, y, by = merge_by_cols, all = TRUE), merged_list)
  return(output)

}

# merge predictions with experimental data
combineSingleExptPred <- function(expt, pred, pred_name, config, cellMapping, outdir,
                                  fill_missing) {
  
  # Step 1: merging overlapping enhancer - gene pairs ----------------------------------------------
  
  message("Overlapping predictions for predictor: ", pred_name)
  
  # subset config to columns that appear in predictions, otherwise code will fail
  config_filt <- subset(config, pred.col %in% colnames(pred))
  missing_pred <- setdiff(config$pred.col, config_filt$pred.col)
  if (length(missing_pred) > 0) {
    warning("Following predictor(s) specified in config file not found for ", pred_name, ": ",
            paste(missing_pred, collapse = ", "), call. = FALSE)
  }
  
  # map cell types
  if (cellMapping != "") {
    message("Mapping cell types of predictions and experimental data")
    pred <- applyCellTypeNameMapping(pred, cellMapping)
  }
  
  # create GenomicRanges for both predictions and experimental data
  pred_gr <- with(pred, GRanges(seqnames = paste0(CellType,":",chrElement,":",GeneSymbol),
                                ranges = IRanges(startElement, endElement)))
  expt_gr <- with(expt, GRanges(seqnames = paste0(CellType,":",chrPerturbationTarget,":",GeneSymbol),
                                ranges = IRanges(startPerturbationTarget, endPerturbationTarget)))
  
  # find overlaps between predictions and experimental data
  ovl <- GenomicRanges::findOverlaps(expt_gr, pred_gr)
  
  # merge predictions with experimental data
  pred_merge_cols <- c("chrTSS", "startTSS", "endTSS", config_filt$pred.col)
  merged <- cbind(expt[queryHits(ovl)], pred[subjectHits(ovl), pred_merge_cols, with = FALSE])

  # Step 2: aggregating pairs with multiple overlaps -----------------------------------------------
  
  # sometimes perturbed elements will overlap multiple predicted elements (eg in the case of a large
  # deletion). in these cases need to summarize, Eg sum ABC.Score across model elements overlapping
  # the deletion. this requires a config file describing how each prediction column should be
  # aggregated
  agg_cols <- setdiff(colnames(merged), config_filt$pred.col)
  merged <- collapseEnhancersOverlappingMultiplePredictions(merged, config_filt, agg_cols)
  
  # Step 3: experimental data missing predictions --------------------------------------------------
  
  # A tested enhancer element may not have a prediction. For ABC this is typically the case if the
  # tested element does not overlap a DHS peak. In this case we need to fill the predictions table
  
  # add 'Prediction' column to merged with value 1, indicating pairs that overlap with predictions
  merged$Prediction <- 1
  
  # get pairs from experimental data that are missing from predictions 
  expt_missing_predictions <- expt[setdiff(seq_len(nrow(expt)), queryHits(ovl)), ]
  
  # write these to a text file in the output directory
  dir.create(file.path(outdir, "experimentalDataMissingPredictions", pred_name), recursive = TRUE,
             showWarnings = FALSE)
  write.table(expt_missing_predictions,
              file = file.path(outdir, "experimentalDataMissingPredictions", pred_name,
                               "expt_missing_predictions.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # fill in missing predictions if specified
  if (fill_missing == TRUE) {

    # fill in missing values
    expt_missing_predictions <- fillMissingPredictions(expt_missing_predictions,
                                                       config = config_filt, agg_cols)
    
    # add 'Prediction' column with value 0, indicating pairs that were not found in predictions
    expt_missing_predictions$Prediction <- 0
    
    # merge filled data with merged data
    cols_we_want <- colnames(merged)
    merged <- rbind(merged, expt_missing_predictions[, ..cols_we_want])
    
    message("Experimental data missing predictions filled. Will be in merged output!")
    
  } else {
    
    message("Experimental data missing predictions ignored. Will not be in merged output!")
    
  }
  
  # rename new columns based on prediction dataset name
  pred_cols <- colnames(merged) %in% c(config_filt$pred.col, "Prediction")
  colnames(merged)[pred_cols] <- paste0(pred_name, ".", colnames(merged)[pred_cols])
  
  ### WHY IS THIS NECESSARY??
  #score_colname <- colnames(merged)[grep("Score", colnames(merged))]
  #merged[[paste0(score_colname,'.Prediction')]] <- merged[[score_colname]] * merged$Prediction
  
  # get distance for predictions that overlap experiments
  #dist_colname <- colnames(merged)[grep("distance", colnames(merged))]
  #merged[[paste0(dist_colname,'.Prediction')]] <- merged[[dist_colname]] * merged$Prediction
  
  return(merged)
}

# aggregate experimental enhancers overlapping multiple predicted enhancers
collapseEnhancersOverlappingMultiplePredictions <- function(df, config, agg_cols) {
  
  # summarize columns as defined in config
  list_for_agg <- as.list(df[, ..agg_cols])
  all_list <- mapply(FUN = function(pred_col, agg_func) {
      aggregate(df[, ..pred_col], by = list_for_agg, FUN = agg_func)
    }, pred_col = config$pred.col, agg_fun = config$agg.func, SIMPLIFY = FALSE)
  
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
  output <- Reduce(function(df1, df2) merge(df1, df2, by = agg_cols), all_list)
  return(output)
  
}

# fill in prediction values for experimental pairs missing from prediction data
fillMissingPredictions <- function(df, config, agg_cols) {
  
  # fill in missing predictions as described in the config file
  for (i in seq_len(nrow(config))) {
    df[, config$pred.col[[i]]] <- config$fill.val[i]
  }
  
  # fill in unknown columns (columns appearing in predictions data, but not experiment or agg_cols)
  unk_cols <- setdiff(c("class", agg_cols), unique(c(colnames(df), config$pred.cols)))
  df[, unk_cols] <- "Merge:UNKNOWN"
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

writeExptSummary <- function(df, outdir) {
  df.summary <- as.data.frame(list(
    numConnections = nrow(df),
    numIncludeInModel = sum(df$IncludeInModel),
    numIncludeInModelRegulated = sum(df$IncludeInModel & df$Regulated)
  ))
  
  write.table(df.summary, file.path(outdir, "expt.summary.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
}

fread_ignore_empty <- function(f, ...) {
  tryCatch({
    return(fread(f, fill = TRUE, ...))
  }, error = function(e){
    print(paste0("Could not open file: ", f))
    return()
  })
}

fread_gz_ignore_empty <- function(f, ...) {
  tryCatch({
    return(fread(paste0("gunzip -c ", f), ...))
  }, error = function(e){
    print(paste0("Could not open file: ", f))
    return()
  })
}

smart_fread <- function(f, ...) {
  if (summary(file(f))$class == "gzfile") {
    out <- fread_gz_ignore_empty(f, ...)
  } else {
    out <- fread_ignore_empty(f, ...)
  }
  
  #con <- file(f)
  #on.exit(close(con), add = TRUE)
  tryCatch({
    closeAllConnections()
  }, error = function(e) {
    print(e)
  }
  )
  
  return(out)
}

loadDelimited <- function(file.list) {
  data.list <- lapply(file.list, smart_fread)
  return(rbindlist(data.list, fill = TRUE))
}

loadFileString <- function(file.str, delim = ",") {
  file.list <- strsplit(file.str, split = delim)[[1]]
  return(loadDelimited(file.list))
}

loadPredictions <- function(pred.table) {
  pred.list <- lapply(pred.table$path, function(s) {
    message("Loading dataset: ", s)
    df <- loadFileString(s)
    message("\tDataset loaded with ", nrow(df), " rows")
    return(df)
    })
  names(pred.list) <- pred.table$name
  return(pred.list)
}
