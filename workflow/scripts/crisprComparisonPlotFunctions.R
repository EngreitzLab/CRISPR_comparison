## Plotting functions for CRISPR comparison

library(tidyverse)
library(data.table)

## MAIN FUNCTIONS ==================================================================================

# make scatter plots of one column (y) against all predictors (x) (default: combined == all cells)
predScatterPlots <- function(df, predictors, y_col, ncol_facet = 3, point_size = 2, text_size = 13,
                             cell_type = "combined") {
  
  # get data for specified cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # convert to long format for faceting
  df_ct <- df_ct %>% 
    pivot_longer(cols = all_of(predictors), names_to = "pred_name", values_to = "pred_value")
  
  # plot each predictor against effect size
  ggplot(df_ct, aes(x = pred_value, y = get(y_col), color = scatterplot_color)) +
    facet_wrap(~pred_name, scales = "free", ncol = ncol_facet) +
    geom_point(size = point_size) +
    scale_color_manual(values = c("Activating" = "red", "Repressive" = "blue", 
                                  "Not Significant" = "gray")) +
    labs(title = "Predictors vs EffectSize", x = "Predictor value", y = y_col, color = "") +
    scale_x_continuous(labels = scales::scientific) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = text_size))
  
}

# compute PR curves for a given cell type (default: combined == all cells)
calcPRCurves <- function(df, predictors, pos_col, cell_type = "combined") {
  
  # extract data for the provided cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # create list of PR curves with one element per predictor
  pr <- lapply(predictors, function(p) {
    performance(prediction(as.numeric(unlist(merged[, ..p])), unlist(merged[, ..pos_col])), 
                measure = "prec", x.measure = "rec")
  })
  
  return(pr)
  
}

# create a summary for a precision-recall curve table
makePRSummaryTable <- function(pr, min_sensitivity = 0.7, thresholds = numeric(0)) {
  
  # compute AUC
  # the head() calls here remove the last element of the vector. 
  # The point is that performance objects produced by ROCR always include a Recall=100% point even
  # if the predictor cannot achieve a recall of 100%. This results in a straight line ending at
  # (1,0) on the PR curve. This should not be included in the AUC computation.
  auc <- lapply(pr, function(s) {
    computeAUC(head(s@x.values[[1]], -1), head(s@y.values[[1]], -1))
    }) 
  
  # get cutoff for the specified minimum sensitivity
  cutoff <- lapply(pr, function(s) computeCutoffGivenDesiredSensitivity(s, min_sensitivity) )
  max_F1 <- lapply(pr, function(s) {
    max(2 / ((1/s@x.values[[1]]) + (1/s@y.values[[1]])), na.rm = TRUE)
    })
  
  # create performance summary list
  perf_summary_list <- list(cutoff = as.data.table(cutoff), AUC = as.data.table(auc),
                            maxF1 = as.data.table(max_F1))
  
  # convert to table
  perf_summary <- t(rbindlist(perf_summary_list))
  colnames(perf_summary) <- c("cutoff", "AUPRC","maxF1")
  
  # add predictor names and min sensitivity
  perf_summary <- data.table(predictor = rownames(perf_summary), min_sensitivity = min_sensitivity,
                             perf_summary)
  
  # add predictor cutoff if specified in config for given predictor
  perf_summary <- perf_summary %>% 
    mutate(isScore = grepl(predictor, pattern = ".Score$")) %>% 
    mutate(pred_name = sub(predictor, pattern = "\\..*$", replacement = "")) %>% 
    mutate(cutoff = if_else(isScore == TRUE & pred_name %in% names(thresholds),
                            true = thresholds[pred_name], false = cutoff)) %>% 
    select(-c(isScore, pred_name))
  
  ## CHECK WHAT THIS IS DOING
#  for (name in pred.table$name){
#    col.name <- perf_summary$predictor[grepl(paste0(name), perf_summary$predictor) & grepl("Score", perf_summary$predictor)]
#    print(col.name)
#    if (!is.na(pred.table[pred.table$name == name, "threshold"])){
#      perf_summary[perf_summary$predictor == col.name, "cutoff"] <- getThresholdValue(pred.table, name)
#    }
#  }
  
  return(perf_summary)
  
}

# make a PR curve plot for a set of provided predictors
makePRCurvePlot <- function(pr_df, pred_cols, pct_pos, min_sensitivity = 0.7,
                            plot_name = "PRC full experimental data", point_size = 3,
                            text_size = 15) {
  
  # only retain pred_cols specified via pred_cols
  pr_df <- subset(pr_df, pred_col %in% pred_cols)
  
  # separate boolean predictors from continuous predictors. TODO: move to predConfig file??
  pr_cutoff <- by(pr_df, pr_df$pred_col, function(df) unique(df$alpha))
  boolean_predictors <- names(pr_cutoff)[unlist(lapply(pr_cutoff, function(s) identical(s, c(Inf, 1, 0))), use.names = FALSE)]
  cont_pred <- subset(pr_df, !(pred_col %in% boolean_predictors))
  bool_pred <- subset(pr_df, pred_col %in% boolean_predictors)
  bool_pred <- subset(bool_pred, alpha == 1)
  
  # get thresholded values
  pred_cols <- unique(pr_df$pred_col)
  threshold_pred <- getPrecisionAtRecall(pr_df, col_list = pred_cols,
                                         min_sensitivity = min_sensitivity)
  
  # convert pred_cols to factors with same levels to ensure proper plotting
  cont_pred$pred_col <- factor(cont_pred$pred_col, levels = pred_cols)
  bool_pred$pred_col <- factor(bool_pred$pred_col, levels = pred_cols)
  threshold_pred$pred_col <- factor(threshold_pred$pred_col, levels = pred_cols)
  
  # create plot starting from continuous predictors then add boolean predictors and fluff
  ggplot(cont_pred, aes(x = recall, y = precision, color = pred_col)) + 
    geom_line(size = 1) +
    geom_point(data = bool_pred, size = point_size) +
    geom_point(data = threshold_pred, size = point_size) +
    geom_hline(yintercept = pct_pos, linetype = "dashed", color = "black") +
    labs(title = plot_name, x  = "Recall", y = "Precision", color = "") + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    theme_bw() +
    theme(text = element_text(size = text_size))
  
}

## HELPER FUNCTIONS ================================================================================

# add label for each pair based on whether it's significant and activates or represses it's target
labelPairs <- function(df, sig_col = "Significant") {
  
  df$scatterplot_color <- with(df,
                               ifelse(get(sig_col) == FALSE,
                                      yes = "Not Significant",
                                      no = ifelse(EffectSize > 0,
                                                  yes = "Repressive",
                                                  no  = "Activating")
                                      )
                               )
  
  return(df)
  
}

# get data for a specified cell type or all if cell_type == "combined"
getCellTypeData <- function(df, cell_type) {
  
  if (cell_type != "combined") {
    df <- subset(df, CellType == cell_type)
  } 
  
  return(df)
}

# convert a list of ROCR performance objects into a table and calculate F1 metric
pr2df <- function(pr, calc_f1 = TRUE) {

  # function to convert one performance object to a table
  convert_pr2df <- function(this_pr) {
    df <- as.data.frame(list(
      alpha = this_pr@alpha.values[[1]],
      precision = this_pr@y.values[[1]],
      recall = this_pr@x.values[[1]]
    ))
    return(df)
  }
  
  # apply to input list
  pr_list <- lapply(pr, convert_pr2df)
  
  # convert list of tables into one table
  pr_df <- rbindlist(pr_list, idcol = "pred_col")
  
  # split predictor into prediction dataset name and predictor name within dataset
  pr_df$pred_name <- sub(pr_df$pred_col, pattern = "^([[:alnum:]]+)\\..+$", replacement = "\\1")
  pr_df$predictor <- sub(pr_df$pred_col, pattern = "^[[:alnum:]]+\\.(.+)$", replacement = "\\1")
  
  # calculate F1 metric if specified
  if (calc_f1 == TRUE) {
    pr_df$F1 <- with(pr_df, 2 / ((1/precision) + (1/recall)))
  }
           
  return(pr_df)
  
}

# compute AUC
computeAUC <- function(x_vals, y_vals) {
  good.idx <- which(!is.na(x_vals) & !is.na(y_vals))
  return(trapz(x_vals[good.idx], y_vals[good.idx]))
}

# compute a predictor cutoff for provided sensitivity
computeCutoffGivenDesiredSensitivity <- function(pr, min_sensitivity) {
  sens <- pr@x.values[[1]]
  prec <- pr@y.values[[1]]
  cutoff.sensitivity <- min(sens[sens >= min_sensitivity])
  idx <- which.max(sens == cutoff.sensitivity)
  idx2 <- idx[which.max(prec[idx])]
  
  cutoff <- pr@alpha.values[[1]][idx2]
  return(cutoff)
}

# assign prediction class label
addPredictionClassLabels <- function(df, perf_summary, pos_col = "Regulated") {
  
  for (i in seq_len(nrow(perf_summary))) {
    df <- addOneLabel(df, cutoff = as.numeric(perf_summary$cutoff[[i]]),
                      score_col = perf_summary$predictor[[i]], pos_col = pos_col)
  }
  
  return(df)
}

# assign prediction class labels for one predictor
addOneLabel <- function(df, cutoff, score_col, pos_col) {
  
  label_name <- paste0(score_col, ".pred.class")
  df[, label_name] <- "NA"
  
  df[which(!is.na(df[, ..score_col]) & df[, ..score_col] > cutoff & df[, ..pos_col]), label_name] <- "TP"
  df[which(!is.na(df[, ..score_col]) & df[, ..score_col] <= cutoff & !df[, ..pos_col]), label_name] <- "TN"
  df[which(!is.na(df[, ..score_col]) & df[, ..score_col] > cutoff & !df[, ..pos_col]), label_name] <- "FP"
  df[which(!is.na(df[, ..score_col]) & df[, ..score_col] <= cutoff & df[, ..pos_col]), label_name] <- "FN"
  
  return(df)
}

# compute the number of positives for a given cell type (default: combined == all cells)
calcPctPos <- function(df, pos_col, cell_type = "combined") {
  
  # extract data for the provided cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # compute percentage of positives
  mean(unlist(df_ct[, ..pos_col]))

}

# get precision at a specific sensitivity (recall)
getPrecisionAtRecall <- function(pr_df, col_list, min_sensitivity) {
  
  thresholds <- data.frame(matrix(ncol=3, nrow=length(col_list)))
  colnames(thresholds) <- c("precision", "recall", "pred_col")
  thresholds$pred_col <- col_list
  thresholds$recall <- min_sensitivity
  
  prec <- lapply(col_list, function(s) {
    recall <- as.numeric(unlist(pr_df[pr_df$pred_col==s, "recall"]))
    precision <- as.numeric(unlist(pr_df[pr_df$pred_col==s, "precision"]))
    approxfun(x = recall, y = precision)(min_sensitivity)
  })
  
  thresholds$precision = as.numeric(unlist(prec))
  return(thresholds)
  
}
