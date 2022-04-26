## Plotting functions for CRISPR comparison

library(tidyverse)
library(data.table)

## MAIN FUNCTIONS ==================================================================================

# plot CRISPR E-G pairs overlapping prediction E-G pairs
plotOverlaps <- function(df, cell_type = "combined"){
  
  # extract data for the provided cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # count number of total CRISPR E-G pairs and overlapping pairs per predictor
  n_pairs <- df_ct %>% 
    group_by(pred_uid, pred_name_long) %>% 
    summarize(`Overlaps predictions` = sum(Prediction == 1),
              `Not in predictions` = sum(Prediction == 0),
              .groups = "drop") %>% 
    pivot_longer(cols = c(`Overlaps predictions`, `Not in predictions`), names_to = "Overlaps",
                 values_to = "pairs")
  
  # plot number of CRISPR E-G pairs overlapping E-G pairs in predictions
  ggplot(n_pairs, aes(x = pred_name_long, y = pairs, fill = Overlaps)) +
    geom_bar(stat = "identity") +
    labs(y = "E-G pairs", x = "Predictor", title = "CRISPR E-G pairs overlapping predictions") +
    scale_fill_manual(values = c("Overlaps predictions" = "steelblue",
                                 "Not in predictions" = "darkgray")) +
    coord_flip() +
    theme_bw()
  
}

# compute PR curves for a given cell type (default: combined == all cells)
calcPRCurves <- function(df, pred_config, pos_col, cell_type = "combined") {
  
  # extract data for the provided cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)

  # split into list for lapply
  df_ct_split <- split(df_ct, f = df_ct$pred_uid)
  
  # get inverse predictors
  inverse_predictors <- pred_config %>% 
    select(pred_uid, inverse_predictor) %>% 
    deframe()
  
  # multiply inverse predictors by -1 so that higher value corresponds to higher score
  inverse_predictors <- inverse_predictors[names(df_ct_split)]  # same as predictors for cell type
  df_ct_split <- mapply(FUN = function(pred, inv_pred) {
    inv_multiplier <- ifelse(inv_pred, -1, 1)
    pred$pred_value <- pred$pred_value * inv_multiplier
    return(pred)
  }, df_ct_split, inverse_predictors, SIMPLIFY = FALSE)
  
  # compute precision-recall performance for each predictor
  pr <- lapply(df_ct_split, FUN = function(p){
    performance(prediction(p$pred_value, p[[pos_col]]), measure = "prec", x.measure = "rec")
  })
  
  # convert to table and calculate F1
  pr_df <- pr2df(pr, calc_f1 = TRUE)
  
  return(pr_df)
  
}

# create performance summary table for all predictors in a PR table
makePRSummaryTable <- function(pr_df, pred_config, min_sensitivity = 0.7) {
  
  # remove any boolean predictors since the following metrics don't make sense for them
  bool_preds <- pull(filter(pred_config, boolean == TRUE), pred_uid)
  pr_df <- filter(pr_df, !pred_uid %in% bool_preds)
  
  # compute performance summary for each predictor
  perf_summary <- pr_df %>% 
    group_split(pred_uid) %>% 
    lapply(calcPerfSummaryOnePred, pred_config = pred_config, min_sensitivity = min_sensitivity) %>% 
    bind_rows()
    
  # add predictor information from pred_config
  perf_summary <- pred_config %>% 
    select(pred_uid, pred_id, pred_col, inverse_predictor, pred_name_long) %>% 
    left_join(x = perf_summary, y = ., by = "pred_uid")
  
  return(perf_summary)
  
}

# make a PR curve plot for a set of provided predictors
makePRCurvePlot <- function(pr_df, pred_config, pct_pos, min_sensitivity = 0.7,
                            plot_name = "PRC full experimental data", line_width = 1, 
                            point_size = 3, text_size = 15, colors = NULL) {
  
  # create performance summary
  perf_summary <- makePRSummaryTable(pr_df, pred_config = pred_config,
                                     min_sensitivity = min_sensitivity)
  
  # get PR values at threshold
  pr_threshold <- perf_summary %>% 
    select(pred_name_long, precision = precision_at_cutoff, recall = sensitivity_at_cutoff)
  
  # add pretty predictor names to pr_df for plotting
  pr_df <- left_join(pr_df, select(pred_config, pred_uid, pred_name_long), by = "pred_uid")
  
  # separate pr data into quantitative and boolean predictors
  bool_preds <- pull(filter(pred_config, boolean == TRUE), pred_uid)
  pr_quant <- filter(pr_df, !pred_uid %in% bool_preds)
  pr_bool <- filter(pr_df, pred_uid %in% bool_preds)
  
  # get precision and recall for boolean predictor at alpha 1
  pr_bool <- filter(pr_bool, alpha == 1)
  
  # create PRC plot (caution, this assumes that there at least 1 quant and 1 bool predictor!)
  p <- ggplot(pr_quant, aes(x = recall, y = precision, color = pred_name_long)) +
    geom_line(size = line_width) +
    geom_point(data = pr_threshold, size = point_size) +
    geom_point(data = pr_bool, size = point_size) +
    geom_hline(yintercept = pct_pos, linetype = "dashed", color = "black") +
    labs(title = plot_name, x  = "Recall", y = "Precision", color = "Predictor") + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    theme_bw() +
    theme(text = element_text(size = text_size))
  
  # add custom colors if any are provided
  if (!is.null(colors)) {
    p <- p + 
      scale_color_manual(values = colors)
  }
  
  # print plot
  p
  
}

# make scatter plots of one column (y) against all predictors (x) (default: combined == all cells)
predScatterPlots <- function(df, y_col, pred_names_col = "pred_uid", point_size = 2,
                             text_size = 13, alpha_value = 1, cell_type = "combined", ncol = NULL,
                             nrow = NULL) {
  
  # get data for specified cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # plot each predictor against effect size
  ggplot(df_ct, aes(x = pred_value, y = get(y_col), color = scatterplot_color)) +
    facet_wrap(~get(pred_names_col), scales = "free", ncol = ncol, nrow = nrow) +
    geom_point(size = point_size, alpha = alpha_value) +
    scale_color_manual(values = c("Activating" = "red", "Repressive" = "blue", 
                                  "Not Significant" = "gray")) +
    labs(title = paste("Predictors vs", y_col), x = "Predictor value", y = y_col, color = "") +
    scale_x_continuous(labels = scales::scientific) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = text_size))
  
}

# make violin plots showing scores for each predictor as a function of experimental outcome
plotPredictorsVsExperiment <- function(df, pos_col = "Regulated", pred_names_col = "pred_uid",
                                       text_size = 13, cell_type = "combined", ncol = NULL,
                                       nrow = NULL) {
  
  # get data for specified cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # plot scores for each predictor as a function of experimental outcome
  ggplot(df_ct, aes(x = get(pos_col), y = pred_value, color = get(pos_col), fill = get(pos_col))) +
    facet_wrap(~get(pred_names_col), scales = "free", ncol = ncol, nrow = nrow) +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "NA") +
    scale_color_manual(values = c("darkgray", "steelblue")) +
    scale_fill_manual(values = c("darkgray", "steelblue")) +
    labs(title = "Predictors vs experimental outcome", x = "Experimental E-G pair",
         y = "Predictor value (sqrt)") +
    scale_y_sqrt() +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = text_size))
  
}

# plot distance distributions for experimental positives and negatives
plotDistanceDistribution <- function(df, dist = "baseline.distToTSS", pos_col = "Regulated",
                                     cell_type = "combined", text_size = 13) {
  
  # get data for specified cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # get data for distance and add label for faceting
  dist_data <- df_ct %>% 
    filter(pred_uid == dist) %>% 
    mutate(label = if_else(get(pos_col) == TRUE, true = "Positives", false = "Negatives")) %>% 
    mutate(label = factor(label, levels = c("Positives", "Negatives"))) %>% 
    mutate(dist_kb = pred_value / 1000)
  
  # plot distance distribution for all pairs in experimental data
  ggplot(dist_data, aes(x = dist_kb, fill = label)) +
    facet_wrap(~label, ncol = 1, scales = "free_y") +
    geom_histogram(binwidth = 10) +
    labs(x = "Distance (kb)") +
    scale_fill_manual(values = c("Positives" = "steelblue", "Negatives" = "darkgray")) +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = text_size))
  
}

# make an upset plot of overlapping features for a given cell type
plotOverlappingFeatures <- function(df, feature_cols, cell_type = "combined") {
  
  # get data for specified cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # create table with unique enhancers and overlapping features
  enh_features <- df_ct %>% 
    mutate(enh_id = paste0(chrom, ":", chromStart, "-", chromEnd)) %>% 
    select(enh_id, all_of(feature_cols)) %>% 
    distinct() %>% 
    mutate(across(all_of(feature_cols), ~ .x * 1))
  
  # set new column names
  colnames(enh_features) <- sub("enh_feature_", "", colnames(enh_features))
  
  # create upset plot with overlapping features
  upset(enh_features, nsets = 10, order.by = "freq", number.angles = 30, point.size = 3.5,
        line.size = 2, mainbar.y.label = "Enhancers overlapping features",
        sets.x.label = "Overlapping sites",
        text.scale = c(1.5, 1.5, 1.2, 1.2, 1.5, 1.5))
  
}

## Make plots for subsets of the data based on gene or enhancer features ---------------------------

# calculate and plot PR curves for a given subset
makePRCurveSubset <- function(df, subset_col, pred_config, pos_col, min_sensitivity = 0.7,
                              line_width = 1, point_size = 3, text_size = 15, nrow = 1,
                              colors = NULL) {
  
  # split df into subsets based on provided column
  df_split <- split(df, f = df[[subset_col]])
  
  # compute PR curve
  prc <- lapply(df_split, FUN = calcPRCurves, pred_config = pred_config, pos_col = pos_col)
  
  # calculate percentage of positives for all splits
  pct_pos <- lapply(df_split, FUN = calcPctPos, pos_col = pos_col)
  
  # create plot titles based on feature name and number of pairs for that feature
  feature_name <- sub(".+_feature_", "", subset_col)
  pairs <- vapply(df_split, FUN = function(x) length(unique(x$name)), FUN.VALUE = integer(1))
  titles <- paste0(feature_name, " = " , names(prc), " (", pairs, " pairs)")
  
  # plot PR curves
  pr_plots <- mapply(FUN = makePRCurvePlot, pr_df = prc, pct_pos = pct_pos, plot_name = titles, 
                     MoreArgs = list(pred_config = pred_config, min_sensitivity = min_sensitivity,
                                     line_width = line_width, point_size = point_size,
                                     text_size = text_size, colors = colors),
                     SIMPLIFY = FALSE)
  
  # create title for this comparison
  title <- ggdraw() + 
    draw_label(feature_name, fontface = "bold", x = 0, hjust = 0, size = text_size * 1.5) +
    theme(plot.margin = margin(0, 0, 0, 7)) # align title with left edge of first plot
  
  # combine into one plot
  pr_plots_row <- plot_grid(plotlist = pr_plots, nrow = nrow)
  plot_grid(title, pr_plots_row, ncol = 1, rel_heights = c(0.1, 1))
  
}

# calculate and plot PR curves for several subset columns and arrange plots into one figure for one
# cell type
makePRCurveSubsets <- function(df, subset_cols, pred_config, pos_col, cell_type = "combined",
                               min_sensitivity = 0.7, line_width = 1, point_size = 3,
                               text_size = 15, nrow = 1, colors = NULL) {
  
  # return NULL if no subsets are available
  if (length(subset_cols) == 0) {
    return(NULL)
  }
  
  # extract data for the provided cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # create PR curve plots for each subset
  pr_plots <- lapply(subset_cols, FUN = makePRCurveSubset, df = df_ct,
                     pred_config = pred_config, pos_col = pos_col,
                     min_sensitivity = min_sensitivity, line_width = line_width,
                     point_size = point_size, text_size = text_size, nrow = nrow, colors = colors)
  
  # create one figure with all plots
  plot_grid(plotlist = pr_plots, nrow = length(subset_cols))
  
}

# make violin plots showing scores for each predictor vs experimental outcome for a given subset
plotPredVsExperimentSubset <- function(df, subset_col, pos_col, pred_names_col = "pred_uid",
                                      text_size = 13) {
  
  # get feature name for title
  feature_name <- sub(".+_feature_", "", subset_col)
  
  # plot scores for each predictor as a function of experimental outcome
  ggplot(df, aes(x = get(pos_col), y = pred_value, color = get(pos_col), fill = get(pos_col))) +
    facet_grid(get(pred_names_col) ~ get(subset_col), scales = "free") +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "NA") +
    scale_color_manual(values = c("darkgray", "steelblue")) +
    scale_fill_manual(values = c("darkgray", "steelblue")) +
    labs(title = feature_name, x = "Experimental E-G pair",
         y = "Predictor value (sqrt)") +
    scale_y_sqrt() +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = text_size))
  
}

# make violin plots showing scores for each predictor vs experimental outcome for a given subset
plotPredVsExperimentSubsets <- function(df, subset_cols, pos_col, cell_type = "combined",
                                        pred_names_col = "pred_uid", text_size = 13) {
  
  # return NULL if no subsets are available
  if (length(subset_cols) == 0) {
    return(NULL)
  }
  
  # extract data for the provided cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # create predictor vs experiment plots for each subset
  pe_plots <- lapply(subset_cols, FUN = plotPredVsExperimentSubset, df = df_ct, pos_col = pos_col,
                     pred_names_col = pred_names_col, text_size = text_size)
  
  # create one figure with all plots
  plot_grid(plotlist = pe_plots, nrow = length(subset_cols))
  
}

# make violin plots showing scores for each predictor vs experimental outcome for a given subset
predScatterPlotsSubset <- function(df, subset_col, y_col, pred_names_col = "pred_uid",
                                   point_size = 2, text_size = 13, alpha_value = 1) {
  
  # get feature name for title
  feature_name <- sub(".+_feature_", "", subset_col)
  
  # plot each predictor against effect size for the given subset
  ggplot(df, aes(x = pred_value, y = get(y_col), color = scatterplot_color)) +
    facet_grid(get(subset_col) ~ get(pred_names_col), scales = "free") +
    geom_point(size = point_size, alpha = alpha_value) +
    scale_color_manual(values = c("Activating" = "red", "Repressive" = "blue", 
                                  "Not Significant" = "gray")) +
    labs(title = feature_name, x = "Predictor value", y = y_col, color = "") +
    scale_x_continuous(labels = scales::scientific) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = text_size))
  
}

# make violin plots showing scores for each predictor vs experimental outcome for a given subset
predScatterPlotsSubsets <- function(df, subset_cols, y_col, cell_type = "combined",
                                    pred_names_col = "pred_uid", point_size = 2, text_size = 13,
                                    alpha_value = 1) {
  
  # return NULL if no subsets are available
  if (length(subset_cols) == 0) {
    return(NULL)
  }
  
  # extract data for the provided cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # create predictor vs experiment scatter plots for each subset
  sc_plots <- lapply(subset_cols, FUN = predScatterPlotsSubset, df = df_ct, y_col = y_col,
                     pred_names_col = pred_names_col, point_size = point_size,
                     text_size = text_size, alpha_value = alpha_value)
  
  # create one figure with all plots
  plot_grid(plotlist = sc_plots, nrow = length(subset_cols))
  
}
 
## HELPER FUNCTIONS ================================================================================

# get data for a specified cell type or all if cell_type == "combined"
getCellTypeData <- function(df, cell_type) {
  
  if (cell_type != "combined") {
    df <- subset(df, ExperimentCellType == cell_type)
  } 
  
  return(df)
}

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
  pr_df <- rbindlist(pr_list, idcol = "pred_uid")
  
  # calculate F1 metric if specified
  if (calc_f1 == TRUE) {
    pr_df$F1 <- with(pr_df, 2 / ((1 / precision) + (1 / recall)))
  }
           
  return(pr_df)
  
}

# create performace summary for one predictor
calcPerfSummaryOnePred <- function(pr_df, pred_config, min_sensitivity) {
  
  # check that pr_df contains data on only one predictor
  predictor <- unique(pr_df$pred_uid)
  if (length(predictor) > 1) {
    stop("Input pr_df contains data for more than one unique predictor.", call. = FALSE)
  }
  
  # make sure that input is sorted according to recall and precision
  pr_df <- arrange(pr_df, recall, desc(precision))
  
  # compute AUC and maximum F1
  # the head() calls here remove the last element of the vector. 
  # The point is that performance objects produced by ROCR always include a Recall = 100% point even
  # if the predictor cannot achieve a recall of 100%. This results in a straight line ending at
  # (1,0) on the PR curve. This should not be included in the AUC computation.
  auprc <- pr_df %>% 
    head(-1) %>% 
    summarize(AUPRC = computeAUC(x_vals = recall, y_vals = precision),
              max_F1 = computeMaxF1(.))
  
  # compute performance at min sensitivity
  perf_min_sens <- computePerfGivenSensitivity(pr_df, min_sensitivity = min_sensitivity)

  # get cutoff specified in pred_config for given predictor
  cutoff <- pred_config %>% 
    filter(pred_uid == predictor) %>% 
    mutate(alpha = if_else(inverse_predictor == TRUE, true = alpha * -1, false = alpha)) %>% 
    pull(alpha)
  
  # compute performance at specified alpha cutoff
  perf_alpha_cutoff <- computePerfGivenCutoff(pr_df, alpha_cutoff = cutoff)
  
  # create output table
  perf_summary <- data.frame(pred_uid = predictor, auprc, perf_min_sens, perf_alpha_cutoff)
  
  # if AUC couldn't be computed, set all performance metrics to NA, since this indicates a PRC
  # without any real points except the starting and end points added by the ROCR package
  # TODO: find better solution for this
  if (all(is.na(auprc))) {
    perf_cols <- !colnames(perf_summary) %in% c("pred_uid", "min_sensitivity", "alpha_cutoff")
    perf_summary[, perf_cols] <- NA_real_ 
  }
  
  return(perf_summary)
  
}

# try to compute AUC
computeAUC <- function(x_vals, y_vals) {
  good.idx <- which(!is.na(x_vals) & !is.na(y_vals))
  if (length(good.idx) > 0) {
    auc <- trapz(x_vals[good.idx], y_vals[good.idx])
  } else {
    auc <- NA_real_
  }
  return(auc)
}

# try to compute maximum F1 from a prc table
computeMaxF1 <- function(pr_df) {
  if (any(!is.na(pr_df$F1))) {
    maxF1 <- max(pr_df$F1, na.rm = TRUE)
  } else {
    maxF1 <- NA_real_
  }
  return(maxF1)
}

# compute performance given a minimum sensitivity (recall)
computePerfGivenSensitivity <- function(pr_df, min_sensitivity) {
  
  # get required values from pr_df
  prec <- pr_df$precision
  sens <- pr_df$recall
  alpha <- pr_df$alpha
  
  # get sensitivity value that is closest (at least) to minimum sensitivity
  cutoff_sens <- min(sens[sens >= min_sensitivity])
  
  # get indices of all sensitivities equal to that cutoff (can be more than 1)
  idx <- which(sens == cutoff_sens)
  
  # pick the one with the maximum precision
  idx2 <- idx[which.max(prec[idx])]
  
  # get cutoff (alpha), precision and senitivity for that point on the PR curve
  alpha_at_sensitivity <- alpha[idx2]
  prec_at_sensitivity <- prec[idx2]
  sens_at_sensitivity <- sens[idx2]
  
  # create output data.frame row
  output <- data.frame(min_sensitivity = min_sensitivity,
                       alpha_at_min_sensitivity = alpha_at_sensitivity,
                       precision_at_min_sensitivity = prec_at_sensitivity, 
                       sensitivity_at_min_sensitivity = sens_at_sensitivity,
                       stringsAsFactors = FALSE)
  
  return(output)
  
}

# compute a performance given a specified alpha cutoff
computePerfGivenCutoff <- function(pr_df, alpha_cutoff) {
  
  # get sensitivity (recall), precision and alpha values
  sens <- pr_df$recall
  prec <- pr_df$precision
  alpha <- pr_df$alpha
  
  # set alpha cutoff to min alpha if it was NA
  if (is.na(alpha_cutoff)) alpha_cutoff <- min(alpha)
  
  # get index of highest alpha value that is larger or equal to alpha_cutoff
  idx <- sum(alpha >= alpha_cutoff)
  
  # get sensitivity and precision at that cutoff
  sens_at_cutoff <- sens[idx]
  prec_at_cutoff <- prec[idx]
  alpha_at_cutoff <- alpha[idx]
  
  # create output data.frame row
  output <- data.frame(alpha_cutoff = alpha_cutoff,
                       alpha_at_cutoff = alpha_at_cutoff,
                       sensitivity_at_cutoff = sens_at_cutoff,
                       precision_at_cutoff = prec_at_cutoff,
                       stringsAsFactors = FALSE)
  
  return(output)
  
}

# assign prediction class label
addPredictionClassLabels <- function(df, perf_summary, pos_col = "Regulated") {
  
  # add alpha thresholds to df
  df <- perf_summary %>% 
    select(cell_type, pred_uid, alpha_cutoff, inverse_predictor) %>% 
    left_join(x = df, y = ., by = c("ExperimentCellType" = "cell_type", "pred_uid"))
  
  # invert iverse predictor values and generate labels (TP, FP, TN, FN)
  pos_col <- sym(pos_col)  ## create symbol from column name for tidy evaluation
  df <- df %>% 
    mutate(pred_value_inv = if_else(inverse_predictor == TRUE, true = pred_value * -1,
                                    false = pred_value)) %>% 
    mutate(pred_class = case_when(
      pred_value_inv > alpha_cutoff & !!pos_col == TRUE ~ "TP",
      pred_value_inv > alpha_cutoff & !!pos_col == FALSE ~ "FP",
      pred_value_inv <= alpha_cutoff & !!pos_col == TRUE ~ "FN",
      pred_value_inv <= alpha_cutoff & !!pos_col == FALSE ~ "TN"
    )) %>% 
    as.data.table()
  
  return(df)
  
}

# DEPRECATED:assign prediction class labels for one predictor
# addOneLabel <- function(df, cutoff, score_col, pos_col) {
#   
#   label_name <- paste0(score_col, ".pred.class")
#   df[, label_name] <- "NA"
#   
#   df[which(!is.na(df[, ..score_col]) & df[, ..score_col] > cutoff & df[, ..pos_col]), label_name] <- "TP"
#   df[which(!is.na(df[, ..score_col]) & df[, ..score_col] <= cutoff & !df[, ..pos_col]), label_name] <- "TN"
#   df[which(!is.na(df[, ..score_col]) & df[, ..score_col] > cutoff & !df[, ..pos_col]), label_name] <- "FP"
#   df[which(!is.na(df[, ..score_col]) & df[, ..score_col] <= cutoff & df[, ..pos_col]), label_name] <- "FN"
#   
#   return(df)
# }

# compute the number of positives for a given cell type (default: combined == all cells)
calcPctPos <- function(df, pos_col, cell_type = "combined") {
  
  # extract data for the provided cell type
  df_ct <- getCellTypeData(df, cell_type = cell_type)
  
  # compute percentage of positives
  mean(df_ct[[pos_col]])

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

# get the number of facet rows and columns from a ggplot object
get_row_col <- function(p) {
  n <- length(unique(ggplot_build(p)$data[[1]]$PANEL))
  par <- ggplot_build(p)$layout$facet$params
  wrap_dims(n, par$nrow, par$ncol)
}

# get default alpha values for all predictors in pred_config that have alpha == NA
getDefaultAlpha <- function(pred_config, merged) {
  
  # get predictors with missing alpha
  missing_alpha <- filter(pred_config, is.na(alpha))
  
  # get default alpha for these (if there are any)
  if (nrow(missing_alpha) > 0) {
    alpha <- mapply(FUN = get_alpha_min, predictor = missing_alpha$pred_uid,
                    inverse = missing_alpha$inverse_predictor, MoreArgs = list(merged = merged),
                    SIMPLIFY = FALSE)
    alpha <- unlist(alpha)
    
    # set alpha in pred_config to default alpha for these predictors
    pred_config$alpha <- replace(pred_config$alpha, list = is.na(pred_config$alpha), values = alpha)
    
  }
  
  return(pred_config)
  
}

# get minimum (or maximum) predictor value from merged data for default alpha values
get_alpha_min <- function(merged, predictor, inverse = FALSE) {
  merged %>% 
    filter(pred_uid == predictor, Prediction == 1) %>% 
    summarize(alpha = if_else(inverse == TRUE, max(pred_value), min(pred_value))) %>% 
    pull(alpha)
}

## DEPRECATED ======================================================================================

# # make jitter plots showing scores for each predictor as a function of experimental outcome
# plotPredictorsVsExperiment <- function(df, pos_col = "Regulated", pred_names_col = "pred_uid",
#                                        point_size = 2, text_size = 13, cell_type = "combined") {
#   
#   # get data for specified cell type
#   df_ct <- getCellTypeData(df, cell_type = cell_type)
#   
#   # plot scores for each predictor as a function of experimental outcome
#   ggplot(df_ct, aes(x = get(pos_col), y = pred_value, color = get(pos_col))) +
#     facet_wrap(~get(pred_names_col), scales = "free") +
#     geom_jitter(size = point_size) +
#     geom_boxplot(fill = NA, outlier.shape = NA, color = "black") +
#     scale_color_manual(values = c("darkgray", "steelblue")) +
#     labs(title = "Predictors vs experimental outcome", x = "Experimental E-G pair",
#          y = "Predictor value") +
#     scale_y_log10() +
#     theme_bw() +
#     theme(legend.position = "none", text = element_text(size = text_size))
#   
# }
