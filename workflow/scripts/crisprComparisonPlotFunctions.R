## Plotting functions for CRISPR comparison

library(tidyverse)
library(data.table)
library(ggpubr)
library(ggcorrplot)

## WORK IN PROGRESS CODE ===========================================================================

# calculate performance metrics using ROCR
calculatePerformance <- function(merged, pos_col, pred_config, measure, x.measure = NULL) {
  
  # split into list for lapply
  merged_split <- split(merged, f = merged$pred_uid)
  
  # get inverse predictors
  inverse_predictors <- pred_config %>% 
    select(pred_uid, inverse_predictor) %>% 
    deframe()
  
  # multiply inverse predictors by -1 so that higher value corresponds to higher score
  inverse_predictors <- inverse_predictors[names(merged_split)]  # same as predictors for cell type
  merged_split <- mapply(FUN = function(pred, inv_pred) {
    inv_multiplier <- ifelse(inv_pred, -1, 1)
    pred$pred_value <- pred$pred_value * inv_multiplier
    return(pred)
  }, merged_split, inverse_predictors, SIMPLIFY = FALSE)
  
  # compute precision-recall performance for each predictor
  perf <- lapply(merged_split, FUN = function(p){
    performance(prediction(p$pred_value, p[[pos_col]]), measure = measure, x.measure = x.measure)
  })
  
  return(perf)
  
}

# convert a list of ROCR performance objects into a table
# TODO: get rid of manual names and use internal names for ROCR performance metrics
perfToTable <- function(perf_list, measure_name, x.measure_name) {
  
  # function to convert one performance object to a table
  convert_perfToTable <- function(perf, measure_name, x.measure_name) {
    df <- data.frame(list(
      alpha = perf@alpha.values[[1]],
      measure = perf@y.values[[1]],
      x.measure = perf@x.values[[1]]
    ))
    return(df)
  }
  
  # apply to input list
  perf_list <- lapply(perf_list, convert_perfToTable)
  
  # convert list of tables into one table
  perf_table <- rbindlist(perf_list, idcol = "pred_uid")
  colnames(perf_table)[colnames(perf_table) == "measure"] <- measure_name
  colnames(perf_table)[colnames(perf_table) == "x.measure"] <- x.measure_name
  
  return(perf_table)
  
}

# make a eceiver operating characteristic (ROC) curve
plotROC <- function(merged, pos_col, pred_config, colors, thresholds = NULL,
                    plot_name = "ROC curve full experimental data", line_width = 1, 
                    point_size = 3, text_size = 15) {
  
  # compute ROC curve
  roc <- calculatePerformance(merged, pos_col = pos_col, pred_config = pred_config, measure = "tpr",
                              x.measure = "fpr")
  
  # make ROC table
  roc <- perfToTable(roc, measure_name = "TPR", x.measure_name = "FPR")
  
  # add pretty predictor names to pr_df for plotting
  roc <- left_join(roc, select(pred_config, pred_uid, pred_name_long), by = "pred_uid")
  
  # separate pr data into quantitative and boolean predictors
  bool_preds <- pull(filter(pred_config, boolean == TRUE), pred_uid)
  roc_quant <- filter(roc, !pred_uid %in% bool_preds)
  roc_bool  <- filter(roc, pred_uid %in% bool_preds)
  
  # get precision and recall for boolean predictor at alpha 1
  roc_bool <- filter(roc_bool, alpha == 1)

  # create PRC plot (caution, this assumes that there at least 1 quant and 1 bool predictor!)
  ggplot(roc_quant, aes(x = FPR, y = TPR, color = pred_name_long)) +
    geom_line(size = line_width) +
    geom_point(data = roc_bool, size = point_size) +
    labs(title = plot_name, x  = "False postitive rate", y = "True postitive rate",
         color = "Predictor") + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    scale_color_manual(values = colors, breaks = names(colors)) +
    theme_bw() +
    theme(text = element_text(size = text_size))
  
}

# make performance (AUPRC) per subset plot
plotPerfSubsets <- function(perf, pred_config, subset_name = NULL, title = NULL, order = NULL) {
  
  # create prettier labels for subsets in plots
  perf <- perf %>% 
    mutate(subset_label = paste0(subset, "\n", pos, "\n", neg)) %>% 
    mutate(subset_label = fct_inorder(subset_label))
  
  # add pretty names for predictors
  perf <- left_join(perf, select(pred_config, pred_uid, pred_name_long), by = c("id" = "pred_uid"))
  
  # order predictors according AUPRC on whole dataset or based on provided order
  if (is.null(order)) {
    perf <- perf %>% 
      filter(subset == "All") %>% 
      select(id, full_all = full) %>% 
      left_join(perf, .,  by = "id") %>%
      mutate(pred_name_long = fct_reorder(pred_name_long, .x = full_all, .desc = TRUE))
  } else {
    perf <- perf %>% 
      mutate(pred_name_long = factor(pred_name_long, levels = order))
  }
  
  # get color for each predictor
  pred_colors <- deframe(select(pred_config, pred_name_long, color))
  
  # create title and x axis label
  if (is.null(subset_name)) subset_name <- unique(perf$subset_col)
  x_label <- paste0(subset_name, "\nCRISPRi positives\nCRISPRi negatives")
  if (is.null(title)) title <- paste(unique(perf$metric), "vs.", subset_name)
  
  # make performance as function of distance plot
  ggplot(perf, aes(x = subset_label, y = full, fill = pred_name_long)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.9),
                  color = "black", width = 0.25) +
    labs(fill = "Predictor", x = x_label, y = unique(perf$metric), title = title) +
    scale_fill_manual(values = pred_colors[levels(perf$pred_name_long)]) + 
    scale_y_continuous(limits = c(0, 1)) +
    theme_bw() +
    theme(legend.position = "bottom")
  
}

# compute bootstrapped performance on subsets (and whole dataset if all == TRUE)
computePerfSubsets <- function(merged, pred_config, subset_col, metric = c("auprc", "precision"),
                               thresholds = NULL, pos_col = "Regulated",
                               bs_iter = 1000, all = TRUE) {
  
  # process metric argument
  metric <- match.arg(metric)
  
  # compute and bootstrap performance on subsets
  perf_subsets <- merged %>%
    split(., f = merged[[subset_col]]) %>% 
    lapply(FUN = convertMergedForBootstrap, pred_config = pred_config, pos_col = pos_col) %>% 
    lapply(FUN = bootstrapPerformanceIntervals, metric = metric, thresholds = thresholds,
           R = bs_iter) %>% 
    bind_rows(.id = "subset")
  
  # compute and bootstrap performance on entire dataset
  if (all == TRUE) {
    perf_subsets <- merged %>% 
      convertMergedForBootstrap(pred_config = pred_config, pos_col = pos_col) %>% 
      bootstrapPerformanceIntervals(metric = metric, thresholds = thresholds, R = bs_iter) %>% 
      mutate(subset = "All", .before = 1) %>% 
      bind_rows(., perf_subsets)
  }
  
  # count the number of positive and negative E-G pair for each subset and add to performance table
  pairs_subsets <- count_pairs_subset(merged, subset_col = subset_col, pos_col = pos_col, all = all)
  perf_subsets <- left_join(perf_subsets, pairs_subsets, by = "subset")
  
  # add used subset column to table
  perf_subsets <- mutate(perf_subsets, subset_col = subset_col, .before = 1)
  
  return(perf_subsets)
  
}

# count the number of crispr positive and negatives in subsets (and whole dataset if all == TRUE)
count_pairs_subset <- function(merged, subset_col, pos_col, all = TRUE) {
  
  # unique experimentally tested E-G pairs
  crispr_pairs <- merged %>% 
    select(name, subset = all_of(subset_col), positive = all_of(pos_col)) %>% 
    distinct()
  
  # count number of positive and negatives in each subset
  pairs_subsets <- crispr_pairs %>% 
    group_by(subset) %>% 
    summarize(pos = sum(positive == TRUE),
              neg = sum(positive == FALSE))
  
  # add the number of positive and negatives in the entire dataset if specified
  if (all == TRUE) {
    
    pairs_subsets <- crispr_pairs %>% 
      summarize(pos = sum(positive == TRUE),
                neg = sum(positive == FALSE)) %>% 
      mutate(subset = "All", .before = 1) %>% 
      bind_rows(., pairs_subsets)
    
  }
  
  return(pairs_subsets)
  
}


## MAIN FUNCTIONS ==================================================================================

# process merged data for benchmarking analyses
processMergedData <- function(merged, pred_config, filter_valid_connections = TRUE,
                              include_missing_predictions = TRUE, distToTSS_as_kb = TRUE) {
  
  # add unique identifiers to merged data
  merged$pred_uid <- paste(merged$pred_id, merged$pred_col, sep = ".")
  
  # only retain predictors that are specified to be included in plots
  plot_preds <- pred_config[pred_config$plot_crispr == TRUE, ][["pred_uid"]]
  merged <- merged[merged$pred_uid %in% plot_preds, ]
  
  # add long names and whether a predictor is boolean from pred_config to merged data
  merged <- pred_config %>% 
    select(pred_uid, boolean, pred_name_long) %>% 
    left_join(x = merged, y = ., by = "pred_uid")
  
  # filter merged data for valid connections if specified
  if (filter_valid_connections == TRUE) {
    merged <- subset(merged, ValidConnection == "TRUE")
  }
  
  # filter out CRE - gene pairs with missing predictions if specified
  if (include_missing_predictions == FALSE) {
    merged <- merged[merged$Prediction == 1, ]
  }
  
  # convert distance to TSS baseline predictor to kb if specified
  if (distToTSS_as_kb == TRUE) {
    merged <- merged %>% 
      mutate(pred_value = if_else(pred_uid == "baseline.distToTSS", true = pred_value / 1000,
                                  false = pred_value)) %>% 
      mutate(pred_name_long = if_else(pred_uid == "baseline.distToTSS",
                                      true = "Distance to TSS (kb)", false = pred_name_long))
  }
  
  return(merged)
  
}

# process pred_config file for benchmarking analyses
processPredConfig <- function(pred_config, merged) {
  
  # pred_config column names relevant for crispr benchmarking
  config_cols <- c("pred_id", "pred_col", "boolean", "alpha", "aggregate_function", "fill_value",
                   "inverse_predictor", "pred_name_long", "color", "plot_crispr")
  
  # initialize optional plot_crispr column with 'TRUE', if not present in config table
  if (!"plot_crispr" %in% colnames(pred_config)) {
    pred_config$plot_crispr <- TRUE
  }
  
  # only retain columns relevant for CRISPR benchmarking
  pred_config <- pred_config[, ..config_cols]
  
  # add unique predictor identifier to pred_config
  pred_config$pred_uid <- paste(pred_config$pred_id, pred_config$pred_col, sep = ".")
  
  # filter pred_config for predictors also occurring in merged data
  merged_pred_uid <- unique(paste(merged$pred_id, merged$pred_col, sep = "."))
  pred_config <- filter(pred_config, pred_uid %in% merged_pred_uid)
  
  # create default baseline predictor configuration
  baseline_preds <- c("distToTSS", "nearestTSS", "nearestGene", "within100kbTSS", "nearestExprTSS",
                      "nearestExprGene", "within100kbExprTSS")
  baseline_preds_uid <- paste("baseline", baseline_preds, sep = ".")
  baseline_pred_config <- data.table(
    pred_uid = baseline_preds_uid,
    pred_id = "baseline",
    pred_col = baseline_preds,
    boolean = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    alpha = c(NA, 1, 1, 1, 1, 1, 1),
    aggregate_function = c("mean", "max", "max", "max", "max", "max", "max"),
    fill_value = c(Inf, 0, 0, 0, 0, 0, 0),
    inverse_predictor = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
    pred_name_long = c("Distance to TSS", "Nearest TSS", "Nearest Gene", "Within 100kb of TSS",
                       "Nearest expr. TSS", "Nearest expr. Gene", "Within 100kb of expr. TSS"),
    color = c("#ffa600", "#595959", "#bebebe", "#000000", "#595959", "#bebebe", "#000000"),
    plot_crispr = TRUE
  )
  
  # only retain default values for baseline predictors also found in merged data
  baseline_pred_config <- filter(baseline_pred_config, pred_uid %in% merged_pred_uid)
  
  # if no colors were set for any of the other predictors, set colors of baseline predictors to NA
  if (all(is.na(pred_config$color))) {
    baseline_pred_config$color <- NA_character_
  }
  
  # config for baseline predictors can be set in pred_config, so only add the default values for
  # those missing from pred_config
  baseline_preds_to_add <- setdiff(baseline_pred_config$pred_uid, pred_config$pred_uid)
  
  # if both nearest TSS/gene and nearest expr. TSS/gene baseline predictors where created, only
  # include nearest expressed versions in plots by default. this can be overridden by specifying
  # predictors to plot using the plot_crispr column
  if (all(c("baseline.nearestTSS", "baseline.nearestExprTSS") %in% baseline_preds_to_add)) {
    baseline_preds_to_add <- setdiff(baseline_preds_to_add, "baseline.nearestTSS")
  }
  
  if (all(c("baseline.nearestGene", "baseline.nearestExprGene") %in% baseline_preds_to_add)) {
    baseline_preds_to_add <- setdiff(baseline_preds_to_add, "baseline.nearestGene")
  }
  
  if (all(c("baseline.within100kbTSS", "baseline.within100kbExprTSS") %in% baseline_preds_to_add)) {
    baseline_preds_to_add <- setdiff(baseline_preds_to_add, "baseline.within100kbTSS")
  }
  
  # add baseline predictors to pred_config
  pred_config <- baseline_pred_config %>% 
    filter(pred_uid %in% baseline_preds_to_add) %>% 
    rbind(pred_config, .)
  
  # only retain predictors that are specified to be included in plots
  pred_config <- pred_config[pred_config$plot_crispr == TRUE, ]
  
  # check that long predictor names are unique
  if (any(table(pred_config$pred_name_long) > 1)) {
    stop("'pred_name_long' in pred_config is not a unique identifier.", call. = FALSE)
  }
  
  return(pred_config)
  
}

# plot CRISPR E-G pairs overlapping prediction E-G pairs
plotOverlaps <- function(merged, title = "E-G pairs in predictions part of CRISPR E-G universe") {
  
  # count number of total CRISPR E-G pairs and overlapping pairs per predictor
  n_pairs <- merged %>% 
    group_by(pred_uid, pred_name_long) %>% 
    summarize(`Overlaps predictions` = sum(Prediction == 1),
              `Not in predictions` = sum(Prediction == 0),
              .groups = "drop") %>% 
    pivot_longer(cols = c(`Overlaps predictions`, `Not in predictions`), names_to = "Overlaps",
                 values_to = "pairs")
  
  # plot number of CRISPR E-G pairs overlapping E-G pairs in predictions
  ggplot(n_pairs, aes(x = pred_name_long, y = pairs, fill = Overlaps)) +
    geom_bar(stat = "identity") +
    labs(y = "E-G pairs", x = "Predictor", title = title) +
    scale_fill_manual(values = c("Overlaps predictions" = "steelblue",
                                 "Not in predictions" = "darkgray")) +
    coord_flip() +
    theme_bw()
  
}

# compute PR curves for merged data
calcPRCurves <- function(df, pred_config, pos_col) {
  
  # split into list for lapply
  df_split <- split(df, f = df$pred_uid)
  
  # get inverse predictors
  inverse_predictors <- pred_config %>% 
    select(pred_uid, inverse_predictor) %>% 
    deframe()
  
  # multiply inverse predictors by -1 so that higher value corresponds to higher score
  inverse_predictors <- inverse_predictors[names(df_split)]  # same as predictors for cell type
  df_split <- mapply(FUN = function(pred, inv_pred) {
    inv_multiplier <- ifelse(inv_pred, -1, 1)
    pred$pred_value <- pred$pred_value * inv_multiplier
    return(pred)
  }, df_split, inverse_predictors, SIMPLIFY = FALSE)
  
  # compute precision-recall performance for each predictor
  pr <- lapply(df_split, FUN = function(p){
    performance(prediction(p$pred_value, p[[pos_col]]), measure = "prec", x.measure = "rec")
  })
  
  # convert to table and calculate F1
  pr_df <- pr2df(pr, calc_f1 = TRUE)
  
  return(pr_df)
  
}

# create bootstrapped performance summary table for all predictors in a PR table
makePRSummaryTableBS <- function(merged, pred_config, pos_col, min_sensitivity = 0.7, R = 1000,
                                 conf = 0.95, ncpus = 1) {
  
  # convert merged to wide format for bootstrapping
  merged_bs <- convertMergedForBootstrap(merged, pred_config = pred_config, pos_col = pos_col)
  
  # extract defined thresholds
  thresholds <- deframe(select(pred_config, pred_uid, alpha))
  thresholds <- thresholds[!is.na(thresholds) & names(thresholds) %in% colnames(merged_bs)]
  
  # bootstrap overall performance (AUPRC) and reformat for performance summary table
  perf <- bootstrapPerformanceIntervals(merged_bs, metric = "auprc", R = R, conf = conf,
                                        ci_type = "perc", ncpus = ncpus)
  perf <- select(perf, pred_uid = id, AUPRC = full, AUPRC_lowerCi = lower, AUPRC_upperCi = upper)
  
  # bootstrap precision at defined thresholds (if there are any)
  if (length(thresholds) > 0) {
    
    # get data on predictors with thresholds
    merged_bs_thresh <- select(merged_bs, name, Regulated, any_of(names(thresholds)))
    
    # run precision bootstraps using defined thresholds 
    prec_thresh <- bootstrapPerformanceIntervals(merged_bs_thresh, metric = "precision",
                                                 thresholds = thresholds, R = R, conf = conf,
                                                 ci_type = "perc", ncpus = ncpus)
    # add to performance table
    perf <- prec_thresh %>% 
      select(pred_uid = id, PrecThresh = full, PrecThresh_lowerCi = lower,
             PrecThresh_upperCi = upper) %>% 
      left_join(enframe(thresholds, name = "pred_uid", value = "threshold"),
                by = "pred_uid") %>% 
      left_join(perf, ., by = "pred_uid")
    
  }
  
  # get performance at minimum sensitivity
  pr <- calcPRCurves(merged, pred_config = pred_config, pos_col = pos_col)
  perf_min_sens <- pr %>% 
    arrange(pred_uid, recall, desc(precision)) %>% 
    split(., f = .$pred_uid) %>% 
    lapply(FUN = computePerfGivenSensitivity, min_sensitivity = min_sensitivity) %>% 
    bind_rows(.id = "pred_uid")
  
  # extract threshold at min sensitivity
  thresholds_min_sens <- deframe(select(perf_min_sens, pred_uid, alpha_at_min_sensitivity))
  
  # bootstrap precison at minimum sensitivity
  prec_min_sens <- bootstrapPerformanceIntervals(merged_bs, metric = "precision",
                                                 thresholds = thresholds_min_sens, R = R,
                                                 conf = conf, ci_type = "perc", ncpus = ncpus)
  
  # add to performance table
  perf <- prec_min_sens %>% 
    select(pred_uid = id, PrecMinSens = full, PrecMinSens_lowerCi = lower,
           PrecMinSens_upperCi = upper) %>% 
    left_join(enframe(thresholds_min_sens, name = "pred_uid", value = "thresholdMinSens"),
              by = "pred_uid") %>% 
    left_join(perf, ., by = "pred_uid")
  
  # sort accoring to overall performance for output
  perf <- arrange(perf, desc(AUPRC))
    
  return(perf)
  
}

# make a PR curve plot for a set of provided predictors
makePRCurvePlot <- function(pr_df, pred_config, n_pos, pct_pos, min_sensitivity = 0.7,
                            plot_name = "PRC full experimental data", line_width = 1, 
                            point_size = 3, text_size = 15, colors = NULL, plot_thresholds = TRUE,
                            na_color = "gray66", na_size_factor = 0.5) {
  
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
  
  # create default colors if none were specified
  if (is.null(colors)) {
    colors <- scales::hue_pal()(nrow(pred_config))
  }
  
  # separate predictors into those with a set color and those with NA (useful against overplotting)
  na_col <- names(colors)[is.na(colors)]
  pr_quant_col <- filter(pr_quant, !pred_name_long %in% na_col)
  pr_quant_na_col <- filter(pr_quant, pred_name_long %in% na_col)
  pr_bool_col <- filter(pr_bool, !pred_name_long %in% na_col)
  pr_bool_na_col <- filter(pr_bool, pred_name_long %in% na_col)
  pr_threshold_col <- filter(pr_threshold, !pred_name_long %in% na_col)
  pr_threshold_na_col <- filter(pr_threshold, pred_name_long %in% na_col)
  
  # set NA colors to a lighter gray than the ggplot default
  colors[is.na(colors)] <- na_color
  
  # create PRC plot (caution, this assumes that there at least 1 quant and 1 bool predictor!)
  p <- ggplot(pr_quant, aes(x = recall, y = precision, color = pred_name_long)) +
    geom_line(data = pr_quant_na_col, linewidth = line_width * na_size_factor)
  
  if (plot_thresholds) p <- p + geom_point(data = pr_threshold_na_col,
                                           size = point_size * na_size_factor)
  
  p <- p + 
    geom_point(data = pr_bool_na_col, size = point_size * na_size_factor) +
    geom_line(data = pr_quant_col, linewidth = line_width)
  
  if (plot_thresholds) p <- p + geom_point(data = pr_threshold_col, size = point_size)
  
  p + 
    geom_point(data = pr_bool_col, size = point_size) +
    geom_hline(yintercept = pct_pos, linetype = "dashed", color = "black") +
    labs(title = plot_name, x  = paste0("Recall (n=", n_pos, ")"), y = "Precision",
         color = "Predictor") + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + 
    scale_color_manual(values = colors, breaks = names(colors)) +
    theme_bw() +
    theme(text = element_text(size = text_size))
  
}

# plot predictor scores vs effect sizes
plotPredictorsVsEffectSize <- function(merged, pos_col = "Regulated", pred_names_col = "pred_uid",
                                       corr_groups = pos_col, point_size = 2, text_size = 13,
                                       title = "Predictor scores vs. CRISPRi effect sizes",
                                       alpha_value = 1, include_boolean_preds = FALSE,
                                       cor_method = "spearman", ncol = NULL, nrow = NULL,
                                       label.x.npc = 0.75) {
  
  # remove boolean predictors if specified
  if (include_boolean_preds == FALSE) {
    if ("boolean" %in% colnames(merged)) {
      merged <- subset(merged, boolean == FALSE)
    } else {
      warning("'boolean' column not in 'merged' data frame, can't filter for boolean predictors",
              call. = FALSE)
    }
  }
  
  # set color based on column identifying positive CRISPRi pairs
  if (!is.null(pos_col)) {
    values <- sort(unique(merged[[pos_col]]))
    colors <- structure(c("darkgray", "steelblue"), names = as.character(values))
  } else {
    colors <- NULL
  }
  
  
  # create basic effect size vs predictors scatter plots
  p <- makeScatterPlots(merged, x_col = "pred_value", y_col = "EffectSize", color_col = pos_col,
                        x_lab = "Predictor score", y_lab = "CRISPRi effect size",
                        title = title, colors = colors,
                        pred_names_col = pred_names_col, point_size = point_size,
                        text_size = text_size, alpha_value = alpha_value, ncol = ncol, nrow = nrow)
  
  # add linear model fit and correlation coefficient
  if (is.null(corr_groups)) {
    p + 
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black") +
      stat_cor(aes(color = NULL, label = ..r.label..), method = cor_method, cor.coef.name = "rho",
               label.y.npc = "top", label.x.npc = label.x.npc, r.accuracy = 0.01,
               show.legend = FALSE)
  } else {
    p + 
      geom_smooth(aes(color = !!sym(corr_groups)), method = "lm", formula = y ~ x, se = FALSE) +
      stat_cor(aes(label = ..r.label..), method = cor_method, cor.coef.name = "rho",
               label.y.npc = "top", label.x.npc = label.x.npc, r.accuracy = 0.01,
               show.legend = FALSE)
  }
  
}

# make violin plots showing scores for each predictor as a function of experimental outcome
plotPredictorsVsExperiment <- function(merged, pos_col = "Regulated", pred_names_col = "pred_uid",
                                       text_size = 13, include_boolean_preds = FALSE, ncol = NULL,
                                       nrow = NULL) {
  
  # remove boolean predictors if specified
  if (include_boolean_preds == FALSE) {
    if ("boolean" %in% colnames(merged)) {
      merged <- subset(merged, boolean == FALSE)
    } else {
      warning("'boolean' column not in 'merged' data frame, can't filter for boolean predictors",
              call. = FALSE)
    }
  }
  
  # set color based on column identifying positive CRISPRi pairs
  values <- sort(unique(merged[[pos_col]]))
  colors <- structure(c("darkgray", "steelblue"), names = as.character(values))
  
  # convert column identifying positive CRISPR hits from string to symbol for ggplot
  pos_col <- sym(pos_col)
  
  # plot scores for each predictor as a function of experimental outcome
  ggplot(merged, aes(x = !!pos_col, y = pred_value, color = !!pos_col, fill = !!pos_col)) +
    facet_wrap(~get(pred_names_col), scales = "free", ncol = ncol, nrow = nrow) +
    geom_violin() +
    geom_boxplot(width = 0.1, outlier.shape = NA, color = "black", fill = "NA") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(title = "Predictors vs experimental outcome", x = "Experimental E-G pair",
         y = "Predictor value (sqrt)") +
    scale_y_sqrt() +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = text_size))
  
}

# plot distance distributions for experimental positives and negatives
plotDistanceDistribution <- function(merged, dist = "baseline.distToTSS", pos_col = "Regulated",
                                     text_size = 13) {
  
  # get data for distance and add label for faceting
  dist_data <- merged %>% 
    filter(pred_uid == dist) %>% 
    mutate(label = if_else(get(pos_col) == TRUE, true = "Positives", false = "Negatives")) %>% 
    mutate(label = factor(label, levels = c("Positives", "Negatives")))
  
  # plot distance distribution for all pairs in experimental data
  ggplot(dist_data, aes(x = pred_value, fill = label)) +
    facet_wrap(~label, ncol = 1, scales = "free_y") +
    geom_histogram(binwidth = 10) +
    labs(x = "Distance to TSS (kb)") +
    scale_fill_manual(values = c("Positives" = "steelblue", "Negatives" = "darkgray")) +
    theme_bw() +
    theme(legend.position = "none", text = element_text(size = text_size))
  
}

# make an upset plot of overlapping features for a given cell type
plotOverlappingFeatures <- function(merged, feature_cols) {
  
  # create table with unique enhancers and overlapping features
  enh_features <- merged %>% 
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

# plot correlation between predictors
plotPredCorMatrix <- function(merged, pred_names_col = "pred_uid", include_boolean_preds = FALSE,
                              method = c("spearman", "pearson", "kendall"),
                              title = "Correlation of predictor scores") {
  
  # parse method argument
  method <- match.arg(method)
  
  # remove boolean predictors if specified
  if (include_boolean_preds == FALSE) {
    if ("boolean" %in% colnames(merged)) {
      merged <- subset(merged, boolean == FALSE)
    } else {
      warning("'boolean' column not in 'merged' data frame, can't filter for boolean predictors",
              call. = FALSE)
    }
  }
  
  # extract predictor scores and transform to wide format
  pred_scores <- merged %>% 
    select(name, all_of(pred_names_col), pred_value) %>% 
    pivot_wider(names_from = all_of(pred_names_col), values_from = pred_value)
  
  # compute correlation coefficients between predictors
  cor_matrix <- cor(select(pred_scores, -name), method = method)
  
  # plot correlation matrix for scores of different predictors
  ggcorrplot(cor_matrix, hc.order = TRUE, type = "lower", lab = TRUE, title = title)
  
}

## Make plots for subsets of the data based on gene or enhancer features ---------------------------

# count the number of positive and negative pairs per gene and enhancer feature
countPairsFeatures <- function(merged, feature_cols_pattern = "^(gene|enh)_feature_.+$",
                               pos_col = "Regulated") {
  
  # get all feature columns
  feature_cols <- grep(colnames(merged), pattern = feature_cols_pattern, value = TRUE)
  
  # count the number of pairs per feature column
  output <- merged %>% 
    select(name, all_of(c(pos_col, feature_cols))) %>% 
    distinct() %>% 
    pivot_longer(cols = all_of(feature_cols), names_to = "feature", values_to = "value",
                 values_transform = as.character) %>% 
    group_by(feature, value) %>%
    summarize("TRUE" = sum(get(pos_col) == TRUE),
              "FALSE" = sum(get(pos_col) == FALSE),
              .groups = "drop") %>% 
    pivot_longer(cols = c("TRUE", "FALSE"), names_to = "outcome", values_to = "pairs") %>% 
    separate(col = feature, into = c("feature_type", "feature"), sep = "_feature_")
  
  return(output)
  
}

# plot the number of pairs grouped by feature types (if feature_types has names, the names will be
# used for pretty plots)
plotPairsFeatures <- function(n_pairs_features,
                              feature_types = c("Gene" = "gene", "Enhancer" = "enh")) {
  
  # split n_pairs_features by feature types
  n_pairs_features <- split(n_pairs_features, f = n_pairs_features$feature_type)
  
  # order and rename feature types (if specified)
  n_pairs_features <- n_pairs_features[feature_types]
  if (!is.null(names(feature_types))) names(n_pairs_features) <- names(feature_types)
  
  # create titles for plots
  titles <- paste(names(n_pairs_features), "features")
  
  # plot number of pairs per feature type
  feat_pair_plots <- mapply(FUN = function(set, title) {
    ggplot(set, aes(x = value, y = pairs, fill = as.logical(outcome))) +
      facet_wrap(~feature) +
      geom_bar(stat = "identity") +
      labs(title = title, y = "Number of E-G pairs", fill = pos_col) +
      scale_fill_manual(values = c("FALSE" = "darkgray", "TRUE" = "steelblue")) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1))
  }, set = n_pairs_features, title = titles, SIMPLIFY = FALSE)
  
  return(feat_pair_plots)
  
}

# calculate and plot PR curves for a given subset
makePRCurveSubset <- function(merged, subset_col, pred_config, pos_col, min_sensitivity = 0.7,
                              line_width = 1, point_size = 3, text_size = 15, nrow = 1,
                              colors = NULL) {
  
  # split df into subsets based on provided column
  merged_split <- split(merged, f = merged[[subset_col]])
  subsets <- names(merged_split)  # used later
  
  # get subsets that do not contain both positives and negatives
  one_class_only <- vapply(merged_split, FUN.VALUE = logical(1), FUN = function(x){
    n_distinct(x$Regulated) != 2 
  })
  
  # filter out any subsets that do not have both
  if (all(one_class_only)) {
    stop("No subsets contain both positive and negatives.", call. = FALSE)
  } else if (any(one_class_only)) {
    warning("Not all subsets contain both positive and negatives. These will be removed from plots.",
            call. = FALSE)
    merged_split <- merged_split[!one_class_only]
  } else {
    merged_split <- merged_split
  }
  
  # compute PR curve
  prc <- lapply(merged_split, FUN = calcPRCurves, pred_config = pred_config, pos_col = pos_col)
  
  # calculate number and percentage of positives for all splits
  n_pos <- lapply(merged_split, FUN = calcNPos, pos_col = pos_col)
  pct_pos <- lapply(merged_split, FUN = calcPctPos, pos_col = pos_col)
  
  # create plot titles based on feature name and number of pairs for that feature
  feature_name <- sub(".+_feature_", "", subset_col)
  pairs <- vapply(merged_split, FUN = function(x) length(unique(x$name)), FUN.VALUE = integer(1))
  titles <- paste0(feature_name, " = " , names(prc), " (", pairs, " pairs)")
  
  # plot PR curves
  pr_plots <- mapply(FUN = makePRCurvePlot, pr_df = prc, n_pos = n_pos, pct_pos = pct_pos, plot_name = titles, 
                     MoreArgs = list(pred_config = pred_config, min_sensitivity = min_sensitivity,
                                     line_width = line_width, point_size = point_size,
                                     text_size = text_size, colors = colors),
                     SIMPLIFY = FALSE)
  
  # add empty plots for subsets without both positives and negatives
  empty <- setdiff(subsets, names(pr_plots))
  for (i in empty) pr_plots[[i]] <- ggplot(NULL)
  
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
makePRCurveSubsets <- function(merged, subset_cols, pred_config, pos_col, cell_type = "combined",
                               min_sensitivity = 0.7, line_width = 1, point_size = 3,
                               text_size = 15, nrow = 1, colors = NULL) {
  
  # return NULL if no subsets are available
  if (length(subset_cols) == 0) {
    return(NULL)
  }
  
  # create PR curve plots for each subset
  pr_plots <- lapply(subset_cols, FUN = makePRCurveSubset, merged = merged,
                     pred_config = pred_config, pos_col = pos_col,
                     min_sensitivity = min_sensitivity, line_width = line_width,
                     point_size = point_size, text_size = text_size, nrow = nrow, colors = colors)
  
  # create one figure with all plots
  plot_grid(plotlist = pr_plots, nrow = length(subset_cols))
  
}

# make violin plots showing scores for each predictor vs experimental outcome for a given subset
plotPredVsExperimentSubset <- function(merged, subset_col, pos_col, pred_names_col = "pred_uid",
                                       text_size = 13) {
  
  # plot scores for each predictor as a function of experimental outcome
  p <- plotPredictorsVsExperiment(merged, pos_col = pos_col, pred_names_col = pred_names_col,
                                  text_size = text_size)
  
  # change title and add faceting based on subset column
  p +
    labs(title = sub(".+_feature_", "", subset_col)) +
    facet_grid(get(pred_names_col) ~ get(subset_col), scales = "free")
  
}

# make violin plots showing scores for each predictor vs experimental outcome for a given subset
plotPredVsExperimentSubsets <- function(merged, subset_cols, pos_col, pred_names_col = "pred_uid",
                                        text_size = 13, label.x.npc = 0.5) {
  
  # return NULL if no subsets are available
  if (length(subset_cols) == 0) {
    return(NULL)
  }
  
  # create predictor vs experiment plots for each subset
  pe_plots <- lapply(subset_cols, FUN = plotPredVsExperimentSubset, merged = merged,
                     pos_col = pos_col, pred_names_col = pred_names_col, text_size = text_size)
  
  # create one figure with all plots
  plot_grid(plotlist = pe_plots, nrow = length(subset_cols))
  
}

# make violin plots showing scores for each predictor vs experimental outcome for a given subset
predVsEffectSizeSubset <- function(merged, subset_col, pos_col, pred_names_col = "pred_uid",
                                   corr_groups = pos_col, point_size = 2, text_size = 13,
                                   alpha_value = 1, label.x.npc = 0.75) {
  
  # plot each predictor against effect size for the given subset
  p <- plotPredictorsVsEffectSize(merged, pos_col = pos_col, pred_names_col = pred_names_col,
                                  corr_groups = corr_groups, point_size = point_size,
                                  text_size = text_size, alpha_value = alpha_value,
                                  label.x.npc = label.x.npc)
  
  # change title and add faceting based on subset column
  p +
    labs(title = sub(".+_feature_", "", subset_col)) +
    facet_grid(get(subset_col) ~ get(pred_names_col), scales = "free")
  
}

# make violin plots showing scores for each predictor vs experimental outcome for a given subset
predVsEffectSizeSubsets <- function(merged, subset_cols, pos_col, pred_names_col = "pred_uid",
                                    corr_groups = pos_col, point_size = 2, text_size = 13,
                                    alpha_value = 1, label.x.npc = 0.75) {
  
  # return NULL if no subsets are available
  if (length(subset_cols) == 0) {
    return(NULL)
  }
  
  # create predictor vs experiment scatter plots for each subset
  es_plots <- lapply(subset_cols, FUN = predVsEffectSizeSubset, merged = merged, pos_col = pos_col,
                     pred_names_col = pred_names_col, corr_groups = corr_groups,
                     point_size = point_size, text_size = text_size, alpha_value = alpha_value,
                     label.x.npc = label.x.npc)
  
  # create one figure with all plots
  plot_grid(plotlist = es_plots, nrow = length(subset_cols))
  
}

## HELPER FUNCTIONS ================================================================================

# apply analyses across all cell types in merged data ----------------------------------------------

# get all cell types in merged data. if there are multiple cell types and if combined = TRUE add
# "combined" to include a set of all cell types combined
getCellTypes <- function(merged, cell_types_col, combined) {
  
  # get all cell types in merged data (remove any NAs)
  cell_types <- unique(merged[[cell_types_col]])
  cell_types <- cell_types[!is.na(cell_types)]
  cell_types <- structure(as.list(cell_types), names = cell_types)
  
  # add combined set if specified
  if (length(cell_types) > 1 & combined == TRUE) {
    cell_types <- c(cell_types, list(combined = as.character(cell_types)))
  }
  
  return(cell_types)
  
}

# get data for a specified cell type or all if cell_type == "combined"
getCellTypeData <- function(df, cell_type, cell_types_col) {
  df[df[[cell_types_col]] %in% cell_type, ]
}

# helper function to apply a function to one cell type in merged data
applyCellType <- function(cell_type, merged, .fun, ..., cell_types_col) {
  
  # get data for the given cell type
  merged_cell_type <- getCellTypeData(merged, cell_type = cell_type, cell_types_col = cell_types_col)
  
  # try to apply the specified function and capture any errors and warnings
  output <- tryCatch(
    withCallingHandlers({
      .fun(merged_cell_type, ...)
    }, warning = function(w) {
      message("For cell type ", cell_type, ": ", w)
      invokeRestart("muffleWarning")
    }), error = function(e) {
      message("For cell type ", cell_type, ": ", e)
      return(NULL)
    })
  
  return(output)
  
}

# simple wrapper to apply a function to all cell types in merged data
applyCellTypes <- function(merged, .fun, ..., cell_types_col = "ExperimentCellType",
                           combined = TRUE, remove_failed = TRUE) {
  
  # get all cell types in merged data
  cell_types <- getCellTypes(merged, cell_types_col = cell_types_col, combined = combined)
  
  # apply function to all cell types
  output <- lapply(cell_types, FUN = applyCellType, merged = merged, .fun = .fun, ...,
                   cell_types_col = cell_types_col)
  
  # remove output for any cell types where function failed
   if (remove_failed == TRUE) {
     output <- output[!vapply(output, FUN = is.null, FUN.VALUE = logical(1))]
   }
  
  return(output)
  
}

# plotting functions -------------------------------------------------------------------------------

# make scatter plots of two columns (e.g. effect size vs predictor score) across all predictors
makeScatterPlots <- function(df, x_col, y_col, color_col = NULL, x_lab = NULL, y_lab = NULL,
                             title = NULL, colors = NULL, pred_names_col = "pred_uid",
                             point_size = 2, text_size = 13, alpha_value = 1, ncol = NULL,
                             nrow = NULL) {
  
  # create labels and title based on input parameters
  if (is.null(x_lab)) x_lab <- x_col
  if (is.null(y_lab)) y_lab <- y_col
  if (is.null(title)) title <- paste(x_col, "vs", y_col)
  
  # convert column names to symbols for ggplot
  x_col <- sym(x_col)
  y_col <- sym(y_col)
  if (!is.null(color_col)) color_col <- sym(color_col)
  
  # plot x_col vs y_col across predictors
  p <- ggplot(df, aes(x = !!x_col, y = !!y_col, color = !!color_col)) +
    facet_wrap(as.formula(paste("~", pred_names_col)), scales = "free", ncol = ncol, nrow = nrow) +
    geom_point(size = point_size, alpha = alpha_value) +
    labs(title = title, x = x_lab, y = y_lab) +
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = text_size))
  
  # set colors if provided
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }
  
  return(p)
  
}

# save a list of plots to output files (... can be any other argument for ggsave)
savePlotList <- function(plot_list, basename, path = ".", ...) {
  
  # output filenames for plots
  outfiles <- paste(rep(path, times = length(plot_list)), names(plot_list),
                    rep(basename, times = length(plot_list)), sep = "/")
  
  # create required directories if needed
  for (i in dirname(outfiles)) {
    dir.create(i, recursive = TRUE, showWarnings = FALSE)
  }
  
  # save list of plots to output files
  invisible(mapply(FUN = ggsave, outfiles, plot_list, MoreArgs = list(...)))
  
}

# print plot a plot list with a tab per plot
printTabbedPlots <- function(plots, section_level = "#", plot_function = plot) {
  for (i in names(plots)){
    cat(paste0(section_level, "#"), i, '{.unlisted .unnumbered}', '\n', '<br>', '\n')
    plot_function(plots[[i]])
    cat('\n', '<br>', '\n\n')
  }
  cat(section_level, "{.unlisted .unnumbered}")
}

# other helper functions ---------------------------------------------------------------------------

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

# compute the fraction of experimental positives
calcPctPos <- function(df, pos_col = "Regulated") {
  mean(df[[pos_col]])
}

# compute number and fractions of experimental positives
calcNPos <- function(df, pos_col = "Regulated") {
  
  expt_pairs <- df %>%
    select(chrom, chromStart, chromEnd, measuredGeneSymbol, all_of(pos_col)) %>% 
    distinct()
  
  sum(expt_pairs[[pos_col]])
  
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
  
  # set alpha to alpha at minimum sensitivity if threshold was set to NA
  perf_summary <- perf_summary %>% 
    left_join(select(pred_config, pred_uid, alpha), by = "pred_uid") %>% 
    mutate(alpha_cutoff = if_else(is.na(alpha), alpha_at_min_sensitivity, alpha_cutoff),
           alpha_at_cutoff = if_else(is.na(alpha), alpha_at_min_sensitivity, alpha_at_cutoff),
           sensitivity_at_cutoff = if_else(is.na(alpha), sensitivity_at_min_sensitivity, sensitivity_at_cutoff),
           precision_at_cutoff = if_else(is.na(alpha), precision_at_min_sensitivity, precision_at_cutoff)) %>% 
    select(-alpha)
  
  return(perf_summary)
  
}

# create performance summary for one predictor
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

