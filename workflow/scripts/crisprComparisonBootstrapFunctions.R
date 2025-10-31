## Functions to compute bootstrapped performance metrics (AUPRC and precision @ threshold) and
## significance for pairwise differences in performance between predictors

library(boot)
library(yardstick)
library(caTools)
library(data.table)
library(dplyr)
library(BiocParallel)

#' Reformat CRISPR benchmarking data for bootstrapping
#'
#' Convert merged data from CRISPR benchmarking pipeline into format required for bootstrapping
#' 
#' @param merged Merged data from CRISPR benchmarking pipeline
#' @param pred_config pred_config table used to run CRISPR benchmarking pipeline
#' @param pos_col Name of column specifying positive CRISPR E-G pairs. Needs to be a logical with 
#'   TRUE/FALSE, where TRUE are positives
#' @param weight_col Optional name for column containing case weights for calculating weighted
#'   performance
#' @param rm_boolean (logical) Should boolean predictors be filtered out? (default: TRUE)
convertMergedForBootstrap <- function(merged, pred_config, pos_col = "Regulated", weight_col = NULL,
                                      rm_boolean = TRUE) {
  
  # get inverse predictors
  inverse_preds <- pred_config[pred_config$inverse_predictor == TRUE, ][["pred_uid"]]
  
  # extract relevant columns from merged data
  merged <- merged %>% 
    select(name, Regulated = all_of(pos_col), weight = all_of(weight_col), pred_uid, pred_value)
  
  # filter out boolean predictors
  if (rm_boolean == TRUE) {
    boolean_preds <- pred_config[pred_config$boolean == TRUE, ][["pred_uid"]]
    merged <- filter(merged, !pred_uid %in% boolean_preds)
  }
  
  # multiply inverse predictors by -1
  merged <- merged %>% 
    mutate(pred_value = if_else(pred_uid %in% inverse_preds,
                                true = pred_value * -1,
                                false = pred_value))
  
  # convert to wide format to create output and convert Regulated to factor
  output <- merged %>% 
    pivot_wider(names_from = pred_uid, values_from = pred_value) %>% 
    mutate(Regulated = factor(Regulated))
  
  return(output)
  
}

#' Get predictor threshold values
#' 
#' Extract predictor thresholds from pred_config file and invert for inverse predictors.
#' 
#' @param pred_config pred_config table used to run CRISPR benchmarking pipeline.
#' @param predictors Vector with predictors for which thresholds should be extracted. Default is
#'  NULL, which will extract thresholds for all predictors if possible.
#' @param threshold_col Column name in pred_config containing threshold values (Default: alpha).
getThresholdValues <- function(pred_config, predictors = NULL, threshold_col = "alpha") {
  
  # filter for subset of predictors only if specified
  if (!is.null(predictors)) {
    pred_config <- pred_config[pred_config$pred_uid %in% predictors, ]
  }
  
  # extract defined thresholds and "invert" (*-1) thresholds for inverse predictors
  thresholds <- deframe(select(pred_config, pred_uid, alpha))
  thresholds[pred_config$inverse_predictor] <- thresholds[pred_config$inverse_predictor] * -1
  
  # only return non-NA thresholds that are found in merged data
  thresholds <- thresholds[!is.na(thresholds)]
  
  return(thresholds)
  
}

#' Bootstrapped performance confidence intervals
#' 
#' Run bootstraps to compute confidence intervals for AUPRC, or precision or recall at a given
#' threshold performance metrics.
#' 
#' @param data Data frame containing merged data in wide format. Needs to contain columns named
#'   'name' with a unique identifier for each E-G pair and 'Regulated' identifying positives and
#'   negatives for benchmarking. All other columns are considered to be scores for predictors. Note
#'   that higher scores are assumed to rank higher. Inverse predictors (e.g. distance to TSS) need
#'   to be multiplied by -1. To convert merged data from CRISPR benchmarking to the required format,
#'   use the convertMergedForBootstrap() function.
#' @param metric Performance metric to bootstrap. Can either be 'auprc' for area under the
#'   precision-recall curve, or 'precision' or 'recall' for precision or recall at threshold
#'   at provided thresholds.
#' @param predictors Predictors (columns in data) for which performance should be bootstrapped. A
#'   simple vector or list with names of predictors for which performance should be bootstrapped.
#' @param thresholds Named vector with thresholds for all predictors (e.g. at 70% recall).
#'   Only required if metric is set to 'precision' or 'recall'.
#' @param weighted (logical) Specifies whether weighted performance should be calculated. Requires
#'   that 'data' contains a column called 'weight' containing a weight between 0 and 1 for each
#'   CRISPR E-G pair.
#' @param R Number of bootstrap replicates (default: 10000).
#' @param conf Desired confidence levels for confidence intervals (default: 0.95).
#' @param ci_type Confidence interval type. See ?boot.ci for more information.
#' @param ncpus Specifies how many CPUs should be used for bootstrapping and computing confidence
#'   intervals. If 1 not parallelization is used, if > 1 parallel computing using the specified
#'   number of CPUs will be used. Parts of parallel computing rely in BiocParallel. 
bootstrapPerformanceIntervals <- function(data, metric = c("auprc", "precision", "recall"),
                                          predictors = NULL, thresholds = NULL, weighted = FALSE,
                                          R = 10000, conf = 0.95, ncpus = 1,
                                          ci_type = c("perc", "norm", "basic", "bca"),
                                          auc_left_rectangle = FALSE) {
  
  # parse input arguments
  metric <- match.arg(metric)
  ci_type <- match.arg(ci_type)
  
  # check that thresholds are provided if precision is bootstrapped
  if (metric %in% c("precision", "recall") & is.null(thresholds)) {
    stop("Thresholds required if bootstrapping precision or recall", call. = FALSE)
  }
  
  # check that data contains 'weight' column if weighted performance is computed
  if (weighted == TRUE & !"weight" %in% colnames(data)) {
    stop("'data' needs to have weights column to calculate weighted performancs", call. = FALSE)
  }
  
  # subset data to specified predictors if passed via arguments
  if (!is.null(predictors)) {
    crispr_cols <- intersect(colnames(data), c("name", "Regulated", "weight"))
    data <- data[, c(crispr_cols, predictors)]
  }
  
  # set parallel argument for boot function
  parallel <- ifelse(ncpus > 1, yes = "multicore", no = "no")
  
  # bootstrap performance
  message("Running bootstraps...")
  bs_perf <- boot(data, statistic = calculate_performance, metric = metric, R = R,
                  parallel = parallel, ncpus = ncpus, thresholds = thresholds, weighted = weighted,
                  auc_left_rectangle = auc_left_rectangle)
  
  # set up parallel backend for computing confidence intervals if specified
  if (ncpus > 1) {
    register(MulticoreParam(workers = ncpus))
  } else {
    register(SerialParam())
  }
  
  # compute confidence intervals for all predictors
  message("Computing confidence intervals...")
  pred_indices <- seq_along(bs_perf$t0)[!is.na(bs_perf$t0)]  # indices of non-NA predictors 
  ci <- bplapply(pred_indices, FUN = boot.ci, boot.out = bs_perf, conf = conf,
                 type = ci_type)
  
  # process boot.ci output to make pretty output table
  output <- process_ci(ci, boot = bs_perf, metric = metric)
  
  return(output)
  
}

#' Bootstrapped pairwise performance comparisons
#' 
#' Run bootstraps to compute confidence intervals for delta AUPRC, or delta precision or recall at
#' threshold for specified predictor pairs. delta is simply defined as performance predictor 1 -
#' performance predictor 2.
#' 
#' @param data Data frame containing merged data in wide format. Needs to contain columns named
#'   'name' with a unique identifier for each E-G pair and 'Regulated' identifying positives and
#'   negatives for benchmarking. All other columns are considered to be scores for predictors. Note
#'   that higher scores are assumed to rank higher. Inverse predictors (e.g. distance to TSS) need
#'   to be multiplied by -1. To convert merged data from CRISPR benchmarking to the required format,
#'   use the convertMergedForBootstrap() function.
#' @param metric Performance metric to bootstrap. Can either be 'auprc' for area under the
#'   precision-recall curve, or 'precision' or 'recall' for precision or recall at threshold
#'   at provided thresholds.
#' @param comparisons List containing pairwise comparisons of predictors that should be computed
#'   (one pair per element). If 'NULL', all pairwise comparisons between all predictors will be
#'   tested.
#' @param thresholds Named vector with thresholds for all predictors (e.g. at 70% recall).
#'   Only required if metric is set to 'precision' or 'recall'.
#' @param weighted (logical) Specifies whether weighted performance should be calculated. Requires
#'   that 'data' contains a column called 'weight' containing a weight between 0 and 1 for each
#'   CRISPR E-G pair.
#' @param R Number of bootstrap replicates (default: 10000).
#' @param conf Desired confidence levels for confidence intervals (default: 0.95).
#' @param ci_type Confidence interval type. See ?boot.ci for more information.
#' @param ncpus Specifies how many CPUs should be used for bootstrapping and computing confidence
#'   intervals. If 1 not parallelization is used, if > 1 parallel computing using the specified
#'   number of CPUs will be used. Parts of parallel computing rely in BiocParallel. 
bootstrapDeltaPerformance <- function(data, metric = c("auprc", "precision", "recall"),
                                      comparisons = NULL, thresholds = NULL, weighted = FALSE,
                                      R = 10000, conf = 0.95,
                                      ci_type = c("perc", "norm", "basic", "bca"), ncpus = 1,
                                      auc_left_rectangle = auc_left_rectangle) {
  
  # parse input arguments
  metric <- match.arg(metric)
  ci_type <- match.arg(ci_type)
  
  # check that thresholds are provided if precision is bootstrapped
  if (metric == "precision" & is.null(thresholds)) {
    stop("Thresholds required if bootstrapping precision", call. = FALSE)
  }
  
  # check that data contains 'weight' column if weighted performance is computed
  if (weighted == TRUE & !"weight" %in% colnames(data)) {
    stop("'data' needs to have weights column to calculate weighted performancs", call. = FALSE)
  }
  
  # subset data to predictors in comparisons if specified, else create all pairwise comparisons
  if (!is.null(comparisons)) {
    crispr_cols <- intersect(colnames(data), c("name", "Regulated", "weight"))
    data <- data[, c(crispr_cols, unique(unlist(comparisons)))]
  } else {
    comparisons <- combn(setdiff(colnames(data), c("name", "Regulated", "weight")), m = 2,
                         simplify = FALSE)
  }
  
  # set names for all comparisons
  names(comparisons) <- vapply(comparisons, FUN = paste, FUN.VALUE = character(1), collapse = " | ")
  
  # set parallel argument for boot function
  parallel <- ifelse(ncpus > 1, yes = "multicore", no = "no")
  
  # bootstrap performance
  message("Running bootstraps...")
  bs_delta <- boot(data, statistic = calc_delta_performance, metric = metric, R = R,
                   parallel = parallel, ncpus = ncpus, thresholds = thresholds,
                   comparisons = comparisons, weighted = weighted,
                   auc_left_rectangle = auc_left_rectangle)
  
  # set up parallel backend for computing confidence intervals if specified (useful for 'bca')
  if (ncpus > 1) {
    register(MulticoreParam(workers = ncpus))
  } else {
    register(SerialParam())
  }
  
  # compute confidence intervals for all predictors
  message("Computing confidence intervals...")
  ci <- bplapply(seq_along(bs_delta$t0), FUN = boot.ci, boot.out = bs_delta, conf = conf,
                 type = ci_type)
  
  # process boot.ci output to make pretty output table
  output <- process_ci(ci, boot = bs_delta, metric = paste0("delta_", metric))
  
  # compute p-values under the null hypothesis that delta is 0
  message("Computing p-values...")
  pvalues <- compute_pvalues(bs_delta, type = ci_type, theta_null = 0, pval_precision = NULL)
  output$pvalue <- pvalues[output$id]
  
  return(output)
  
}

#' Bootstrapped performance differences between datasets
#' 
#' Run bootstraps to compute confidence intervals for delta AUPRC, or delta precision or recall at
#' threshold between two benchmarking datasets for specified predictors. delta is simply defined as
#' performance predictor 1 - performance predictor 2.
#' 
#' @param data1,data2 Data frames containing merged data for the two benchmarking datasets in wide
#'   format. Need to contain columns named 'name' with a unique identifier for each E-G pair and
#'   'Regulated' identifying positives and negatives for benchmarking. All other columns are
#'   considered to be scores for predictors. Note that higher scores are assumed to rank higher.
#'   Inverse predictors (e.g. distance to TSS) need to be multiplied by -1. To convert merged data
#'   from CRISPR benchmarking to the required format, use the convertMergedForBootstrap() function.
#' @param metric Performance metric to bootstrap. Can either be 'auprc' for area under the
#'   precision-recall curve, or 'precision' or 'recall' for precision or recall at threshold
#'   at provided thresholds.
#' @param predictors Predictors (columns in data) for which performance differences should be
#'  computed. A simple vector or list with names of predictors to include. If not specified, the
#'  intersect between predictors in data1 and data2 will be used.
#' @param thresholds Named vector with thresholds for all predictors (e.g. at 70% recall).
#'   Only required if metric is set to 'precision' or 'recall'.
#' @param weighted (logical) Specifies whether weighted performance should be calculated. Requires
#'   that 'data' contains a column called 'weight' containing a weight between 0 and 1 for each
#'   CRISPR E-G pair.
#' @param R Number of bootstrap replicates (default: 10000).
#' @param conf Desired confidence levels for confidence intervals (default: 0.95).
#' @param ci_type Confidence interval type. See ?boot.ci for more information.
#' @param ncpus Specifies how many CPUs should be used for bootstrapping and computing confidence
#'   intervals. If 1 not parallelization is used, if > 1 parallel computing using the specified
#'   number of CPUs will be used. Parts of parallel computing rely in BiocParallel.
bootstrapDeltaPerformanceDatasets <- function(data1, data2,
                                              metric = c("auprc", "precision", "recall"),
                                              predictors = NULL, thresholds = NULL,
                                              weighted = FALSE, R = 10000, conf = 0.95, ncpus = 1,
                                              ci_type = c("perc", "norm", "basic", "bca"),
                                              auc_left_rectangle = FALSE) {
  
  # parse input arguments
  metric <- match.arg(metric)
  ci_type <- match.arg(ci_type)
  
  # check that thresholds are provided if precision is bootstrapped
  if (metric == "precision" & is.null(thresholds)) {
    stop("Thresholds required if bootstrapping precision", call. = FALSE)
  }
  
  # check that data contains 'weight' column if weighted performance is computed
  if (weighted == TRUE & !"weight" %in% colnames(data1)) {
    stop("'data1' needs to have weights column to calculate weighted performancs", call. = FALSE)
  }
  if (weighted == TRUE & !"weight" %in% colnames(data2)) {
    stop("'data2' needs to have weights column to calculate weighted performancs", call. = FALSE)
  }
  
  # get predictors shared across datasets unless predictors to analyze are specified
  if (is.null(predictors)) {
    predictors <- setdiff(intersect(colnames(data1), colnames(data2)),
                          c("name", "Regulated", "weight"))
  }
  
  # subset data to predictors to analyze
  crispr_cols1 <- intersect(colnames(data1), c("name", "Regulated", "weight"))
  crispr_cols2 <- intersect(colnames(data2), c("name", "Regulated", "weight"))
  data1 <- data1[, c(intersect(crispr_cols1, crispr_cols2), predictors)]
  data2 <- data2[, c(intersect(crispr_cols1, crispr_cols2), predictors)]
  
  # combine both datasets into one table and add column specifying dataset for each E-G pair
  data <- bind_rows(data1, data2, .id = "dataset")
  
  # set parallel argument for boot function
  parallel <- ifelse(ncpus > 1, yes = "multicore", no = "no")
  
  # bootstrap performance
  message("Running bootstraps...")
  bs_delta <- boot(data, statistic = calc_delta_performance_datasets, metric = metric, R = R,
                   strata = data$dataset, parallel = parallel, ncpus = ncpus,
                   thresholds = thresholds, weighted = weighted,
                   auc_left_rectangle = auc_left_rectangle)
  
  # set up parallel backend for computing confidence intervals if specified (useful for 'bca')
  if (ncpus > 1) {
    register(MulticoreParam(workers = ncpus))
  } else {
    register(SerialParam())
  }
  
  # compute confidence intervals for all predictors
  message("Computing confidence intervals...")
  ci <- bplapply(seq_along(bs_delta$t0), FUN = boot.ci, boot.out = bs_delta, conf = conf,
                 type = ci_type)
  
  # process boot.ci output to make pretty output table
  output <- process_ci(ci, boot = bs_delta, metric = paste0("delta_", metric))
  
  # compute p-values under the null hypothesis that delta is 0
  message("Computing p-values...")
  pvalues <- compute_pvalues(bs_delta, type = ci_type, theta_null = 0, pval_precision = NULL)
  output$pvalue <- pvalues[output$id]
  
  return(output) 
  
}

#' Plot bootstrapped performance / delta performance
#' 
#' @param results Data frame containing bootstrapped performance metrics or delta performance.
#' @param title (optional) Main title for plot.
plotBootstrappedIntervals <- function(results, title = NULL) {
  
  # set id to factor ordered full metric and add column specifiying whether CI is different from 0
  results <- results %>% 
    mutate(id = fct_reorder(id, .x = full)) %>% 
    mutate(diff_zero = if_else(lower <= 0 & upper >= 0, true = FALSE, false = TRUE))
  
  # create default main title
  default_title <- switch(unique(results$metric), "auprc" = "AUPRC", "precision" = "Precision",
                          "delta_auprc" = "Delta AUPRC", "delta_precision" = "Delta Precision")
  
  # create axis title
  if (unique(results$metric) %in% c("delta_auprc", "delta_precision")) {
    axis_title <- paste(default_title, "[first predictor - second predictor]")
  } else {
    axis_title <- default_title
  }
  
  # create title for color/fill legend
  color_title <- paste0(unique(results$conf) * 100, "% interval\ndifferent from 0")
  
  # use default main title unless manually specified otherwise
  if (is.null(title)) title <- default_title
  
  # plot mean delta, 95% intervals and range for all comparisons
  ggplot(results, aes(x = id, y = full)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = diff_zero), linewidth = 1.5, width = 0) +
    geom_errorbar(aes(ymin = min, ymax = max, color = diff_zero), linewidth = 0.5, width = 0) +
    geom_point(aes(color = diff_zero, fill = diff_zero), shape = 23, size = 2) +
    labs(title = title, y = axis_title, color = color_title, fill = color_title) +
    coord_flip() +
    scale_color_manual(values = c("FALSE" = "gray55", "TRUE" = "black")) +
    scale_fill_manual(values = c("FALSE" = "gray55", "TRUE" = "red")) +
    theme_bw() +
    theme(axis.title.y = element_blank())
  
}

## FUNCTIONS TO COMPUTE PERFORMANCE AND DELTA PERFORMANCE ==========================================

# calculate performance AUPRC, or precision or recall at threshold
calculate_performance <- function(data, indices, metric, thresholds, weighted,
                                  auc_left_rectangle = FALSE) {
  
  # select bootstrap sample
  data <- data[indices, ]
  
  # get all predictors in input data
  preds <- setdiff(colnames(data), c("name", "Regulated", "weight", "dataset"))
  
  # if thresholds are provided, get thresholds for these predictors else create NA thresholds
  if (!is.null(thresholds)) {
    thresholds <- thresholds[preds]
  } else {
    thresholds <- rep_len(NA_real_, length(preds))
  }
  
  # calculate performance for all predictors
  performance <- mapply(FUN = calculate_performance_one_pred, pred = preds, threshold = thresholds,
                        MoreArgs = list(data = data, metric = metric, weighted = weighted,
                                        auc_left_rectangle = auc_left_rectangle),
                        SIMPLIFY = TRUE)
  
  return(performance)
  
}

# function to calculate delta AUPRC, or precision or recall at threshold between pairwise
# predictor combinations
calc_delta_performance <- function(data, indices, metric, thresholds, comparisons, weighted,
                                   auc_left_rectangle = FALSE) {
  
  # calculate bootstrapped performance
  perf <- calculate_performance(data, indices = indices, metric = metric, thresholds = thresholds,
                                weighted = weighted,
                                auc_left_rectangle = auc_left_rectangle)
  
  # calculate delta performance for all specified comparisons
  delta_perf <- vapply(comparisons, FUN = function(comp, perf) {
    perf[[comp[[1]]]] - perf[[comp[[2]]]]
  }, perf = perf, FUN.VALUE = numeric(1))
  
  return(delta_perf)
  
}

# function to calculate delta AUPRC, or precision or recall at threshold between 2 datasets
calc_delta_performance_datasets <- function(data, indices, metric, thresholds, weighted,
                                            auc_left_rectangle = FALSE) {
  
  # select bootstrap sample
  data <- data[indices, ]
  
  # calculate bootstrapped performance for both stratifications
  data1 <- data[data$dataset == "1", ]
  data2 <- data[data$dataset == "2", ]
  perf1 <- calculate_performance(data1, indices = seq_len(nrow(data1)), metric = metric,
                                 thresholds = thresholds, weighted = weighted,
                                 auc_left_rectangle = auc_left_rectangle)
  perf2 <- calculate_performance(data2, indices = seq_len(nrow(data2)), metric = metric,
                                 thresholds = thresholds, weighted = weighted,
                                 auc_left_rectangle = auc_left_rectangle)
  
  # calculate delta performance for all predictors
  delta_perf <- perf1 - perf2
    
  return(delta_perf)
  
}


## HELPER FUNCTIONS ================================================================================

# calculate performance auprc, or precision or recall at threshold for one predictor
calculate_performance_one_pred <- function(data, pred, threshold, metric, weighted,
                                           auc_left_rectangle = FALSE) {
  
  # return NA if 'Regulated' column does not contain at least one positive and negative
  if (length(unique(data$Regulated)) != 2) {
    warning("Not both positives and negatives ('Regulated') in bootstrap sample. Returning 'NA'.",
            call. = FALSE)
    return(NA_real_)
  }
  
  # compute precision-recall curve
  if (weighted == TRUE) {
    pr <- pr_curve(data, !!sym(pred), truth = Regulated, event_level = "second",
                   case_weights = weight)
  } else {
    pr <- pr_curve(data, !!sym(pred), truth = Regulated, event_level = "second")
  }
  
  # remove top an bottom rows to make AUPRC calculation consistent with CRISPR benchmarking pipeline
  pr <- head(pr, -1)[-1,]
  
  # rename and reorder columns
  pr <- data.frame(
    alpha = pr$.threshold,
    precision = pr$precision,
    recall = pr$recall
  )
  
  # calculate AUPRC, or precision or recall at threshold performance
  if (metric %in% c("precision", "recall")) {
    performance <- calculate_performance_at_threshold(pr, threshold = threshold, metric = metric)
  } else if (metric == "auprc") {
    performance <- calculate_auprc(pr, auc_left_rectangle = auc_left_rectangle)
  } else {
    stop("Invalid 'metric' argument", call. = FALSE)
  }
  
  return(performance)
  
}

# calculate area-under-the-precision-recall-curve (AUPRC)
calculate_auprc <- function(pr, auc_left_rectangle = FALSE) {
  
  # the head() calls here remove the last element of the vector. 
  # The point is that performance objects produced by ROCR always include a Recall = 100% point even
  # if the predictor cannot achieve a recall of 100%. This results in a straight line ending at
  # (1,0) on the PR curve. This should not be included in the performance computation.
  pr <- head(pr, -1)
  
  # compute auprc
  auprc <- compute_auc(x_vals = pr$recall, y_vals = pr$precision,
                       auc_left_rectangle = auc_left_rectangle)
  
  return(auprc)
  
}

# try to compute area under the curve
compute_auc <- function(x_vals, y_vals, 
                        auc_left_rectangle = FALSE) {
  good.idx <- which(!is.na(x_vals) & !is.na(y_vals))
  if (length(good.idx) > 0) {
    auc <- trapz(x_vals[good.idx], y_vals[good.idx])
    # Add the area of the rectangle formed by the left end 
    # of the curve and the origin.
    if(auc_left_rectangle == "True"){
      auc <- auc + x_vals[good.idx[1]] * y_vals[good.idx[1]]
    }
  } else {
    auc <- NA_real_
  }
  return(auc)
}

# calculate precision at a given threshold
calculate_performance_at_threshold <- function(pr, threshold, metric) {
  
  # get index of highest alpha value that is larger or equal to alpha_cutoff
  idx <- sum(pr$alpha >= threshold)
  
  # get precision at this alpha value
  perf_at_threshold <- pr[[metric]][[idx]]
  
  return(perf_at_threshold)
  
}

# function to extract the full value and upper and lower CI boundaries for a given bootci object
extract_ci <- function(bootci) {
  
  # get full data metric
  full <- bootci$t0
  
  # get CI data for given CI type
  ci <- bootci[[4]]
  
  # assemble output data.frame
  output <- data.frame(id = names(full), full = full, conf = ci[[1]], lower = ci[[length(ci) - 1]],
                       upper = ci[[length(ci)]], row.names = NULL)
  
  return(output)
  
}

# function to extract full range of
extract_range <- function(boot) {
  
  # calculate range of bootstrapped metric
  range <- apply(boot$t, MARGIN = 2, FUN = range, na.rm = TRUE)
  rownames(range) <- c("min", "max")
  colnames(range) <- names(boot$t0)
  
  # transpose and make names a column
  range <- as.data.frame(t(range))
  range$id <- rownames(range)
  rownames(range) <- NULL
  
  return(range)
  
}

# process outpur from boot.ci function and calculate absolute range of values
process_ci <- function(ci, boot, metric) {
  
  # extract relevant information for confidence intervals and create output data.frame
  ci <- lapply(ci, FUN = extract_ci)
  ci <- rbindlist(ci)
  
  # extract full range of bootstrapped metric
  range <- extract_range(boot)
  
  # add range to confidence intervals to create output
  output <- merge(ci, range, by = "id")
  
  # add metric and rearrange columns
  output <- mutate(output, metric = metric, .after = 1)
  
  return(output)
  
}

# compute p-values from bootstrapping results
compute_pvalues <- function(boot, type, theta_null, pval_precision) {
  
  # get indices for different bootstraps
  indices <- structure(seq_along(boot$t0), names = names(boot$t0))
  
  # compute p-values from boot object
  pvals <- bplapply(indices, FUN = boot.pval, boot_res = boot, type = type, theta_null = theta_null,
                    pval_precision = pval_precision)
  
  # convert from list ot simple vector
  pvals <- unlist(pvals)
  
  return(pvals)
  
}

# function to compute p-values from boot objects. taken from the boot.pval package, which is not on
# conda yet: https://github.com/mthulin/boot.pval/blob/main/R/boot.pval.R
boot.pval <- function(boot_res,
                      type = "perc",
                      theta_null = 0,
                      pval_precision = NULL,
                      ...) {
  
  if(is.null(pval_precision)) { pval_precision = 1/boot_res$R }
  
  # create a sequence of alphas:
  alpha_seq <- seq(1e-16, 1-1e-16, pval_precision)
  
  # compute the 1-alpha confidence intervals, and extract their bounds:
  ci <- suppressWarnings(boot::boot.ci(boot_res,
                                       conf = 1- alpha_seq,
                                       type = type,
                                       ...))
  
  bounds <- switch(type,
                   norm = ci$normal[,2:3],
                   basic = ci$basic[,4:5],
                   stud = ci$student[,4:5],
                   perc = ci$percent[,4:5],
                   bca = ci$bca[,4:5])
  
  # find the smallest alpha such that theta_null is not contained in the 1-alpha
  # confidence interval:
  alpha <- alpha_seq[which.min(theta_null >= bounds[,1] & theta_null <= bounds[,2])]
  
  # return the p-value:
  return(alpha)
  
}
