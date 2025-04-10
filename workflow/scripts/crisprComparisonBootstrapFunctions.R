## Functions to compute bootstrapped performance metrics (AUPRC and precision @ threshold) and
## significance for pairwise differences in performance between predictors

library(boot)
#library(boot.pval)
library(ROCR)
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
#' @param rm_boolean (logical) Should boolean predictors be filtered out? (default: TRUE)
convertMergedForBootstrap <- function(merged, pred_config, pos_col = "Regulated",
                                      rm_boolean = TRUE) {
  
  # get inverse predictors
  inverse_preds <- pred_config[pred_config$inverse_predictor == TRUE, ][["pred_uid"]]
  
  # extract relevant columns from merged data
  merged <- select(merged, name, Regulated = all_of(pos_col), pred_uid, pred_value)
  
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
  
  # convert to wide format to create output
  output <- pivot_wider(merged, names_from = pred_uid, values_from = pred_value)
  
  return(output)
  
}

#' Bootstrapped performance confidence intervals
#' 
#' Run bootstraps to compute confidence intervals for AUPRC or precision at a given threshold
#' performance metrics.
#' 
#' @param data Data frame containing merged data in wide format. Needs to contain columns named
#'   'name' with a unique identifier for each E-G pair and 'Regulated' identifying positives and
#'   negatives for benchmarking. All other columns are considered to be scores for predictors. Note
#'   that higher scores are assumed to rank higher. Inverse predictors (e.g. distance to TSS) need
#'   to be multiplied by -1. To convert merged data from CRISPR benchmarking to the required format,
#'   use the convertMergedForBootstrap() function.
#' @param metric Performance metric to bootstrap. Can either be "auprc" or "precision" for precision
#'   at provided thresholds.
#' @param predictors Predictors (columns in data) for which performance should be bootstrapped. A
#'   simple vector or list with names of predictors for which performance should be bootstrapped.
#' @param thresholds Named vector with thresholds for all predictors (e.g. at 70% recall).
#'   Only required if metric is set to 'precision'.
#' @param R Number of bootstrap replicates (default: 10000).
#' @param conf Desired confidence levels for confidence intervals (default: 0.95).
#' @param ci_type Confidence interval type. See ?boot.ci for more information.
#' @param ncpus Specifies how many CPUs should be used for bootstrapping and computing confidence
#'   intervals. If 1 not parallelization is used, if > 1 parallel computing using the specified
#'   number of CPUs will be used. Parts of parallel computing rely in BiocParallel. 
bootstrapPerformanceIntervals <- function(data, metric = c("auprc", "precision"), predictors = NULL,
                                          thresholds = NULL, R = 10000, conf = 0.95,
                                          ci_type = c("perc", "norm", "basic", "bca"), ncpus = 1) {
  
  # parse input arguments
  metric <- match.arg(metric)
  ci_type <- match.arg(ci_type)
  
  # check that thresholds are provided if precision is bootstrapped
  if (metric == "precision" & is.null(thresholds)) {
    stop("Thresholds required if bootstrapping precision", call. = FALSE)
  }
  
  # get function to compute specified metric
  metric_fun <- switch(metric, "auprc" = calculate_auprc, "precision" = calculate_precision)
  
  # subset data to specified predictors if passed via arguments
  if (!is.null(predictors)) {
    data <- data[, c("name", "Regulated", predictors)]
  }
  
  # set parallel argument for boot function
  parallel <- ifelse(ncpus > 1, yes = "multicore", no = "no")
  
  # bootstrap performance
  message("Running bootstraps...")
  bs_perf <- boot(data, statistic = metric_fun, R = R, parallel = parallel, ncpus = ncpus,
                  thresholds = thresholds)
  
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
#' Run bootstraps to compute confidence intervals for delta AUPRC or delta precision (at threshold)
#' for specified predictor pairs. delta is simply defined as predictor 1 - predictor 2.
#' 
#' @param data Data frame containing merged data in wide format. Needs to contain columns named
#'   'name' with a unique identifier for each E-G pair and 'Regulated' identifying positives and
#'   negatives for benchmarking. All other columns are considered to be scores for predictors. Note
#'   that higher scores are assumed to rank higher. Inverse predictors (e.g. distance to TSS) need
#'   to be multiplied by -1. To convert merged data from CRISPR benchmarking to the required format,
#'   use the convertMergedForBootstrap() function.
#' @param metric Performance metric to bootstrap. Can either be "auprc" or "precision" for precision
#'   at provided thresholds.
#' @param comparisons List containing pairwise comparisons of predictors that should be computed
#'   (one pair per element). If 'NULL', all pairwise comparisons between all predictors will be
#'   tested.
#' @param thresholds Named vector with thresholds for all predictors (e.g. at 70% recall).
#'   Only required if metric is set to 'precision'.
#' @param R Number of bootstrap replicates (default: 10000).
#' @param conf Desired confidence levels for confidence intervals (default: 0.95).
#' @param ci_type Confidence interval type. See ?boot.ci for more information.
#' @param ncpus Specifies how many CPUs should be used for bootstrapping and computing confidence
#'   intervals. If 1 not parallelization is used, if > 1 parallel computing using the specified
#'   number of CPUs will be used. Parts of parallel computing rely in BiocParallel. 
bootstrapDeltaPerformance <- function(data, metric = c("auprc", "precision"), comparisons = NULL,
                                      thresholds = NULL, R = 10000, conf = 0.95,
                                      ci_type = c("perc", "norm", "basic", "bca"), ncpus = 1) {
  
  # parse input arguments
  metric <- match.arg(metric)
  ci_type <- match.arg(ci_type)
  
  # check that thresholds are provided if precision is bootstrapped
  if (metric == "precision" & is.null(thresholds)) {
    stop("Thresholds required if bootstrapping precision", call. = FALSE)
  }
  
  # get function to compute specified metric
  metric_fun <- switch(metric, "auprc" = calc_delta_auprc, "precision" = calc_delta_precision)
  
  # subset data to predictors in comparisons if specified, else create all pairwise comparisons
  if (!is.null(comparisons)) {
    data <- data[, c("name", "Regulated", unique(unlist(comparisons)))]
  } else {
    comparisons <- combn(setdiff(colnames(data), c("name", "Regulated")), m = 2, simplify = FALSE)
  }
  
  # set names for all comparisons
  names(comparisons) <- vapply(comparisons, FUN = paste, FUN.VALUE = character(1), collapse = " | ")
  
  # set parallel argument for boot function
  parallel <- ifelse(ncpus > 1, yes = "multicore", no = "no")
  
  # bootstrap performance
  message("Running bootstraps...")
  bs_delta <- boot(data, statistic = metric_fun, R = R, parallel = parallel, ncpus = ncpus,
                   thresholds = thresholds, comparisons = comparisons)
  
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
    geom_errorbar(aes(ymin = lower, ymax = upper, color = diff_zero), size = 1.5, width = 0) +
    geom_errorbar(aes(ymin = min, ymax = max, color = diff_zero), size = 0.5, width = 0) +
    geom_point(aes(color = diff_zero, fill = diff_zero), shape = 23, size = 2) +
    labs(title = title, y = axis_title, color = color_title, fill = color_title) +
    coord_flip() +
    scale_color_manual(values = c("FALSE" = "gray55", "TRUE" = "black")) +
    scale_fill_manual(values = c("FALSE" = "gray55", "TRUE" = "red")) +
    theme_bw() +
    theme(axis.title.y = element_blank())
  
}

## FUNCTIONS TO COMPUTE PERFORMANCE AND DELTA PERFORMANCE ==========================================

# function to calculate AUPRC for all predictors (thresholds is ignored)
calculate_auprc <- function(data, indices, thresholds) {
  
  # select bootstrap sample
  data <- data[indices, ]
  
  # get all predictors in input data
  preds <- setdiff(colnames(data), c("name", "Regulated"))
  
  # calculate performance for all predictors
  auprc <- vapply(preds, FUN = calculate_auprc_one_pred, data = data, FUN.VALUE = numeric(1))
  
  return(auprc)
  
}

# function to calculate precision at threshold for all predictors
calculate_precision <- function(data, indices, thresholds) {
  
  # select bootstrap sample
  data <- data[indices, ]
  
  # get all predictors in input data
  preds <- setdiff(colnames(data), c("name", "Regulated"))
  
  # get thresholds for these predictors
  thresholds <- thresholds[preds]
  
  # calculate performance for all predictors
  precision <- mapply(FUN = calculate_precision_one_pred, pred = preds, threshold = thresholds,
                      MoreArgs = list(data = data), SIMPLIFY = TRUE)
  
  return(precision)
  
}

# function to calculate delta AUPRC between all predictor pairwise combinations
calc_delta_auprc <- function(data, indices, thresholds, comparisons) {
  
  # calculate bootstrapped auprc
  auprc <- calculate_auprc(data, indices = indices, thresholds = thresholds)
  
  # calculate delta auprc for all specified comparisons
  delta_auprc <- vapply(comparisons, FUN = function(comp, perf) {
    perf[[comp[[1]]]] - perf[[comp[[2]]]]
  }, perf = auprc, FUN.VALUE = numeric(1))
  
  return(delta_auprc)
  
}

# function to calculate delta precision between all predictor pairwise combinations
calc_delta_precision <- function(data, indices, thresholds, comparisons) {
  
  # calculate bootstrapped precision  
  precision <- calculate_precision(data, indices = indices, thresholds = thresholds)
  
  # calculate delta precision for all specified comparisons
  delta_precision <- vapply(comparisons, FUN = function(comp, perf) {
    perf[[comp[[1]]]] - perf[[comp[[2]]]]
  }, perf = precision, FUN.VALUE = numeric(1))
  
  return(delta_precision)
  
}

## HELPER FUNCTIONS ================================================================================

# calculate AUPRC for a given predictors
calculate_auprc_one_pred <- function(data, pred) {
  
  # return NA if 'Regulated' column does not contain at least one positive and negative
  if (length(unique(data$Regulated)) != 2) {
    warning("Not both positives and negatives ('Regulated') in bootstrap sample. Returning 'NA'.",
            call. = FALSE)
    return(NA_real_)
  }
  
  # compute precision-recall curve
  pr <- performance(prediction(data[[pred]], data$Regulated), measure = "prec", x.measure = "rec")
  
  # convert to data.frame
  pr <- data.frame(
    alpha = pr@alpha.values[[1]],
    precision = pr@y.values[[1]],
    recall = pr@x.values[[1]]
  )
  
  # the head() calls here remove the last element of the vector. 
  # The point is that performance objects produced by ROCR always include a Recall = 100% point even
  # if the predictor cannot achieve a recall of 100%. This results in a straight line ending at
  # (1,0) on the PR curve. This should not be included in the performance computation.
  pr <- head(pr, -1)
  
  # compute AUPRC
  auprc <- calculate_auc(x_vals = pr$recall, y_vals = pr$precision)
  
  return(auprc)
  
}

# calculate precision at threshold for all predictors
calculate_precision_one_pred <- function(data, pred, threshold) {
  
  # return NA if 'Regulated' column does not contain at least one positive and negative
  if (length(unique(data$Regulated)) != 2) {
    warning("Not both positives and negatives ('Regulated') in bootstrap sample. Returning 'NA'.",
            call. = FALSE)
    return(NA_real_)
  }
  
  # compute precision-recall curve
  pr <- performance(prediction(data[[pred]], data$Regulated), measure = "prec", x.measure = "rec")
  
  # convert to data.frame
  pr <- data.frame(
    alpha = pr@alpha.values[[1]],
    precision = pr@y.values[[1]],
    recall = pr@x.values[[1]]
  )
  
  # calculate precision at threshold
  perc_at_threshold <- calculate_precision_at_threshold(pr, threshold = threshold)
  
}

# try to compute AUC
calculate_auc <- function(x_vals, y_vals) {
  good.idx <- which(!is.na(x_vals) & !is.na(y_vals))
  if (length(good.idx) > 0) {
    auc <- trapz(x_vals[good.idx], y_vals[good.idx])
    # Add the area of the rectangle formed by the left end 
    # of the curve and the origin.
    auc <- auc + x_vals[good.idx[1]] * y_vals[good.idx[1]]
  } else {
    auc <- NA_real_
  }
  return(auc)
}

# calculate precision at a given threshold
calculate_precision_at_threshold <- function(pr, threshold) {
  
  # get index of highest alpha value that is larger or equal to alpha_cutoff
  idx <- sum(pr$alpha >= threshold)
  
  # get precision at this alpha value
  prec_at_threshold <- pr$precision[[idx]]
  
  return(prec_at_threshold)
  
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
