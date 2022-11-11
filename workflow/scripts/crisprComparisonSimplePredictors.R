## A couple of simple baseline predictors

# function to compute specified baseline predictors
computeBaselinePreds <- function(expt, tss_annot = NULL, gene_annot = NULL, expressed_genes = NULL,
                                 preds = c("distToTSS", "distToGene", "nearestTSS", "nearestGene",
                                           "nearestExprTSS", "nearestExprGene")) {
  
  # parse preds argument
  preds <- match.arg(preds, several.ok = TRUE)
  
  # get categories of specified predictors to check input data
  tss_preds  <- intersect(preds, c("distToTSS","nearestTSS", "nearestExprTSS"))
  gene_preds <- intersect(preds, c("distToGene", "nearestGene", "nearestExprGene"))
  expr_preds <- intersect(preds, c("nearestExprTSS", "nearestExprGene"))
  
  # abort if required input data is not provided for all specified predictors
  if (length(tss_preds) > 0 & is.null(tss_annot)) {
    stop("tss_annot required to compute: ", paste(tss_preds, collapse = ", "), call. = FALSE)
  }
  if (length(gene_preds) > 0 & is.null(gene_annot)) {
    stop("gene_annot required to compute: ", paste(gene_preds, collapse = ", "), call. = FALSE)
  }
  if (length(expr_preds) > 0 & is.null(expressed_genes)) {
    stop("expr_genes required to compute: ", paste(expr_preds, collapse = ", "), call. = FALSE)
  }
  
  # extract expressed gene names if expressed_genes is provided
  if (!is.null(expressed_genes)) {
    expr_genes <- expressed_genes[expressed_genes$expressed == TRUE, ][["gene"]]
  }
  
  # compute baseline predictors
  baseline_preds <- lapply(preds, FUN = compute_baseline_predictor, expt = expt,
                           tss_annot = tss_annot, gene_annot = gene_annot, expr_genes = expr_genes)
  
  # combine into one table
  output <- rbindlist(baseline_preds)
  
  return(output)
  
}

## HELPER FUNCTIONS ================================================================================

# function to compute one baseline predictor
compute_baseline_predictor <- function(pred, expt, tss_annot, gene_annot, expr_genes) {
  
  # compute baseline predictors
  output <- switch(
    pred,
    "distToTSS" = computeDistToGene(expt, annot = tss_annot, name = pred, fix_annot = "center"),
    "distToGene" = computeDistToGene(expt, annot = gene_annot, name = pred, fix_annot = "none"),
    "nearestTSS" = nearestFeaturePred(expt, features = tss_annot, name = pred),
    "nearestGene" = nearestFeaturePred(expt, features = gene_annot, name = pred),
    "nearestExprTSS" = nearestFeaturePred(
      expt, features = tss_annot[tss_annot$gene %in% expr_genes, ], name = pred),
    "nearestExprGene" = nearestFeaturePred(
      expt, features = gene_annot[gene_annot$gene %in% expr_genes, ], name = pred)
  )

  return(output)
  
}

# compute baseline distance to any gene annotation (TSS or gene)
computeDistToGene <- function(expt, annot, name, fix_annot = c("none", "center", "start", "end")) {
  
  # parse fix argument
  fix_annot <- match.arg(fix_annot, choices = c("none", "center", "start", "end"))
  
  # assume first three column of annot are chr, start, end and convert annot to GRanges object
  colnames(annot)[1:3] <- c("chr", "start", "end")
  annot <- makeGRangesFromDataFrame(annot, keep.extra.columns = TRUE)
  names(annot) <- annot$gene
  
  # convert gene annotations to 1bp coordinates if specified by fix_annot (e.g. for distance to TSS)
  if (fix_annot != "none") {
    annot <- resize(annot, width = 1, fix = fix_annot)
  }
  
  # GRanges object for E-G pairs in expt using centers of candidate enhancers as coordinates
  eg_pairs <- makeGRangesFromDataFrame(expt, seqnames.field = "chrom", start.field = "chromStart",
                                       end = "chromEnd", keep.extra.columns = TRUE)
  eg_pairs <- resize(eg_pairs, width = 1, fix = "center")
  
  # combine into a paired GRanges object containing enhancer and TSS/gene annotations for each pair
  eg_pairs <- Pairs(first = eg_pairs, second = annot[expt$measuredGeneSymbol])
  
  # compute distance between enhancers and provided TSS / gene annotations
  distance_pred <- distance(eg_pairs)
  
  # create output table
  output <- data.table(
    expt,
    PredictionCellType = NA_character_,
    pred_id = "baseline",
    pred_col = name,
    pred_value = distance_pred,
    Prediction = 1
  )
  
  return(output)
  
}

# nearest genomic feature (e.g. TSS or gene body) to the candidate enhancer. features needs to be a
# data frame in bed style format, i.e chr, start, end, name, etc.
nearestFeaturePred <- function(expt, features, name) {
  
  # create GRanges object containing TSS coordinates
  features <- features[, 1:4]
  colnames(features) <- c("chr", "start", "stop", "name")
  features <- makeGRangesFromDataFrame(features, keep.extra.columns = TRUE)
  
  # get closest TSS to every CRE in expt
  expt <- nearest_genomic_feature(expt, features = features)
  
  # create predictor whether the gene of each pair is the closest TSS
  output <- data.table(
    expt[, -"nearest_feature"],
    PredictionCellType = NA_character_,
    pred_id = "baseline",
    pred_col = name,
    pred_value = as.numeric(expt$nearest_feature == expt$measuredGeneSymbol),
    Prediction = 1
  )
  
}

# find the closest genomic feature for every enhancer in experimental data
nearest_genomic_feature <- function(expt, features, ignore.strand = TRUE) {
  
  # add uniaue identifier for every E-G pair in expt
  expt$uid <- paste0("eg_pair_", seq_len(nrow(expt)))
  
  # create genomic ranges for CREs
  cres <- makeGRangesFromDataFrame(expt, seqnames.field = "chrom", start.field = "chromStart",
                                   end.field = "chromEnd", keep.extra.columns = TRUE)
  
  # find nearest feature for every CRE
  nearest_feat <- nearest(cres, subject = features, ignore.strand = ignore.strand)
  
  # get name of nearest feature for each CRE
  nearest_feat_names <- data.frame(uid = cres$uid, nearest_feature = features[nearest_feat]$name)
  
  # add this to expt data
  output <- merge(expt, nearest_feat_names, by = "uid", all.x = TRUE)
  
  return(output[, -"uid"])
  
}

## DEPRECATED ======================================================================================

# compute baseline 'distance to TSS' predictor
computeDistToTSS <- function(expt) {
  
  # compute distance to TSS
  dist_to_tss <- with(
    expt,
    abs((chromStart + chromEnd) / 2 - (startTSS + endTSS) / 2)
  )
  
  # create output table
  output <- data.table(
    expt,
    PredictionCellType = NA_character_,
    pred_id = "baseline",
    pred_col = "distToTSS",
    pred_value = dist_to_tss,
    Prediction = 1
  )
  
  return(output)
  
}
