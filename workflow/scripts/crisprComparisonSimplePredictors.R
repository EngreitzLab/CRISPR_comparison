## A couple of simple baseline predictors

# function to compute specified baseline predictors
computeBaselinePreds <- function(expt, tss_annot = NULL, gene_annot = NULL, expressed_genes = NULL,
                                 preds = c("distToTSS", "distToGene", "nearestTSS", "nearestGene",
                                           "within100kbTSS", "within100kbGene", "nearestExprTSS",
                                           "nearestExprGene", "within100kbExprTSS",
                                           "within100kbExprGene")) {
  
  # parse preds argument
  preds <- match.arg(preds, several.ok = TRUE)
  
  # get categories of specified predictors to check input data
  tss_preds  <- intersect(preds, c("distToTSS","nearestTSS", "nearestExprTSS", "within100kbTSS",
                                   "within100kbExprTSS"))
  gene_preds <- intersect(preds, c("distToGene", "nearestGene", "nearestExprGene",
                                   "within100kbExprGene"))
  expr_preds <- intersect(preds, c("nearestExprTSS", "nearestExprGene", "within100kbExprTSS",
                                   "within100kbExprGene"))
  
  # abort if required input data is not provided for all specified predictors
  if (length(tss_preds) > 0 & is.null(tss_annot)) {
    stop("tss_annot required to compute: ", paste(tss_preds, collapse = ", "), call. = FALSE)
  }
  if (length(gene_preds) > 0 & is.null(gene_annot)) {
    stop("gene_annot required to compute: ", paste(gene_preds, collapse = ", "), call. = FALSE)
  }
  if (length(expr_preds) > 0 & is.null(expressed_genes)) {
    stop("expressed_genes required to compute: ", paste(expr_preds, collapse = ", "), call. = FALSE)
  }

  # compute baseline predictors for each cell type in experimental data
  cell_types <- unique(expt$CellType)
  baseline_preds <- lapply(cell_types, FUN = compute_baseline_predictors_cell_type, preds = preds, 
                           expt = expt, tss_annot = tss_annot, gene_annot = gene_annot,
                           expr_genes = expressed_genes)
  
  # combine into one table and reformat for output
  baseline_preds <- rbindlist(baseline_preds)
  colnames(baseline_preds)[colnames(baseline_preds) == "CellType"] <- "ExperimentCellType"
  
  return(baseline_preds)
  
}

## HELPER FUNCTIONS ================================================================================

# compute all baseline predictors for one cell type (required for 'expressed genes' predictors)
compute_baseline_predictors_cell_type <- function(x, preds, expt, tss_annot, gene_annot,
                                                  expr_genes) {
  
  # extract experimental data for given cell type
  expt <- expt[CellType == x, ]
  
  # extract expressed gene names if expr_genes is provided
  if (!is.null(expr_genes)) {
    expr_genes <- expr_genes[cell_type == x & expressed == TRUE, ][["gene"]]
  }
  
  # compute baseline predictors
  baseline_preds <- lapply(preds, FUN = compute_baseline_predictor, expt = expt,
                           tss_annot = tss_annot, gene_annot = gene_annot, expr_genes = expr_genes)
  
  # combine into one table
  baseline_preds <- rbindlist(baseline_preds)
  
  return(baseline_preds)
  
}

# function to compute one baseline predictor
compute_baseline_predictor <- function(pred, expt, tss_annot, gene_annot, expr_genes) {
  
  # compute baseline predictors
  output <- switch(
    pred,
    "distToTSS" = computeDistToGene(expt, annot = tss_annot, name = pred, fix_annot = "center"),
    "distToGene" = computeDistToGene(expt, annot = gene_annot, name = pred, fix_annot = "none"),
    "nearestTSS" = nearestFeaturePred(expt, features = tss_annot, name = pred),
    "nearestGene" = nearestFeaturePred(expt, features = gene_annot, name = pred),
    "within100kbTSS" = withinDistFeature(expt, features = tss_annot, dist = 1e+05, name = pred),
    "within100kbGene" = withinDistFeature(expt, features = gene_annot, dist = 1e+05, name = pred),
    "nearestExprTSS" = nearestFeaturePred(
      expt, features = tss_annot[tss_annot$gene %in% expr_genes, ], name = pred),
    "nearestExprGene" = nearestFeaturePred(
      expt, features = gene_annot[gene_annot$gene %in% expr_genes, ], name = pred),
    "within100kbExprTSS" = withinDistFeature(
      expt, features = tss_annot[tss_annot$gene %in% expr_genes, ], dist = 1e+05, name = pred),
    "within100kbExprGene" = withinDistFeature(
      expt, features = gene_annot[gene_annot$gene %in% expr_genes, ], dist = 1e+05, name = pred)
  )

  return(output)
  
}

# compute baseline distance to any gene annotation (TSS or gene)
computeDistToGene <- function(expt, annot, name, fix_annot = c("none", "center", "start", "end")) {
  
  # parse fix argument
  fix_annot <- match.arg(fix_annot)
  
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
    pred_elements = expt$name,
    pred_uid = paste("baseline", name, sep = "."), 
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
  
  # initial column names in expt
  expt_cols <- colnames(expt)
  
  # create GRanges object containing feature coordinates
  features <- features[, 1:4]
  colnames(features) <- c("chr", "start", "stop", "name")
  features <- makeGRangesFromDataFrame(features, keep.extra.columns = TRUE)
  
  # get closest TSS to every CRE in expt
  expt <- nearest_genomic_feature(expt, features = features)
  
  # create predictor whether the gene of each pair is the closest TSS
  output <- data.table(
    expt[, ..expt_cols],
    PredictionCellType = NA_character_,
    pred_elements = expt$name,
    pred_uid = paste("baseline", name, sep = "."), 
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

# candidate enhancers within a certain distance from gene TSS or body baseline predictor
withinDistFeature <- function(expt, features, dist, name) {
  
  # initial column names in expt
  expt_cols <- colnames(expt)
  
  # create GRanges object containing feature coordinates
  features <- features[, 1:4]
  colnames(features) <- c("chr", "start", "stop", "name")
  features <- makeGRangesFromDataFrame(features, keep.extra.columns = TRUE)
  
  # create unique CRE identifiers
  expt$cre_id <- with(expt, paste0(chrom, ":", chromStart, "-", chromEnd))
  
  # get unique CREs create GRanges object containing CRE coordinates
  cres <- unique(expt[, c("chrom" ,"chromStart", "chromEnd", "cre_id")])
  cres <- makeGRangesFromDataFrame(cres, seqnames.field = "chrom", start.field = "chromStart",
                                   end.field = "chromEnd", keep.extra.columns = TRUE)
  
  # extend CREs by dist on both sides to create windows for overlapping
  cres <- resize(cres, width = dist * 2, fix = "center")
  
  # find all features within the specified distance from each CRE
  ovl <- findOverlaps(cres, features)
  
  # get names of CREs and feature pairs within specified distance
  within_dist_pairs <- data.table(cre_id = cres[queryHits(ovl)]$cre_id,
                                  measuredGeneSymbol = features[subjectHits(ovl)]$name,
                                  pred_value = 1)
  
  # add this to expt data
  expt <- merge(expt, within_dist_pairs, by = c("cre_id", "measuredGeneSymbol"), all.x = TRUE)
  expt$pred_value[is.na(expt$pred_value)] <- 0
  
  # reformat for output
  output <- data.table(
    expt[, ..expt_cols],
    PredictionCellType = NA_character_,
    pred_elements = expt$name,
    pred_uid = paste("baseline", name, sep = "."),    
    pred_id = "baseline",
    pred_col = name,
    pred_value = expt$pred_value,
    Prediction = 1
  )
  
}
