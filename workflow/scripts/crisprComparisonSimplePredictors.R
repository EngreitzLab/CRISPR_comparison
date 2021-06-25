## A couple of simple baseline predictors

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

## HELPER FUNCTIONS ================================================================================

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
