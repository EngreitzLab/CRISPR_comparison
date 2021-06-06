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
