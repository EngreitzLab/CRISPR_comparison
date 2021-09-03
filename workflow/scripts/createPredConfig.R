
## Create a pred_config.txt file with default parameters for comparisons where none is provided

message("pred_config file not provided, generating one with default parameters.")

# create pred_config table with default values for each predictor set. this assumes that the column
# with the predictor value is called "Score" and it is treated like a quantitative score, where
# higher values correspond to higher confidence.
pred_config <- data.frame(
  pred_id = snakemake@params$pred_names,
  pred_col = "Score",
  boolean = FALSE,
  alpha = NA,
  aggregate_function = "sum",
  fill_value = 0,
  inverse_predictor = FALSE,
  pred_name_long = paste0(snakemake@params$pred_names, ".Score"),
  color = NA,
  stringsAsFactors = FALSE
  )

# save to output file
write.table(pred_config, file = snakemake@output[[1]], row.names = FALSE, quote = FALSE, sep = "\t")
