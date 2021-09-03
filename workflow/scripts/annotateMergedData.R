## Overlap merged data with genomic features in .bed format. Adds one additional column per
## overlapped feature to merged data. Added columns start with an "overlaps_" prefix.


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

## Define functions --------------------------------------------------------------------------------

# annotate merged data (data.frame) with overlapping features (list of data.frames)
annotate_pairs <- function(merged, features) {
  
  # create GRanges objects from features
  features <- lapply(features, FUN = makeGRangesFromDataFrame, seqnames.field = "V1",
                     start.field = "V2", end.field = "V3", strand.field = "V6",
                     starts.in.df.are.0based = TRUE)
  
  # convert to GRangesList
  features <- GRangesList(features)
  
  # create identifier for CREs in merged data
  merged <- unite(merged, col = "cre_id", chrom, chromStart, chromEnd, sep = "_", remove = FALSE)
  
  # extract CRE coordinates from merged data and create GRanges object
  cre_coords <- merged %>% 
    select(chr = chrom, start = chromStart, end = chromEnd, cre_id) %>% 
    distinct() %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  
  # get CREs that overlap with any features and convert to data frame
  overlaps <- lapply(features, FUN = findOverlaps, query = cre_coords, ignore.strand = TRUE)
  
  # create data frame with all CREs with "empty" overlaps column
  all_cres <- data.table(cre_id = cre_coords$cre_id, overlap = FALSE)
  
  # fill in overlaps for every feature
  cre_overlaps <- lapply(overlaps, FUN = function(ovl, all_cres) {
    all_cres[unique(queryHits(ovl)), "overlap"] <- TRUE
    return(all_cres)
  }, all_cres = all_cres)
  
  # convert to one data frame
  cre_overlaps <- rbindlist(cre_overlaps, idcol = "feature")
  
  # convert to wide format and add prefix to feature columns
  cre_overlaps <- dcast(cre_overlaps, cre_id ~ feature, value.var = "overlap")
  colnames(cre_overlaps)[-1] <- paste0("overlaps_", colnames(cre_overlaps)[-1])
  
  # add overlaps to original merged data to create output
  output <- merge(merged, cre_overlaps, by = "cre_id")
  output <- output[, -"cre_id"]
  
  return(output)
  
}

## Annotate merged data ----------------------------------------------------------------------------

# load merged data
merged <- fread(snakemake@input$merged)

# get annotation feature files
config <- snakemake@config$comparisons[[snakemake@wildcards$comparison]]
feature_files <- config$annotation_features

# if no annotation features are provided, simply write merged data to output file else annotate
if (is.null(feature_files)) {
  
  message("No annotation features provided.")
  readr::write_tsv(merged, file = snakemake@output[[1]])
  
} else {

  # load annotation features
  message("Loading annotation features...")
  features <- lapply(feature_files, FUN = fread)

  # overlap merged data with features
  message("Annotating merged data...")
  output <- annotate_pairs(merged, features = features)

  # write to new output file
  message("Writing to output file...")
  readr::write_tsv(output, file = snakemake@output[[1]])

  message("Done!")

}
