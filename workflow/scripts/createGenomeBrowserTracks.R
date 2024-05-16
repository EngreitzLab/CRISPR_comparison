## Create genome browser tracks for EPbenchmarking CRISPR dataset

# save.image("tracks.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

## Define functions --------------------------------------------------------------------------------

# function to create crispr tracks for one cell type
create_crispr_tracks <- function(merged, outdir, cell_type) {
  
  # extract all crispr E-G pairs for the given cell type
  crispr <- merged %>% 
    filter(ExperimentCellType == cell_type) %>% 
    select(chrom, chromStart, chromEnd, chrTSS, startTSS, endTSS, name, EffectSize, Regulated) %>% 
    distinct()
  
  # create bed track with all tested crispr elements
  elements_bed <- crispr %>% 
    mutate(strand = ".", name = paste0(chrom, ":", chromStart, "-", chromEnd), score = 0) %>% 
    select(chrom, chromStart, chromEnd, name, score, strand) %>% 
    distinct()
  
  # create bedpe track for all E-G pairs
  bedpe <- crispr %>% 
    mutate(strand1 = ".", strand2 = ".") %>% 
    select(chrom1 = chrom, start1 = chromStart, end1 = chromEnd, chrom2 = chrTSS, start2 = startTSS,
           end2 = endTSS, name, score = EffectSize, strand1, strand2, Regulated)
  
  # split into tracks for positives and negatives
  pos_bedpe <- select(filter(bedpe, Regulated == TRUE),  -Regulated)
  neg_bedpe <- select(filter(bedpe, Regulated == FALSE), -Regulated)
  
  # output file paths
  elements_bed_outfile <- file.path(outdir, cell_type, "crispr_elements.bed")
  pos_bedpe_outfile <- file.path(outdir, cell_type, "crispr_positives.bedpe")
  neg_bedpe_outfile <- file.path(outdir, cell_type, "crispr_negatives.bedpe")
  
  # write tracks to output files
  dir.create(file.path(outdir, cell_type), recursive = TRUE, showWarnings = FALSE)
  fwrite(elements_bed, file = elements_bed_outfile, sep = "\t", col.names = FALSE, quote = FALSE)
  fwrite(pos_bedpe, file = pos_bedpe_outfile, sep = "\t", col.names = FALSE, quote = FALSE)
  fwrite(neg_bedpe, file = neg_bedpe_outfile, sep = "\t", col.names = FALSE, quote = FALSE)
  
}

## Create genome browser tracks --------------------------------------------------------------------

# load merged data
merged <- fread(snakemake@input$merged)

# create CRISPR genome browser tracks for all cell types
for (i in unique(merged$ExperimentCellType)) {
  create_crispr_tracks(merged, outdir = snakemake@output[[1]], cell_type = i)
}
