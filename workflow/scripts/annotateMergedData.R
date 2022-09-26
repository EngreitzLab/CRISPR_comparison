## Overlap merged data with genomic features in .bed format and chromatin assays in .bam files.
## Adds one additional column per overlapped feature to merged data. Added columns start with an
## "overlaps_" prefix

# save.image("annot.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(rtracklayer)
})

## Define functions --------------------------------------------------------------------------------

# annotate EG pairs in merged data with gene features  (list of data.frames)
annotate_gene_features <- function(merged, gene_features) {
  
  # check that gene features have distinct feature names
  all_colnames <- unlist(lapply(gene_features, FUN = colnames))
  all_feature_names <- all_colnames[all_colnames != "gene"]
  if (any(table(all_feature_names) > 1)) {
    stop("Gene feature files do not have unique feature names.", call. = FALSE)
  }
  
  # merge all gene features in to one table
  all_gene_features <- purrr::reduce(gene_features, full_join, by = "gene")
  
  # set nee column names
  colnames(all_gene_features)[-1] <- paste0("gene_feature_", colnames(all_gene_features)[-1])
  
  # add gene features to merged data
  output <- left_join(merged, all_gene_features, by = c("measuredGeneSymbol" = "gene"))
  
  return(output)
  
}

# annotate merged data (data.frame) with overlapping enhancer features (list of data.frames)
annotate_enh_features <- function(merged, enh_features) {
  
  # create GRanges objects from features
  enh_features <- lapply(enh_features, FUN = makeGRangesFromDataFrame, seqnames.field = "V1",
                     start.field = "V2", end.field = "V3", strand.field = "V6",
                     starts.in.df.are.0based = TRUE)
  
  # convert to GRangesList
  enh_features <- GRangesList(enh_features)
  
  # create identifier for CREs in merged data
  merged <- unite(merged, col = "cre_id", chrom, chromStart, chromEnd, sep = "_", remove = FALSE)
  
  # extract CRE coordinates from merged data and create GRanges object
  cre_coords <- merged %>% 
    select(chr = chrom, start = chromStart, end = chromEnd, cre_id) %>% 
    distinct() %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  
  # get CREs that overlap with any features and convert to data frame
  overlaps <- lapply(enh_features, FUN = findOverlaps, query = cre_coords, ignore.strand = TRUE)
  
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
  colnames(cre_overlaps)[-1] <- paste0("enh_feature_", colnames(cre_overlaps)[-1])
  
  # add overlaps to original merged data to create output
  output <- merge(merged, cre_overlaps, by = "cre_id")
  output <- output[, -"cre_id"]
  
  return(output)
  
}

# annotate enhancers by overlapping them with enhancer assay bam files
annotate_enh_assays <- function(merged, enh_assays, normalize = TRUE) {
  
  # create identifier for CREs in merged data
  merged <- unite(merged, col = "cre_id", chrom, chromStart, chromEnd, sep = "_", remove = FALSE)
  
  # extract CRE coordinates from merged data and create GRanges object
  cre_coords <- merged %>% 
    select(chr = chrom, start = chromStart, end = chromEnd, cre_id) %>% 
    distinct() %>% 
    makeGRangesFromDataFrame(keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  
  # count reads in each CRE
  cre_assay_reads <- enh_assays %>%
    lapply(FUN = countOverlaps, query = cre_coords, ignore.strand = TRUE) %>%
    bind_cols()
  
  # combine with enhancer coordinates and convert to long format
  cre_assay_reads <- cre_coords %>%
    data.frame(stringsAsFactors = FALSE) %>%
    select(-strand) %>%
    dplyr::rename(chrom = seqnames, chromStart = start, chromEnd = end) %>%
    bind_cols(cre_assay_reads) %>%
    pivot_longer(cols = -c(1:5), names_to = "assay", values_to = "reads")
  
  # normalize for sequencing depth and enhancer size
  if (normalize == TRUE) {
    
    # total reads per assay
    total_reads <- vapply(enh_assays, FUN = length, FUN.VALUE = integer(1))
    total_reads <- data.frame(assay = names(total_reads), total_reads = total_reads,
                              row.names = NULL, stringsAsFactors = FALSE)
    
    # add total reads to read counts and normalize by sequencing depth
    cre_assay_reads <- cre_assay_reads %>%
      left_join(total_reads, by = "assay") %>% 
      mutate(reads = reads / (total_reads / 1e6)) %>% 
      select(-total_reads)
    
  }
  
  # convert to wide format add to merged data to create output
  output <- cre_assay_reads %>%
    mutate(assay = paste0("enh_assay_", assay)) %>% 
    select(-c(chrom, chromStart, chromEnd, width)) %>%
    pivot_wider(names_from = "assay", values_from = "reads") %>%
    left_join(x = merged, y = ., by = "cre_id")
  
  # remove cre_id column 
  output <- output[, -"cre_id"]
  
  return(output)
  
}

## Annotate merged data ----------------------------------------------------------------------------

# get annotation features and assay files
config <- snakemake@config$comparisons[[snakemake@wildcards$comparison]]
gene_features_files <- config$gene_features
enh_features_files <- config$enh_features
enh_assays_files <- config$enh_assays

# if no annotation features are provided, simply copy the input file
if (is.null(c(gene_features_files, enh_features_files, enh_assays_files))) {
  
  message("No annotation features or assays provided.")
  invisible(file.copy(from = snakemake@input$merged, to = snakemake@output[[1]], overwrite = TRUE))
  
} else {
  
  # load merged data
  message("Loading merged data...")
  merged <- fread(snakemake@input$merged)
  
  # annotate E-G pairs based on gene features
  if (!is.null(gene_features_files)) {
    
    # load annotation features
    message("Loading gene features...")
    gene_features <- lapply(gene_features_files, FUN = fread)
    
    # annotate genes in merged data with gene features
    message("Annotating merged data with gene features...")
    merged <- annotate_gene_features(merged, gene_features = gene_features)
    
  }

  # annotate enhancers with overlapping genomic features if provided
  if (!is.null(enh_features_files)) {
    
    # load annotation features
    message("Loading enhancer features...")
    enh_features <- lapply(enh_features_files, FUN = fread)
    
    # overlap merged data with features
    message("Annotating merged data with enhancer features...")
    merged <- annotate_enh_features(merged, enh_features = enh_features)
    
  }
  
  # annotate enhancers with assay reads if provided
  if (!is.null(enh_assays_files)) {
    
    # load assay bam files
    message("Loading enhancer assay bam files...")
    assays <- lapply(enh_assays_files, FUN = readGAlignments)
    
    # annotate merged data with read counts from assay bam files
    message("Annotating merged data with enhancer assay reads...")
    merged <- annotate_enh_assays(merged, enh_assays = assays, normalize = TRUE)
    
  }

  # write to new output file
  message("Writing to output file...")
  fwrite(merged, file = snakemake@output[[1]], sep = "\t", na = "NA")

  message("Done!")

}
