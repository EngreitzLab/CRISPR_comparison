# rules to perform comparisons of CRE predictions to CRISPR data

# get pred_config file if specified in config, else create name for default pred_config file
def get_pred_config(wildcards):
  pred_config = config["comparisons"][wildcards.comparison]["pred_config"]
  if pred_config is None:
    comparison = wildcards.comparison
    pred_config = "results/" + comparison + "/pred_config.txt"
  return pred_config

# get optional input files if they are specified in config
def get_optional_file(wildcards, file):
  try:
    file = config["comparisons"][wildcards.comparison][file]
  except KeyError:
    file = None
  if file is None:
    file = []
  else:
    if type(file) is dict:
      file = file.values()
  return file

## RULES -------------------------------------------------------------------------------------------

# create pred_config file with default values
rule createPredConfig:
  output:
    "results/{comparison}/pred_config.txt"
  params:
    pred_names = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].keys()
  conda: "../envs/r_crispr_comparison.yml"
  script:
    "../../workflow/scripts/createPredConfig.R"

# merge predictions with experimental data
rule mergePredictionsWithExperiment:
  input:
    predictions = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].values(),
    experiment  = lambda wildcards: config["comparisons"][wildcards.comparison]["expt"],
    tss_universe = lambda wildcards: config["comparisons"][wildcards.comparison]["tss_universe"],
    gene_universe = lambda wildcards: config["comparisons"][wildcards.comparison]["gene_universe"],
    pred_config = get_pred_config,
    cell_type_mapping = lambda wildcards: get_optional_file(wildcards, "cell_type_mapping"),
    expressed_genes = lambda wildcards: get_optional_file(wildcards, "expressed_genes")
  output:
    merged = temp("results/{comparison}/expt_pred_merged.txt.gz")
  params:
    pos_col = "Regulated"
  log: "results/{comparison}/logs/mergePredictionsWithExperiment.log"
  conda: "../envs/r_crispr_comparison.yml"
  resources:
    mem_mb = 48000
  script:
   "../../workflow/scripts/mergePredictionsWithExperiment.R"
   
# annotate enhancers in merged data with overlapping genomic features and assays
rule annotateEnhFeatures:
  input:
    merged = "results/{comparison}/expt_pred_merged.txt.gz",
    gene_features = lambda wildcards: get_optional_file(wildcards, "gene_features"),
    enh_features = lambda wildcards: get_optional_file(wildcards, "enh_features"),
    enh_assays = lambda wildcards: get_optional_file(wildcards, "enh_assays")
  output:
    "results/{comparison}/expt_pred_merged_annot.txt.gz"
  conda: "../envs/r_crispr_comparison.yml"
  resources:
    mem_mb = 32000
  script:
    "../../workflow/scripts/annotateMergedData.R"
   
# perform comparisons of predictions to experimental data
rule comparePredictionsToExperiment:
  input:
    merged = "results/{comparison}/expt_pred_merged_annot.txt.gz",
    pred_config = get_pred_config
  output: "results/{comparison}/{comparison}_crispr_comparison.html"
  params:
     pred_names = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].keys(),
     include_missing_predictions = True,
     pos_col = "Regulated",
     min_sensitivity = 0.7,
     dist_bins_kb = lambda wildcards: config["comparisons"][wildcards.comparison]["dist_bins_kb"]
  conda: "../envs/r_crispr_comparison.yml"
  resources:
    mem_mb = 8000
  script:
    "../../workflow/scripts/comparePredictionsToExperiment.Rmd"
    
# create genome browser tracks
rule createGenomeBrowserTracks:
  input: 
    merged = "results/{comparison}/expt_pred_merged_annot.txt.gz"
  output: directory("results/{comparison}/genome_browser_tracks")
  conda: "../envs/r_crispr_comparison.yml"
  resources:
    mem_mb = 8000  
  script:
    "../../workflow/scripts/createGenomeBrowserTracks.R"    
