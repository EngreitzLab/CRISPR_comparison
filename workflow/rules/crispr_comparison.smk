# rules to perform comparisons of CRE predictions to CRISPR data

# get pred_config file if specified in config, else create name for default pred_config file
def get_pred_config(wildcards):
  pred_config = config["comparisons"][wildcards.comparison]["pred_config"]
  if pred_config is None:
    comparison = wildcards.comparison
    pred_config = "results/" + comparison + "/pred_config.txt"
  return pred_config

# get cell type mapping files if they are specified in config
def get_cell_type_mappings(wildcards):
  ct_map_config = config["comparisons"][wildcards.comparison]["cell_type_mapping"]
  if ct_map_config is None:
    ct_map = []
  else:
    ct_map = ct_map_config.values()
  return ct_map 

# get annotation features if they are specified in config
def get_annot_features(wildcards):
  feat_config = config["comparisons"][wildcards.comparison]["annotation_features"]
  if feat_config is None:
    feat = []
  else:
    feat = feat_config.values()
  return feat

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
    cell_type_mapping = get_cell_type_mappings
  output:
    merged = "results/{comparison}/expt_pred_merged.txt.gz"
  log: "results/{comparison}/logs/mergePredictionsWithExperiment.log"
  conda: "../envs/r_crispr_comparison.yml"
  script:
   "../../workflow/scripts/mergePredictionsWithExperiment.R"
   
# annotate merged data with overlapping genomic features
rule annotateMergedData:
  input:
    merged = "results/{comparison}/expt_pred_merged.txt.gz",
    annot = get_annot_features
  output:
    "results/{comparison}/expt_pred_merged_annot.txt.gz"
  conda: "../envs/r_crispr_comparison.yml"
  script:
    "../../workflow/scripts/annotateMergedData.R"
   
# perform comparisons of predictions to experimental data
rule comparePredictionsToExperiment:
  input:
    merged = "results/{comparison}/expt_pred_merged_annot.txt.gz",
    pred_config = get_pred_config
  output:
    "results/{comparison}/crispr_comparison.html"
  params:
     pred_names = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].keys(),
     include_missing_predictions = True,
     min_sensitivity = 0.7
  conda: "../envs/r_crispr_comparison.yml"
  script:
    "../../workflow/scripts/comparePredictionsToExperiment.Rmd"
