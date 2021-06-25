# rules to perform comparisons of CRE predictions to CRISPR data

def get_pred_config(wildcards):
  pred_config = config["comparisons"][wildcards.comparison]["pred_config"]
  if pred_config is None:
    comparison = wildcards.comparison
    pred_config = "results/" + comparison + "/pred_config.txt"
  return pred_config

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
    pred_config = get_pred_config
  output:
    merged = "results/{comparison}/expt_pred_merged.txt"
  log: "results/{comparison}/logs/mergePredictionsWithExperiment.log"
  params:
    pred_names = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].keys(),
    cell_type_mapping = lambda wildcards: config["comparisons"][wildcards.comparison]["cell_type_mapping"]
  conda: "../envs/r_crispr_comparison.yml"
  script:
   "../../workflow/scripts/mergePredictionsWithExperiment.R"
   
# perform comparisons of predictions to experimental data
rule comparePredictionsToExperiment:
  input:
    merged = "results/{comparison}/expt_pred_merged.txt",
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
    
