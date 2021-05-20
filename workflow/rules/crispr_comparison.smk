# rules to perform comparisons of CRE predictions to CRISPR data

# python helper function to get pred_config file for a given comparison
#def get_pred_config(wildcards):
#  try:
#    pred_config = config["comparisons"][wildcards.comparison]["pred_config"]
#  except:
#    print ("pred_config file not provided, generating one with default parameters.")
#    comparison = wildcards.comparison
#    pred_config = "results/" + comparison + "/pred_config.txt"
#  return pred_config

def get_pred_config(wildcards):
  pred_config = config["comparisons"][wildcards.comparison]["pred_config"]
  if pred_config is None:
    print ("pred_config file not provided, generating one with default parameters.")
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
    gene_universe = lambda wildcards: config["comparisons"][wildcards.comparison]["gene_universe"],
    pred_config = get_pred_config
  output:
    merged = "results/{comparison}/expt_pred_merged.txt"
  log: "results/{comparison}/logs/mergePredictionsWithExperiment.log"
  params:
    pred_names = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].keys(),
    expt_positive_column = lambda wildcards: config["comparisons"][wildcards.comparison]["expt_positive_column"],
    cell_name_mapping = lambda wildcards: config["comparisons"][wildcards.comparison]["cell_name_mapping"]
  conda: "../envs/r_crispr_comparison.yml"
  script:
   "../../workflow/scripts/mergePredictionsWithExperiment.R"
   
# perform comparisons of predictions to experimental data
rule comparePredictionsToExperiment:
  input:
    merged = "results/{comparison}/expt_pred_merged.txt",
    pred_config = lambda wildcards: config["comparisons"][wildcards.comparison]["pred_config"],
  output:
    "results/{comparison}/crispr_comparison.html"
  params:
     pred_names = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].keys(),
     expt_positive_column = lambda wildcards: config["comparisons"][wildcards.comparison]["expt_positive_column"],
     include_missing_predictions = True,
     min_sensitivity = 0.7
  conda: "../envs/r_crispr_comparison.yml"
  script:
    "../../workflow/scripts/comparePredictionsToExperiment.Rmd"
    
