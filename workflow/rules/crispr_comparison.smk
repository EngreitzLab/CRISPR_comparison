# rules to perform comparisons of CRE predictions to CRISPR data

# merge predictions with experimental data
rule mergePredictionsWithExperiment:
  input:
    predictions = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].values(),
    experiment  = lambda wildcards: config["comparisons"][wildcards.comparison]["expt"],
    pred_config = lambda wildcards: config["comparisons"][wildcards.comparison]["pred_config"]
  output:
    merged = "results/{comparison}/expt_pred_merged.txt"
  log: "results/{comparison}/logs/mergePredictionsWithExperiment.log"
  params:
    pred_names = lambda wildcards: config["comparisons"][wildcards.comparison]["pred"].keys(),
    pred_threshold = lambda wildcards: config["comparisons"][wildcards.comparison]["pred_threshold"],
    expt_positive_column = lambda wildcards: config["comparisons"][wildcards.comparison]["expt_positive_column"],
    ignore_expt_missing_predictions = lambda wildcards: config["comparisons"][wildcards.comparison]["ignore_expt_missing_predictions"],
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
     min_sensitivity = lambda wildcards: config["comparisons"][wildcards.comparison]["min_sensitivity"],
     pred_threshold = lambda wildcards: config["comparisons"][wildcards.comparison]["pred_threshold"]
  conda: "../envs/r_crispr_comparison.yml"
  script:
    "../../workflow/scripts/comparePredictionsToExperiment.Rmd"
    
