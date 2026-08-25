# rules to perform comparisons of CRE predictions to CRISPR data

# valid prediction file formats and format assumed if not specified in config
# (keep in sync with 'possible_formats' in crisprComparisonLoadInputData.R)
PRED_FILE_FORMATS = ["generic", "ENCODE", "IGVF"]
DEFAULT_PRED_FILE_FORMAT = "generic"

# get all prediction files and concatenate them into one array
def get_predictions(wildcards):
  preds = config["comparisons"][wildcards.comparison]["pred"]
  preds_array = []
  for value in preds.values():
    if isinstance(value, list):
      preds_array.extend(value)
    else:
      preds_array.append(value)
  return preds_array

# get pred_config file if specified in config, else create name for default pred_config file
def get_pred_config(wildcards):
  pred_config = config["comparisons"][wildcards.comparison]["pred_config"]
  if pred_config is None:
    comparison = wildcards.comparison
    pred_config = "results/" + comparison + "/pred_config.txt"
  return pred_config

# get prediction file format(s) for a comparison. this is either one format for all predictors, or
# a mapping keyed by the names of 'pred' specifying the format of each predictor individually.
# predictors omitted from the mapping (or set to Null) use DEFAULT_PRED_FILE_FORMAT. all files
# listed for one predictor are assumed to have the same format
def get_pred_file_format(wildcards):
  comparison = config["comparisons"][wildcards.comparison]
  formats = comparison.get("pred_file_format")

  # if no formats are specified, all files are loaded using the default format
  if formats is None:
    return DEFAULT_PRED_FILE_FORMAT

  # if formats are specified per predictor, fill in the default for any missing predictors
  if isinstance(formats, dict):
    pred_names = list(comparison["pred"].keys())
    unknown = [name for name in formats if name not in pred_names]
    if unknown:
      raise ValueError(
        "'pred_file_format' in comparison '" + wildcards.comparison + "' contains predictors not "
        "found in 'pred': " + ", ".join(unknown))
    multiple = [name for name, fmt in formats.items() if isinstance(fmt, list)]
    if multiple:
      raise ValueError(
        "Only one 'pred_file_format' can be specified per predictor, but multiple were found in "
        "comparison '" + wildcards.comparison + "' for: " + ", ".join(multiple))
    formats = {name: formats.get(name) or DEFAULT_PRED_FILE_FORMAT for name in pred_names}
    invalid = [fmt for fmt in formats.values() if fmt not in PRED_FILE_FORMATS]
  else:
    invalid = [formats] if formats not in PRED_FILE_FORMATS else []

  # raise an error early if any invalid formats were specified
  if invalid:
    raise ValueError(
      "Invalid 'pred_file_format' in comparison '" + wildcards.comparison + "': " +
      ", ".join(map(str, invalid)) + ". Valid formats: " + ", ".join(PRED_FILE_FORMATS))

  return formats

# get optional input parameter if they are specified in config
def get_optional_parameter(wildcards, param, default=[]):
  try:
    param = config["comparisons"][wildcards.comparison][param]
  except KeyError:
    param = None
  if param is None:
    param = default
  else:
    if type(param) is dict:
      param = param.values()
  return param

## RULES -------------------------------------------------------------------------------------------

# create minimal pred_config file with default values
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
    predictions = get_predictions,
    experiment  = lambda wildcards: config["comparisons"][wildcards.comparison]["expt"],
    tss_universe = lambda wildcards: config["comparisons"][wildcards.comparison]["tss_universe"],
    gene_universe = lambda wildcards: config["comparisons"][wildcards.comparison]["gene_universe"],
    pred_config = get_pred_config,
    cell_type_mapping = lambda wildcards: get_optional_parameter(wildcards, "cell_type_mapping"),
    expressed_genes = lambda wildcards: get_optional_parameter(wildcards, "expressed_genes")
  output:
    merged = temp("results/{comparison}/expt_pred_merged.txt.gz")
  params:
    pos_col = "Regulated",
    include_col = lambda wildcards: get_optional_parameter(wildcards, "include_col", None),
    pred_format = get_pred_file_format,
    filter_tss = lambda wildcards: get_optional_parameter(wildcards, "filter_pred_tss", True),
    filter_include_col = False
  log: "results/{comparison}/logs/mergePredictionsWithExperiment.log"
  conda: "../envs/r_crispr_comparison.yml"
  resources:
    mem_mb = 32000
  script:
   "../../workflow/scripts/mergePredictionsWithExperiment.R"
   
# annotate enhancers in merged data with overlapping genomic features and assays
rule annotateEnhFeatures:
  input:
    merged = "results/{comparison}/expt_pred_merged.txt.gz",
    gene_features = lambda wildcards: get_optional_parameter(wildcards, "gene_features"),
    enh_features = lambda wildcards: get_optional_parameter(wildcards, "enh_features"),
    enh_assays = lambda wildcards: get_optional_parameter(wildcards, "enh_assays")
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
     dist_bins_kb = lambda wildcards: config["comparisons"][wildcards.comparison]["dist_bins_kb"],
     include_col = lambda wildcards: get_optional_parameter(wildcards, "include_col", None)
  conda: "../envs/r_crispr_comparison.yml"
  resources:
    mem_mb = 32000,
    runtime = "6h"
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
