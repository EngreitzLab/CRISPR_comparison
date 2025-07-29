## CRISPR data for benchmarking enhancer-gene predictive models

This directory contains CRISPR ground truth datasets for benchmarking enhancer predictive models.
Currently the following datasets are available:

- **training_K562**: K562 training data from Gschwind et al., 2025 containing 471 CRISPR positive and
9885 negative E-G pairs.
- **heldout_5_cell_types**: Filtered heldout data from 5 cell types from Gschwind et al., 2025
containing 190 CRISPR positive and 4188 negative E-G pairs

### File format

Each file contains the following columns providing the required data to benchmark predictive models
as well as additional genomic feature annotations.

- chrom [str: chr1, chr2, …, chr22, chrX]: chromosome of perturbed element
- chromStart [int]: start coordinate of perturbed element
- chromEnd [int]: end coordinate of perturbed element
- name [str: for example: PPIF|chr10:80712806-80713306:.]: enhancer-gene pair identifier
- EffectSize [float]: measured and signed effect size on target gene expression
- chrTSS [str: chr1, chr2, …, chr22, chrX]: chromosome of target gene gene’s transcription start site (TSS)
- startTSS [int]: start coordinate of target gene TSS
- endTSS [int]: end coordinate of target gene TSS
- measuredGeneSymbol [str]: gene symbol of target gene
- measuredGeneEnsemblId [str: for example: ENSG00000108179]: Ensembl ID of target gene
- Significant [bool]: whether the effect on target gene is statistically significant (pValueAdjusted < 0.05)
- pValueAdjusted [float]: adjusted p-value for effect size significance 
- PowerAtEffectSize10 [float: 0–1]: power to detect effect size 10%
- PowerAtEffectSize15 [float: 0–1]: power to detect effect size 15%
- PowerAtEffectSize20 [float: 0–1]: power to detect effect size 20%
- PowerAtEffectSize25 [float: 0–1]: power to detect effect size 25%
- PowerAtEffectSize50 [float: 0–1]: power to detect effect size 50%
- ValidConnection [str: for example: “TRUE”, “overlaps potential promoter”, “overlaps target gene exon”]: whether the pair is a valid regulatory connection. If not “TRUE” then listing the reason why this is considered an invalid element-gene pairing, e.g., “overlaps potential promoter”.
- CellType [str: K562, GM12878, HCT116, Jurkat, WTC11]: cell type in which measurement was made
- Reference [str: for example: Nasser et al., 2021 from Ulirsch et al., 2016; Gasperini et al., 2019; …]: full reference of data source, including original study
- Regulated [bool]: whether the element has a significant and negative effect on target gene expression (pValueAdjusted < 0.05 AND EffectSize < 0)
- Dataset [str: for example:  Nasser2021, Gasperini2019, …]: identifier of data source
- distanceToTSS [int]: distance from center of perturbed element to target gene TSS
- measuredGeneUbiquitousExpressed [bool]: whether the target gene is classified as ubiquitously and uniformly expressed
- elementChromatinCategory [str: H3K27ac high, H3K27ac, No H3K27ac, CTCF element, H3K27me3 element]: classification of perturbed element based on chromatin marks
- H3K27ac_peak_overlap  [int: 0, 1]: 1 if the perturbed element overlaps a H3K27ac ChIP-seq peak, otherwise 0
- H3K27ac.RPM [float]: H3K27ac ChIP-seq RPM in perturbed element
- H3K27ac.percentile [float: 0–1]: percentile of H3K27ac signal in perturbed element relative to all genome-wide non-promoter elements in the cell type
- H3K27ac.RPM.expandedRegion [float]: H3K27ac RPM in perturbed element (+/- 150 bp)
- H3K27ac.expandedRegion.percentile [float: 0–1]: percentile of H3K27ac signal in perturbed element (+/- 150 bp) relative to all genome-wide non-promoter elements in the cell type
- DHS.RPM [float]: DNase-seq RPM in perturbed element
- DHS.percentile [float: 0–1]: percentile of DNase-seq signal in perturbed element relative to all genome-wide non-promoter elements in the cell type
- CTCF_peak_overlap [int: 0, 1]: 1 if the perturbed element overlaps a CTCF ChIP-seq peak, otherwise 0
- CTCF.RPM [float]: CTCF ChIP-seq RPM in perturbed element
- H3K27me3_peak_overlap  [int: 0, 1]: 1 if the perturbed element overlaps a H3K27me3 ChIP-seq peak, otherwise 0
- H3K27me3.RPM [float]: H3K27me3 ChIP-seq RPM in perturbed element
- H3K4me1_peak_overlap  [int: 0, 1]: 1 if the perturbed element overlaps a H3K4me31 ChIP-seq peak, otherwise 0
- H3K4me1.RPM [float]: H3K4me1 ChIP-seq RPM in perturbed element
- direct_vs_indirect_negative [float: 0–1]: probability that an observed negative effect size represents is the result of a direct cis-regulatory interaction. Computed based on direct and indirect effects rate columns
- direct_rate_negative [float: 0–1]: estimated rate of direct negative effects based on distance to TSS for tested element-gene pairs
- indirect_rate_negative [float: 0–1]: estimated or imputed rate of indirect negative effects for the dataset to which this pair belongs
- imputed_indirect_rate [bool]: whether the indirect effect rate was imputed (TRUE) or dataset-specific (FALSE)
