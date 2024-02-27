# primed-genesis-gwas

Workflow for running GENESIS on PRIMED data

## Data preparation

### Genotype files

The workflow expects genotype data in VCF (gzipped or not). (BCF is currently not accepted but can be added; please submit an issue if you need this option.) The PRIMED data model specifies that genotype data should be in hg38 but this workflow does accept hg19 (see `genome_build` input). The `vcf_files` input is an array, so VCF files may be split by chromosome (this is most efficient as workflow jobs will be scattered by input file). If the genotypes were imputed, set `is_imputed` to `true`.

### Phenotype

The workflow expects a phenotype file in CSV format. A sample identifier column
in this file links the phenotypes and genotypes; the name of this column is controlled
by the `pheno_id` input with default "sample_id" following the PRIMED data model.
If sex is included in the phenotype file, the column should be called "sex" with
values "M" and "F"; GENESIS will use this information to calculate allele frequency
on the sex chromosomes. If age is included, specify the column names with `age_column`
(default "age_at_observation"); this will be used to calculate age-related metadata
for the PRIMED data model.

[This notebook](https://uw-gac.github.io/primed_example_notebooks/analysis/pheno_file_gwas.nb.html) demonstrates how to prepare a phenotype file in R.

### Kinship matrix

For sample sets with recent relatedness, a kinship matrix may be included in either RData or CSV format (see `kinship_matrix` input). The row and column names must match the phenotype and genotype files.

## Workflow steps

1. Convert VCF to GDS using this workflow: https://github.com/manning-lab/vcfToGds
2. Run GENESIS using this workflow: https://github.com/AnalysisCommons/genesis_wdl/tree/v1_5
3. Format output in the [PRIMED GSR data model](https://github.com/UW-GAC/primed_data_models/blob/main/PRIMED_GSR_data_model.json)


## Workflow inputs

input | description
--- | ---
vcf_files | Array of vcf files (google bucket paths)
pheno_file | CSV file with phenotypes (google bucket path)
outcome | Name of the column in pheno_file with the outcome trait
outcome_type | Allowed values: quantitative, binary
outcome_unit | Unit of measurement for the outcome (before transformation)
outcome_definition | A brief description of how the outcome was measured or defined
covariates | Comma-separated string with column names of covariates in pheno_file
kinship_matrix | (optional) Kinship matrix with sample ids as the row and column names. Matricies saved as Rda will load faster, but csv is accepted as well. Rda files should contain a single numeric matrix object.
pheno_id | Name of the column in pheno_file with sample identifiers (should match the header of the VCF files) (default "sample_id")
transform | Boolean for whether to rank-normalize residuals and scale, and re-fit null model (default false). Can be paired with other [optional inputs](https://github.com/AnalysisCommons/genesis_wdl/tree/v1_5?tab=readme-ov-file#inputs)
min_mac | Minimum minor allele count for threshold (default 5)
results_prefix | Prefix of output files (default "gwas")
genome_build | Allowed values: hg38, hg19 (default hg38)
strand | Allowed values: +, -. forward, reverse (default +)
age_column | (optional) Name of the column in pheno_file with age (default: age_of_observation)
consent_code | Consent abbreviation
contributor_contact | Email of the PRIMED contributor who can be contacted for data related questions
genotyping_technology | Allowed values: WGS, WES, genome-wide array, exome array, other array
genotyping_platform | Genotyping platform description including manufacturer, array name, sequencer name
is_imputed | Boolean indicator of whether the analysis was performed using imputed genotypes or dosages
imputation_reference_panel | Allowed values: 1000 Genomes, HRC, TOPMed, Other (required if if_imputed is true)
imputation_reference_panel_detail | Details of the imputation reference panel; e.g. version number or name of panel when imputation_reference_panel = "Other" (required if if_imputed is true)
imputation_quality_filter | Minimum imputation quality value (e.g. Rsq, info) for filtering imputed variants. If no filter, enter value of 0. (required if if_imputed is true)
cohorts | A list of cohorts that collected the samples.
population_descriptor | The concept or classification scheme used to categorize people into populations for this analysis
population_labels | name given to a population that describes or classifies it according to the dimension along which it was identified
population_proportions | proportion of participants from each population in the same order mapping to the values in the population_labels variable
countries_of_recruitment | Reported countries of recruitment
countries_of_birth | (optional) Reported countries of birth
import_tables | A boolean indicating whether data model tables should be imported to the workspace.
overwrite | A boolean indicating whether existing rows in the workspace data tables should be overwritten.
workspace_name | A string with the workspace name. e.g, if the workspace URL is https://anvil.terra.bio/#workspaces/fc-product-demo/Terra-Workflows-Quickstart, the workspace name is "Terra-Workflows-Quickstart"
workspace_namespace | A string with the workspace name. e.g, if the workspace URL is https://anvil.terra.bio/#workspaces/fc-product-demo/Terra-Workflows-Quickstart, the workspace namespace is "fc-product-demo"


See the original [README](https://github.com/AnalysisCommons/genesis_wdl/blob/v1_5/README.md) for description of additional optional inputs.
