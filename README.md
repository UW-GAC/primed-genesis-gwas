# primed-genesis-gwas
Workflow for running GENESIS on PRIMED data

Uses workflow from https://github.com/AnalysisCommons/genesis_wdl/tree/v1_5 with additional step of VCF conversion from https://github.com/manning-lab/vcfToGds

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


See the original [README](https://github.com/AnalysisCommons/genesis_wdl/blob/v1_5/README.md) for description of additional optional inputs
