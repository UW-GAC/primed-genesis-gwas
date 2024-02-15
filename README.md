# primed-genesis-gwas
Workflow for running GENESIS on PRIMED data

Uses workflow from https://github.com/AnalysisCommons/genesis_wdl/tree/v1_5 with additional step of VCF conversion from https://github.com/manning-lab/vcfToGds

input | description
--- | ---
vcf_files | Array of vcf files (google bucket paths)
pheno_file | CSV file with phenotypes (google bucket path)
outcome_name | name of the column in pheno_file with the outcome variable
outcome_type | allowed values: Continuous, Dichotomous
covariates_string | comma-separated string with column names of covariates in pheno_file
pheno_id | name of the column in pheno_file with sample identifiers (should match the header of the VCF files) (default "sample_id")
results_file | prefix of output files (default: "gwas")

See the original [README](https://github.com/AnalysisCommons/genesis_wdl/blob/v1_5/README.md) for description of outputs and additional optional inputs
