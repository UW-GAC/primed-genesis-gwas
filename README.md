# primed-genesis-gwas
Workflow for running GENESIS on PRIMED data

Uses workflow from https://github.com/AnalysisCommons/genesis_wdl/tree/v1_5 with additional step of VCF conversion from https://github.com/manning-lab/vcfToGds

input | description
--- | ---
vcf_files | Array of vcf files (google bucket paths)
pheno_file | CSV file with phenotypes (google bucket path)
outcome | name of the column in pheno_file with the outcome trait
outcome_type | allowed values: quantitative, binary
covariates | comma-separated string with column names of covariates in pheno_file
pheno_id | name of the column in pheno_file with sample identifiers (should match the header of the VCF files) (default "sample_id")
transform | Boolean for whether to rank-normalize residuals and scale, and re-fit null model (default false). Can be paired with other [optional inputs](https://github.com/AnalysisCommons/genesis_wdl/tree/v1_5?tab=readme-ov-file#inputs)
results_prefix | prefix of output files (default "gwas")
genome_build | allowed values: hg38, hg19 (default hg38)
strand | allowed values: +, -. forward, reverse (default +)

See the original [README](https://github.com/AnalysisCommons/genesis_wdl/blob/v1_5/README.md) for description of additional optional inputs
