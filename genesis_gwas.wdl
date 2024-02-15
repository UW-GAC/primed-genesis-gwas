version 1.0

# forked the original repo to update some syntax to WDL version 1.0
import "https://raw.githubusercontent.com/UW-GAC/genesis_wdl/v1_5/genesis_GWAS.wdl" as genesis
import "https://raw.githubusercontent.com/manning-lab/vcfToGds/main/vcfToGds.wdl" as gds

workflow genesis_gwas {
    input {
        Array[File] vcf_files
        File pheno_file
        String outcome_name
        String outcome_type
        String covariates_string
        String pheno_id = "sample_id"
	    String results_file = "gwas"
    }

    call gds.vcfToGds_wf {
        input:
            vcf_files = vcf_files
    }

    call genesis.genesis_gwas_wf {
        input:
            these_genotype_files = vcfToGds_wf.gds_files,
            this_pheno_file = pheno_file,
            this_outcome_name = outcome_name,
            this_outcome_type = outcome_type,
            this_covariates_string = covariates_string,
            this_pheno_id = pheno_id,
            this_results_file = results_file,
            this_test_type = "Single"
    }

    output {
        File null_model = genesis_gwas_wf.null_model
        Array[File] raw_association_files = genesis_gwas_wf.raw_association_files
        File all_summary_statistics = genesis_gwas_wf.all_summary_statistics
        File top_summary_statistics = genesis_gwas_wf.top_summary_statistics
        File summary_plots = genesis_gwas_wf.summary_plots
    }
}