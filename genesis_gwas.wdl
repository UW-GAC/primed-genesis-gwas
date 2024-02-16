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
        String strand = "+"
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

    Array[Array[File]] stat_files = [genesis_gwas_wf.raw_association_files, [genesis_gwas_wf.all_summary_statistics, genesis_gwas_wf.top_summary_statistics]]
    scatter(f in flatten(stat_files)) {
        call file_in_data_model {
            input:
                csv_file = f,
                outcome_type = outcome_type,
                strand = strand
        }
    }

    output {
        File null_model = genesis_gwas_wf.null_model
        Array[File] summary_statistics = file_in_data_model.tsv_file
        File summary_plots = genesis_gwas_wf.summary_plots
    }
}


task file_in_data_model {
    input {
        File csv_file
        String outcome_type
        String strand
    }

    command <<<
        Rscript -e "\
        library(dplyr); \
        library(readr); \
        gsr <- read_csv('~{csv_file}')[,-1]; \
        gsr <- select(gsr, SNPID=snpID, chromosome=chr, position=pos, effect_allele=ref, other_allele=alt, effect_allele_freq=freq, mac=MAC, p_value=ends_with('pval'), beta=Est, se=Est.SE); \
        gsr <- mutate(gsr, strand='~{strand}', beta_ci_lower=(beta + qnorm((1-0.95)*0.05)*se), beta_ci_upper=(beta + qnorm(1-((1-0.95))*0.05)*se), p_value_log10=-log10(p_value)); \
        if ('~{outcome_type}' == 'Dichotomous') { \
            gsr <- mutate(gsr, odds_ratio=exp(beta), OR_ci_lower=exp(beta_ci_lower), OR_ci_upper=exp(beta_ci_upper)); \
        }; \
        out_file <- sub('.csv', '.tsv', csv_file, fixed=TRUE); \
        write_tsv(gsr, out_file); \
        "
    >>>

    output {
        File tsv_file = sub(basename(csv_file), ".csv", ".tsv")
    }

    runtime {
        docker: "ghcr.io/anvilproject/anvil-rstudio-bioconductor:3.18.0"
    }
}
