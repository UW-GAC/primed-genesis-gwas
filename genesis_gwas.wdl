version 1.0

# forked the original repo to update some syntax to WDL version 1.0
import "https://raw.githubusercontent.com/UW-GAC/genesis_wdl/v1_5/genesis_GWAS.wdl" as genesis
import "https://raw.githubusercontent.com/manning-lab/vcfToGds/main/vcfToGds.wdl" as gds
import "https://raw.githubusercontent.com/UW-GAC/primed-file-checks/main/validate_gsr_model.wdl" as validate

workflow genesis_gwas {
    input {
        Array[File] vcf_files
        File pheno_file
        String outcome
        String outcome_type
        String outcome_unit
        String outcome_definition
        String covariates
        File? kinship_matrix
        String pheno_id = "sample_id"
        Boolean transform = false
        Int min_mac = 5
        String results_prefix = "gwas"
        String genome_build = "hg38"
        String strand = "+"
        String age_column = "age_at_observation"
        String model_url = "https://raw.githubusercontent.com/UW-GAC/primed_data_models/main/PRIMED_GSR_data_model.json"
        String workspace_name
        String workspace_namespace
        Boolean import_tables = true
        Boolean overwrite = true
    }

    call gds.vcfToGds_wf {
        input:
            vcf_files = vcf_files
    }

    call genesis.genesis_gwas_wf {
        input:
            these_genotype_files = vcfToGds_wf.gds_files,
            this_pheno_file = pheno_file,
            this_outcome_name = outcome,
            this_outcome_type = if outcome_type == "binary" then "Dichotomous" else "Continuous",
            this_covariates_string = covariates,
            this_kinship_matrix = kinship_matrix,
            this_pheno_id = pheno_id,
            this_transform = if transform then "transform" else "none",
            this_min_mac = min_mac,
            this_results_file = results_prefix,
            this_test_type = "Single",
            this_genome_build = genome_build
    }

    Array[Array[File]] stat_files = [genesis_gwas_wf.raw_association_files, [genesis_gwas_wf.all_summary_statistics, genesis_gwas_wf.top_summary_statistics]]
    scatter(f in flatten(stat_files)) {
        call file_in_data_model {
            input:
                csv_file = f,
                trait_type = outcome_type,
                strand = strand,
                results_prefix = results_prefix
        }
    }

    String method = if defined(kinship_matrix) then (if outcome_type == "binary" then "GLMM" else "LMM") else (if outcome_type == "binary" then "logistic regression" else "linear regression")

    call prepare_gsr_data_model {
        input:
            pheno_file = pheno_file,
            data_files = file_in_data_model.tsv_file,
            md5sum = file_in_data_model.md5sum,
            chromosome = file_in_data_model.chromosome,
            n_variants = file_in_data_model.n_variants,
            trait = outcome,
            trait_type = outcome_type,
            trait_unit = outcome_unit,
            trait_transformation = if transform then "rank-normalized" else "none",
            trait_definition = outcome_definition,
            covariates = covariates,
            reference_assembly = if genome_build == "hg38" then "GRCh38" else "GRCh37",
            min_MAC_filter = min_mac,
            age_column = age_column,
            sex_column = "sex",
            female_value = "F",
            analysis_method = method,
            analysis_software = "GENESIS"
    }

    call validate.validate_gsr_model {
        input: table_files = prepare_gsr_data_model.table_files,
               model_url = model_url,
               workspace_name = workspace_name,
               workspace_namespace = workspace_namespace,
               overwrite = overwrite,
               import_tables = import_tables
    }

    output {
        File null_model = genesis_gwas_wf.null_model
        Array[File] summary_statistics = file_in_data_model.tsv_file
        File summary_plots = genesis_gwas_wf.summary_plots
        File validation_report = validate_gsr_model.validation_report
        Array[File]? tables = validate_gsr_model.tables
    }
}


task file_in_data_model {
    input {
        File csv_file
        String trait_type
        String strand
        String results_prefix
    }

    Int disk_size = ceil(3*(size(csv_file, "GB"))) + 10
    String out_file = sub(basename(csv_file), ".csv", ".tsv")

    command <<<
        Rscript -e "\
        library(dplyr); \
        library(readr); \
        gsr <- read_csv('~{csv_file}')[,-1]; \
        gsr <- select(gsr, SNPID=snpID, chromosome=chr, position=pos, \
            effect_allele=ref, other_allele=alt, effect_allele_freq=freq, mac=MAC, \
            p_value=ends_with('pval'), beta=Est, se=Est.SE); \
        gsr <- mutate(gsr, strand='~{strand}', \
            beta_ci_lower=(beta + qnorm((1-0.95)*0.05)*se), \
            beta_ci_upper=(beta + qnorm(1-((1-0.95))*0.05)*se), \
            p_value_log10=-log10(p_value)); \
        if ('~{trait_type}' == 'binary') { \
            gsr <- mutate(gsr, odds_ratio=exp(beta), \
                OR_ci_lower=exp(beta_ci_lower), \
                OR_ci_upper=exp(beta_ci_upper)); \
        }; \
        write_tsv(gsr, '~{out_file}'); \
        writeLines(as.character(nrow(gsr)), 'n_variants.txt'); \
        chr_string <- sub('^~{results_prefix}\\\\.', '', sub('\\\\.tsv\\\\.gz$', '', '~{out_file}')); \
        if (grepl('variants', chr_string)) chr <- 'ALL' else chr <- chr_string; \
        writeLines(chr, 'chromosome.txt')
        "
        md5sum ~{out_file} | cut -d " " -f 1 | sed 's/ //g' > md5sum.txt
    >>>

    output {
        File tsv_file = out_file
        String md5sum = read_string("md5sum.txt")
        String chromosome = read_string("chromosome.txt")
        Int n_variants = read_int("n_variants.txt")
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.16.0"
        disks: "local-disk " + disk_size + " SSD"
        memory: "16 GB"
    }
}


task prepare_gsr_data_model {
    input {
        Array[String] data_files
        Array[String] md5sum
        Array[String] chromosome
        Array[Int] n_variants
        String consent_code
        String contributor_contact
        String trait
        String trait_type
        String trait_unit
        String trait_transformation
        String trait_definition
        String covariates
        String reference_assembly
        Int min_MAC_filter
        String genotyping_technology
        String genotyping_platform
        Boolean is_imputed
        String? imputation_reference_panel
        String? imputation_reference_panel_detail
        Float? imputation_quality_filter
        String cohorts
        String population_descriptor
        String population_labels
        String population_proportions
        String countries_of_recruitment
        String? countries_of_birth
        String analysis_method
        String analysis_software
        File pheno_file
        String age_column
        String sex_column
        String female_value
    }

    String covariate_string = sub(covariates, ",", "|")
    String tech_string = sub(genotyping_technology, ",", "|")
    String platform_string = sub(genotyping_platform, ",", "|")
    String cohort_string = sub(cohorts, ",", "|")
    String pop_label_string = sub(population_labels, ",", "|")
    String pop_prop_string = sub(population_proportions, ",", "|")
    String recruitment_string = sub(countries_of_recruitment, ",", "|")
    String birth_string = sub(select_first([countries_of_birth, "NA"]), ",", "|")

    command <<<
        Rscript -e "\
        library(dplyr); \
        library(readr); \
        parse_array <- function(x) unlist(strsplit(x, split=' ', fixed=TRUE)); \
        file_path <- parse_array('~{sep=' ' data_files}'); \
        md5sum <- parse_array('~{sep=' ' md5sum}'); \
        chromosome <- parse_array('~{sep=' ' chromosome}'); \
        n_variants <- parse_array('~{sep=' ' n_variants}'); \
        total_variants <- n_variants[grep('all_variants', basename(file_path))]; \
        dat <- tibble(file_path=file_path, \
            md5sum=md5sum, \
            chromosome=chromosome, \
            n_variants=n_variants, \
            file_type='data'); \
        write_tsv(dat, 'gsr_file_table.tsv'); \
        analysis <- list(gsr_source='PRIMED', \
            consent_code='~{consent_code}', \
            upload_date=as.character(Sys.Date()), \
            contributor_contact='~{contributor_contact}', \
            trait='~{trait}', \
            trait_type='~{trait_type}', \
            trait_unit='~{trait_unit}', \
            trait_transformation='~{trait_transformation}', \
            trait_definition='~{trait_definition}', \
            covariates='~{covariate_string}', \
            reference_assembly='~{reference_assembly}', \
            n_variants=total_variants, \
            min_MAC_filter='~{min_MAC_filter}', \
            genotyping_technology='~{tech_string}', \
            genotyping_platform='~{platform_string}', \
            is_imputed='~{true='TRUE' false='FALSE' is_imputed}', \
            imputation_reference_panel='~{default='NA' imputation_reference_panel}', \
            imputation_reference_panel_detail='~{default='NA' imputation_reference_panel_detail}', \
            imputation_quality_filter='~{default='NA' imputation_quality_filter}', \
            is_meta_analysis='FALSE', \
            analysis_method='~{analysis_method}', \
            analysis_software='~{analysis_software}', \
            cohorts='~{cohorts}', \
            population_descriptor='~{population_descriptor}', \
            population_labels='~{pop_label_string}', \
            population_proportions='~{pop_prop_string}', \
            countries_of_recruitment='~{recruitment_string}', \
            countries_of_birth='~{birth_string}'
            ); \
        phen <- read_csv('~{pheno_file}'); \
        n_samp <- nrow(phen); \
        if ('~{trait_type}' == 'binary') { \
            n_case <- sum(na.omit(phen[['~{trait}']] == 1)); \
            n_ctrl <- sum(na.omit(phen[['~{trait}']] == 0)); \
            n_eff <- 4/(1/n_case + 1/n_ctrl); \
            analysis <- c(analysis, list( \
                n_samp=n_samp, \
                n_case=n_case, \
                n_ctrl=n_ctrl, \
                n_effective=n_eff) \
                ); \
        } else { \
            analysis <- c(analysis, list( \
                n_samp=n_samp, \
                n_effective=n_samp) \
                ); \
        }; \
        if (is.element('~{sex_column}', names(phen))) { \
            prop_f <- sum(phen[['~{sex_column}']] == '~{female_value}') / n_samp; \
            analysis <- c(analysis, list( \
                proportion_female=prop_f) \
                ); \
        }; \
        if (is.element('~{age_column}', names(phen))) { \
            age <- na.omit(phen[['~{age_column}']]); \
            analysis <- c(analysis, list( \
                age_mean=mean(age), \
                age_median=median(age), \
                age_min=min(age), \
                age_max=max(age)) \
                ); \
            if ('~{trait_type}' == 'binary') { \
                age_case <- na.omit(phen[['~{age_column}']][phen[['~{trait}']] == 1]); \
                age_ctrl <- na.omit(phen[['~{age_column}']][phen[['~{trait}']] == 0]); \
                analysis <- c(analysis, list( \
                    age_mean_case=mean(age_case), \
                    age_median_case=median(age_case), \
                    age_min_case=min(age_case), \
                    age_max_case=max(age_case), \
                    age_mean_ctrl=mean(age_ctrl), \
                    age_median_ctrl=median(age_ctrl), \
                    age_min_ctrl=min(age_ctrl), \
                    age_max_ctrl=max(age_ctrl)) \
                ); \
            }; \
        }; \
        analysis_vec <- unlist(analysis)
        analysis_vec <- analysis_vec[analysis_vec != 'NA']
        write_tsv(tibble(field=names(analysis_vec), value=analysis_vec), 'analysis_table.tsv'); \
        "
    >>>

    output {
        Map[String, File] table_files = {
            "analysis": "analysis_table.tsv",
            "gsr_file": "gsr_file_table.tsv"
        }
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-gcr-public/anvil-rstudio-bioconductor:3.16.0"
    }
}
