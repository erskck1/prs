#!/usr/bin/env Rscript

# Script to take base GWAS catalogue summary statistics and format them so that they match SAIGE statistics and can be used for input to PRSice input:

###################
# Import packages #
###################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})

calculate_maf <- function(eaf) {
  if(!is.na(eaf)) {
    if(eaf < 0.5) {
      return(eaf)
    } else {
      return(1-eaf)
    }
  } else {
    return(eaf)
  }
}

#####################
# Parsing arguments #
#####################

option_list <- list(make_option(c("--input_gwas_cat"), action="store", type='character',help="String containing input GWAS catalogue file (base cohort)."))

args = parse_args(OptionParser(option_list = option_list))

# Arg to variable
input_gwas_cat = '/home/ersoykocak/Downloads/27989323-GCST004420-EFO_0008082.h.tsv.gz'#args$input_gwas_cat

######################################
# Importing data and transforming it #
######################################

gwas_catalogue_file <- as_tibble(fread(input_gwas_cat))

#### Keep harmonized data only ####

base <- gwas_catalogue_file %>% select(starts_with("hm_"), starts_with("p_"),starts_with("effect_allele_"))

#### Remove SNPs with no beta or OR - these cannot be used by PRSice ####

base <- filter(base, !(is.na(hm_beta) == TRUE & is.na(hm_odds_ratio) == TRUE))

#### For SNPs with at least a beta or OR, alternatively use the beta or OR to calculate the other ####

base$hm_beta <- as.numeric(base$hm_beta)
base$hm_odds_ratio <- as.numeric(base$hm_odds_ratio)

base <- base %>%
  mutate(hm_beta = if_else(is.na(hm_beta), log(hm_odds_ratio), hm_beta), 
         hm_odds_ratio = if_else(is.na(hm_odds_ratio), exp(hm_beta), hm_odds_ratio))

#### For SNPs with no p-value, replace the NA by a 1 ####

base <- base %>% mutate(p_value = if_else(is.na(p_value), 1, p_value))

#### Remove duplicate SNPs - these cannot be used by PRSice (an error will be thrown) ####

base <- distinct(base, hm_rsid, .keep_all = TRUE)
base <- distinct(base, hm_variant_id, .keep_all = TRUE)

#### Calculate MAF using effect allele frequency ####

base$MAF <- sapply(base$effect_allele_frequency, calculate_maf)
number_of_na <-sum(is.na(base$effect_allele_frequency)) 

if(number_of_na/nrow(base) > 0.25) {
  base <- subset(base, is.na(MAF) | MAF > 0.01)
} else {
  base <- subset(base, !is.na(MAF) & MAF > 0.01)
}

#### Change column names to match SAIGE output ####

colnames(base)[which(names(base) == "hm_rsid")] <- "SNPID"
colnames(base)[which(names(base) == "hm_chrom")] <- "CHR"
colnames(base)[which(names(base) == "hm_pos")] <- "POS"
colnames(base)[which(names(base) == "hm_other_allele")] <- "Allele2"
colnames(base)[which(names(base) == "hm_effect_allele")] <- "Allele1"
colnames(base)[which(names(base) == "hm_beta")] <- "BETA"
colnames(base)[which(names(base) == "hm_odds_ratio")] <- "OR"
colnames(base)[which(names(base) == "p_value")] <- "p.value"

#### Save data ####

write.table(base, "base.data", quote = F, row.names =F, sep = " ")


