#!/usr/bin/env Rscript

# script to take in the majority voting results and all the calls from the normal and tumour DNA samples from LENS and produce the output files for the HLAtyping workflow:
# ${dataset}-${pat_name}-${run}.hla_alleles_support.tsv and ${dataset}-${pat_name}-${run}.supported_hla_alleles
# ${dataset}-${pat_name}-${run}.hla_alleles_support.tsv looks like the following:
# allele	hla_allele_support
# HLA-A02:01	optitype:tumor_dna;seq2hla:tumor_rna;optitype:normal_dna
# HLA-B44:02	optitype:tumor_dna;seq2hla:tumor_rna;optitype:normal_dna
# HLA-C05:01	optitype:tumor_dna;seq2hla:tumor_rna;optitype:normal_dna
# ${dataset}-${pat_name}-${run}.supported_hla_alleles looks like the following:
# HLA-A02:01,HLA-B44:02,HLA-C05:01
# usage: Rscript parse_outputs_majority_vote.R --dataset ${dataset} --pat_name ${pat_name} --run ${run} --all_calls_ad /path/to/all_calls_ad.tsv --all_calls_nd /path/to/all_calls_nd.tsv --majority_voting_ad /path/to/majority_voting_ad.tsv --majority_voting_nd /path/to/majority_voting_nd.tsv 
# load libraries
library(optparse)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(vroom)
library(tidytable)

# take in command line arguments
option_list = list(
  make_option(c("--dataset"), type="character", default=NULL, 
              help="name of dataset", metavar="character"),
  make_option(c("--pat_name"), type="character", default=NULL, 
              help="patient ID", metavar="character"),
  make_option(c("--run"), type="character", default=NULL,
              help="LENS run", metavar="character"),
  make_option(c("--all_calls_nd"), type="character", default=NULL,
              help="path to all calls normal DNA", metavar="character"),
  make_option(c("--all_calls_ad"), type="character", default=NULL,
              help="path to all calls tumour DNA", metavar="character"),
  make_option(c("--majority_voting_nd"), type="character", default=NULL,
              help="path to majority vote normal DNA output", metavar="character"),
  make_option(c("--majority_voting_ad"), type="character", default=NULL,
              help="path to majority vote tumour DNA output", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 7){
  print_help(opt_parser)
  stop("7 arguments must be supplied", call.=FALSE)
}

# assign arguments to variables
opt$dataset -> dataset
opt$pat_name -> pat_name
opt$run -> run
opt$all_calls_nd -> all_calls_nd
opt$all_calls_ad -> all_calls_ad
opt$majority_voting_nd -> majority_voting_nd
opt$majority_voting_ad -> majority_voting_ad

add_hla_prefix <- function(df, cols){
  df_new <- df
  for (i in 1:length(cols)){
    gene <- substr(cols[i], 1, 1)
    df_new[[cols[i]]] <- ifelse(is.na(df_new[[cols[i]]]), 
                                NA, 
                                paste0("HLA-", gene, df_new[[cols[i]]]))
  }
  return(df_new)
}

create_hla_alleles_support_table <- function(nd_df, ad_df){
  all_alleles_nd <- unique(sort(unlist(nd_df[,3:8])))
  all_alleles_ad <- unique(sort(unlist(ad_df[,3:8])))
  all_alleles <- unique(sort(c(all_alleles_nd, all_alleles_ad)))
  long_format_calls_nd <- nd_df %>%
    pivot_longer(
      cols = starts_with(c("A", "B", "C")), # Select all allele columns
      names_to = "locus",                   # New column for the original column names (A1, A2, etc.)
      values_to = "allele",            # New column for the HLA types
      values_drop_na = TRUE                 # Automatically remove rows where the allele call was NA
    )
  result_nd <- long_format_calls_nd %>%
    # Keep only the rows where the allele_call is in your 'alleles' vector
    filter(allele %in% all_alleles) %>%
    # Group by the allele call
    group_by(allele) %>%
    # Summarize to get a unique list of tools for each allele
    summarise(tools = list(unique(tool)))
  final_result_nd <- result_nd %>%
    mutate(hla_allele_support = map_chr(tools, ~ paste(., collapse = ":normal_dna;"))) %>% 
    dplyr::select(-tools) %>% 
    mutate(hla_allele_support_normal = paste(hla_allele_support, ":normal_dna", sep = "")) %>% 
    dplyr::select(-hla_allele_support)
  
  long_format_calls_ad <- ad_df %>%
    pivot_longer(
      cols = starts_with(c("A", "B", "C")), # Select all allele columns
      names_to = "locus",                   # New column for the original column names (A1, A2, etc.)
      values_to = "allele",            # New column for the HLA types
      values_drop_na = TRUE                 # Automatically remove rows where the allele call was NA
    )
  result_ad <- long_format_calls_ad %>%
    # Keep only the rows where the allele_call is in your 'alleles' vector
    filter(allele %in% all_alleles) %>%
    # Group by the allele call
    group_by(allele) %>%
    # Summarize to get a unique list of tools for each allele
    summarise(tools = list(unique(tool)))
  final_result_ad <- result_ad %>%
    mutate(hla_allele_support = map_chr(tools, ~ paste(., collapse = ":tumor_dna;"))) %>% 
    dplyr::select(-tools) %>% 
    mutate(hla_allele_support_tumor = paste(hla_allele_support, ":tumor_dna", sep = "")) %>% 
    dplyr::select(-hla_allele_support)
  final_result <- full_join(final_result_nd, final_result_ad, by = "allele") %>% 
    mutate(hla_allele_support_tools = paste(hla_allele_support_normal, hla_allele_support_tumor, sep = ";")) %>% 
    dplyr::select(-c(hla_allele_support_normal, hla_allele_support_tumor))
  return(final_result)
}

create_hla_alleles_support_table_majority_vote <- function(majority_voting_table_normal, majority_voting_table_tumor){
  support_table_normal <- data.frame(allele = unlist(unique(majority_voting_table_normal[,2:7])), hla_allele_support_normal = "majority_vote_normal_dna")
  support_table_tumor <- data.frame(allele = unique(unlist(majority_voting_table_tumor[,2:7])), hla_allele_support_tumor = "majority_vote_tumor_dna")
  support_table <- full_join(support_table_normal, support_table_tumor) %>% 
    mutate(hla_allele_support_majority_vote = paste(hla_allele_support_normal, hla_allele_support_tumor, sep = ";")) %>% 
    dplyr::select(-c(hla_allele_support_normal, hla_allele_support_tumor))
  return(support_table)
}

# These lines are used to test out the code
# dataset <- "hlamajority"
# pat_name <- "Pt01"
# run <- "ad-1-nd2"
# all_calls_nd <- "testdata/Pt01_nd-2_hlamajority_all_calls_mhci.tsv"
# majority_voting_nd <- "testdata/Pt01_nd-2_hlamajority_majority_vote_mhci.tsv"
# all_calls_ad <- "testdata/Pt01_ad-1_hlamajority_all_calls_mhci.tsv"
# majority_voting_ad <- "testdata/Pt01_ad-1_hlamajority_majority_vote_mhci.tsv"

# read in data
all_calls_nd_in <- vroom(all_calls_nd, col_types = "cccccccc")
majority_voting_nd_in <- vroom(majority_voting_nd, col_types = "ccccccc")
all_calls_ad_in <- vroom(all_calls_ad, col_types = "cccccccc")
majority_voting_ad_in <- vroom(majority_voting_ad,  col_types = "cccccccc")
# add HLA- prefix to all data
hla_cols <- c("A1", "A2", "B1", "B2", "C1", "C2")
all_calls_nd_in_hla <- add_hla_prefix(all_calls_nd_in, cols = hla_cols)
majority_voting_nd_in_hla <- add_hla_prefix(majority_voting_nd_in, cols = hla_cols)
all_calls_ad_in_hla <- add_hla_prefix(all_calls_ad_in, cols = hla_cols)
majority_voting_ad_in_hla <- add_hla_prefix(majority_voting_ad_in, cols = hla_cols)

# create table with HLA support across all tools and tumor/normal DNA
alleles_support_table_tools <- create_hla_alleles_support_table(all_calls_nd_in_hla, all_calls_ad_in_hla)
# create table telling whether each allele was chosen by majority voting using tumor and/or normal DNA
alleles_support_table_majority_vote <- create_hla_alleles_support_table_majority_vote(majority_voting_nd_in_hla, majority_voting_ad_in_hla)
# join these two tables so that they have a column hla_allele_support telling the user which tools called the allele and whether they were included in the majority voting output
alleles_support_table <- full_join(alleles_support_table_majority_vote, alleles_support_table_tools) %>% 
                         mutate(hla_allele_support = paste(hla_allele_support_majority_vote, hla_allele_support_tools, sep= ";")) %>% 
                        dplyr::select(-c(hla_allele_support_majority_vote, hla_allele_support_tools))                                                                                
alleles_support_table$hla_allele_support <- gsub(";?NA;?", "", alleles_support_table$hla_allele_support)

# create output file for majority voting
majority_votes_all <- sort(unique(c(unlist(majority_voting_nd_in_hla[,2:7]), unlist(majority_voting_ad_in_hla[,2:7]))))
majority_votes_df <- t(data.frame(majority_votes_all))
majority_votes_df_outfile <- paste(dataset, pat_name, run, "supported_hla_alleles.csv", sep = "-")
write.table(majority_votes_df, file = majority_votes_df_outfile, row.names = F, col.names = F, quote = F, sep = ",")                           
alleles_support_table_outfile <- paste(dataset, pat_name, run, "hla_alleles_support.csv", sep = "-")
write.table(alleles_support_table, file = alleles_support_table_outfile, sep = "\t", quote = F, row.names = F)
