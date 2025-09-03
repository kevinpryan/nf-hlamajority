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

# take in command line arguments
option_list = list(
  make_option(c("-d", "--dataset"), type="character", default=NULL, 
              help="name of dataset", metavar="character"),
  make_option(c("-p", "--pat_name"), type="character", default=NULL, 
              help="patient ID", metavar="character"),
  make_option(c("-r", "--run"), type="character", default=NULL,
              help="LENS run", metavar="character"),
  make_option(c("-cnd", "--all_calls_nd"), type="character", default=NULL,
              help="path to all calls normal DNA", metavar="character"),
  make_option(c("-cad", "--all_calls_ad"), type="character", default=NULL,
              help="path to all calls tumour DNA", metavar="character"),
  make_option(c("-mnd", "--majority_voting_nd"), type="character", default=NULL,
              help="path to majority vote normal DNA output", metavar="character"),
  make_option(c("-mad", "--majority_voting_ad"), type="character", default=NULL,
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

dataset <- "hlamajority"
pat_name <- "Pt01"
run <- "ad-1-nd2"
all_calls_nd <- "testdata/Pt01_nd-2_hlamajority_all_calls_mhci.tsv"
majority_voting_nd <- "testdata/Pt01_nd-2_hlamajority_majority_vote_mhci.tsv"
all_calls_ad <- "testdata/Pt01_ad-1_hlamajority_all_calls_mhci.tsv"
majority_voting_ad <- "testdata/Pt01_ad-1_hlamajority_majority_vote_mhci.tsv"

all_calls_nd_in <- vroom(all_calls_nd, col_types = "cccccccc")
majority_voting_nd_in <- vroom(majority_voting_nd, col_types = "ccccccc")
all_calls_ad_in <- vroom(all_calls_ad, col_types = "cccccccc")
majority_voting_ad_in <- vroom(majority_voting_ad,  col_types = "cccccccc")

add_hla_prefix <- function(df, cols){
  df_new <- df
  for (i in 1:length(cols)){
    gene <- substr(cols[i], 1, 1)
    df_new[[cols[i]]] <- ifelse(df_new[[cols[i]]] == "NA", NA, paste("HLA-", gene, df_new[[cols[i]]], sep = ""))
  }
  return(df_new)
}
hla_cols <- c("A1", "A2", "B1", "B2", "C1", "C2")
all_calls_nd_in_hla <- add_hla_prefix(all_calls_nd_in, cols = hla_cols)
majority_voting_nd_in_hla <- add_hla_prefix(majority_voting_nd_in, cols = hla_cols)
all_calls_ad_in_hla <- add_hla_prefix(all_calls_ad_in, cols = hla_cols)
majority_voting_ad_in_hla <- add_hla_prefix(majority_voting_ad_in, cols = hla_cols)

hla_alleles_support <- data.frame(alleles = unique(unlist(majority_voting_ad_in_hla[,hla_cols]), unlist(majority_voting_nd_in_hla[,hla_cols])))
