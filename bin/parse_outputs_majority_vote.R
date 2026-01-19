#!/usr/bin/env Rscript

# script to take in the outputs of 4 HLA typing tools and carry out majority voting to get the MHC I calls
# usage: Rscript parse_outputs_majority_vote.R --samplename meta.id --optitype /path/to/optitype/ --polysolver /path/to/polysolver/ --kourami /path/to/kourami/ --hlala /path/to/hlala/ --benchmark /path/to/benchmark/results.txt
# load libraries
library(optparse)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(vroom)
# source functions
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/HLA-LA_conversion.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/Optitype_conversion.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/Polysolver_conversion.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/kourami_conversion.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/majority_voting.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/df_to_list.R")
# source("~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/bin/are_vectors_identical.R")
box::purge_cache()
box::use(lib/HLA_LA_conversion[...])
box::use(lib/Optitype_conversion[...])
box::use(lib/Polysolver_conversion[...])
box::use(lib/kourami_conversion[...])
box::use(lib/majority_voting[...])
box::use(lib/df_to_list[...])
box::use(lib/are_vectors_identical[...])
# take in command line arguments

option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL, 
              help="name of sample", metavar="character"),
  make_option(c("-o", "--optitype"), type="character", default=NULL, 
              help="path to optitype output", metavar="character"),
  make_option(c("-p", "--polysolver"), type="character", default=NULL,
              help="path to polysolver output", metavar="character"),
  make_option(c("-l", "--hlala"), type="character", default=NULL,
              help="path to hlala output", metavar="character"),
  make_option(c("-k", "--kourami"), type="character", default=NULL,
              help="path to kourami output", metavar="character"),
  make_option(c("-b", "--benchmark"), type="character", default=NULL,
             help="path to benchmark rankings", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 6){
  print_help(opt_parser)
  stop("6 arguments must be supplied", call.=FALSE)
}

# assign arguments to variables
opt$samplename -> samplename
opt$optitype -> optitype_in
opt$polysolver -> polysolver_in
opt$hlala -> hlala_in
opt$kourami -> kourami_in
opt$benchmark -> benchmark_in

# Read in files

#samplename <- "3532"
#optitype_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/test_outputs/optitype_calls/3532/optitype_calls/"
#polysolver_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/test_outputs/polysolver_calls/3532/polysolver_calls/"
#hlala_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/test_outputs/hlala_calls/3532/hlala_calls/"
#kourami_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/test_outputs/kourami_calls/"
#benchmark_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/assets/benchmarking_results_claeys.csv"
#samplename <- "3532"
#optitype_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/test_outputs/optitype_calls/3532/optitype_calls/"
#polysolver_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/test_outputs/polysolver_calls/3532/polysolver_calls/"
#hlala_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/test_outputs/hlala_calls/3532/hlala_calls/"
#kourami_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/test_outputs/kourami_calls/"
#benchmark_in <- "~/Documents/PhD/misc/useful/nextflow/nf-hlatyping/assets/benchmarking_results_claeys.csv"
#samplename <- "3532"
#print("trying to find hlala output dir")
#fileList <- list.files(path = outputFolder, pattern = "*_bestguess_G.txt$")
print("reading in hlala")
hlala <- toolOutputToR.HLA_LA(hlala_in, mhci_only = T, trim = T)
print(hlala)
print("reading in optitype...")
optitype <- toolOutputToR.Optitype(optitype_in)
print(optitype)
print("reading in polysolver...")
polysolver <- toolOutputToR.Polysolver(polysolver_in, trim = T)
print(polysolver)
print("reading in kourami...")
kourami <- toolOutputToR.kourami(kourami_in, mhci_only = T, trim = T)
print(kourami)
print("combining...")
combined <- rbind(hlala, optitype, polysolver, kourami)
rownames(combined) <- c("hlala", "optitype", "polysolver", "kourami")
combined$tool <- rownames(combined)
combined$sample <- rep(samplename, nrow(combined))
print(combined)
# Read in benchmarking (might not be necessary - just using optitype as best)
benchmark <- read.csv(benchmark_in)
benchmark <- benchmark %>% dplyr::filter(tool %in% c("HLA*LA", "Kourami", "Optitype", "Polysolver") & seq_type == "WES")
benchmark$A <- rev(rank(benchmark$A))
benchmark$B <- rev(rank(benchmark$B))
benchmark$C <- rev(rank(benchmark$C))
benchmark <- benchmark %>% dplyr::select(c(tool, A, B, C))
benchmark$tool <- c("hlala", "kourami", "optitype", "polysolver")

# rename cols
colnames(combined) <- c("A1", "A2", "B1", "B2", "C1", "C2", "tool", "sample")

### the commented out code tests out different different scenarios
#calls_all_identical <- list(optitype = "03:01", polysolver = "03:01",  kourami = "03:01", hlala = "03:01")
#calls_all_identical <- data.frame(A1 = rep("03:01",4), A2 = rep("05:01", 4), B1 = rep("01:01", 4), B2 = rep("02:02",4), C1 = rep("03:03",4), C2 = rep("06:01",4))
#rownames(calls_all_identical) <- c("hlala", "kourami", "optitype", "polysolver")
#A_list_test <- df_to_list(calls_all_identical, cols = c("A1", "A2"))
#A_list_test_notna <- A_list_test[not_na(A_list_test)]
#A_test_identical <- outer(A_list_test_notna, A_list_test_notna, FUN = are_vectors_identical_vectorised)
#A_test_identical_vote <- majority_vote_comparison(A_test_identical, A_list_test_notna, benchmark, "A")

# calls_none_identical <- list(optitype = "03:01", polysolver = "04:01",  kourami = "05:01", hlala = "06:01")
# calls_optitype_na <- list(optitype = NA, polysolver = "04:01",  kourami = "05:01", hlala = "06:01")
# calls_kourami_na <- list(optitype = "02:01", polysolver = "04:01",  kourami = NA, hlala = "04:01")
# calls_tie <- list(optitype = "02:01", polysolver = "02:01",  kourami = "05:01", hlala = "05:01")
# calls_tie_comp <- outer(calls_tie, calls_tie, FUN = are_vectors_identical_vectorised)
# majority_vote_comparison(calls_tie_comp, calls_tie_comp, benchmark, "A")
# 
# calls_kourami_na <- calls_kourami_na[not_na(calls_kourami_na)]
# calls_kourami_na_comp <- outer(calls_kourami_na, calls_kourami_na, FUN = are_vectors_identical_vectorised)
# majority_vote_comparison(calls_kourami_na_comp, calls_kourami_na_comp, benchmark, "A")
# 
# calls_optitype_na <- calls_optitype_na[not_na(calls_optitype_na)]
# calls_optitype_na_comp <- outer(calls_optitype_na, calls_optitype_na, FUN = are_vectors_identical_vectorised)
# majority_vote_comparison(calls_kourami_na_comp, calls_optitype_na, benchmark, "A")
# 
# call_only_kourami <- list(optitype = NA, polysolver = NA,  kourami = "06:01", hlala = NA)
# call_only_kourami <- call_only_kourami[not_na(call_only_kourami)]
# call_only_kourami_comp <- outer(call_only_kourami, call_only_kourami, FUN = are_vectors_identical_vectorised)
# majority_vote_comparison(call_only_kourami_comp, call_only_kourami, benchmark, "A")
# 
# call_two_tools_identical <- list(optitype = c(NA,NA), polysolver = c(NA,NA),  kourami = c("06:01","06:01"), hlala = c("06:01", "06:01"))
# call_two_tools_identical <- call_two_tools_identical[not_na(call_two_tools_identical)]
# call_two_tools_identical_comp <- outer(call_two_tools_identical, call_two_tools_identical, FUN = are_vectors_identical_vectorised)
# outfile <- majority_vote_comparison(call_two_tools_identical_comp, call_two_tools_identical, benchmark, "A")

all_na <- list(optitype = c(NA,NA), polysolver = c(NA,NA),  kourami = c(NA,NA), hlala = c(NA, NA))
all_na <- all_na[not_na(all_na)]
all_na_identical <- outer(all_na, all_na, FUN = are_vectors_identical_vectorised)
outfile <- majority_vote_comparison(comparison = all_na_identical, calls = all_na, benchmark = benchmark, hla = "A")

# Run majority voting for HLA-A
A_list <- df_to_list(combined, cols = c("A1", "A2"))
print("A_list")
print(A_list)
A_list_notna <- A_list[not_na(A_list)]
print("A_list_notna")
print(A_list_notna)
A_identical <- outer(A_list_notna, A_list_notna, FUN = are_vectors_identical_vectorised)
print("A_identical")
print(A_identical)
A_vote <- majority_vote_comparison(A_identical, A_list, benchmark, "A")
print("A_vote")
print(A_vote)
# Run majority voting for HLA-B
B_list <- df_to_list(combined, cols = c("B1", "B2"))
B_list_notna <- B_list[not_na(B_list)]
B_identical <- outer(B_list_notna, B_list_notna, FUN = are_vectors_identical_vectorised)
B_vote <- majority_vote_comparison(B_identical, B_list, benchmark, "B")

# Run majority voting for HLA-C
C_list <- df_to_list(combined, cols = c("C1", "C2"))
C_list_notna <- C_list[not_na(C_list)]
C_identical <- outer(C_list_notna, C_list_notna, FUN = are_vectors_identical_vectorised)
C_vote <- majority_vote_comparison(C_identical, C_list, benchmark, "C")

# prepare outputs: table with all calls across all tools, table with majority voting result
rownames(combined) <- NULL
full_output <- combined %>% relocate(., sample, .before = A1) %>% relocate(., tool, .before = A1)
write.table(full_output, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_all_calls_mhci.tsv", sep = ""))
dummy_output <- data.frame(x = c(1,2), y = c(3,4))
write.table(dummy_output, file = "dummy_out8.txt")

majority_output <- data.frame(sample = samplename)
majority_output$A1 <- A_vote["A1"]
majority_output$A2 <- A_vote["A2"]
majority_output$B1 <- B_vote["B1"]
majority_output$B2 <- B_vote["B2"]
majority_output$C1 <- C_vote["C1"]
majority_output$C2 <- C_vote["C2"]
print("majority vote output")
print(majority_output)
write.table(majority_output, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_majority_vote_mhci.tsv", sep = ""))
