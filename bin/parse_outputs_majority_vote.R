#!/usr/bin/env Rscript

# script to take in the outputs of 4 HLA typing tools and carry out majority voting to get the MHC I calls
# usage: Rscript parse_outputs_majority_vote.R --samplename meta.id --optitype /path/to/optitype/ --polysolver /path/to/polysolver/ --kourami /path/to/kourami/ --hlala /path/to/hlala/ --weights /path/to/weights/results.csv --method "majority|weighted" --depth /path/to/depth/file

# load libraries
library(optparse)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)
library(vroom)
box::purge_cache()
box::use(lib/HLA_LA_conversion[...])
box::use(lib/Optitype_conversion[...])
box::use(lib/Polysolver_conversion[...])
box::use(lib/kourami_conversion[...])
box::use(lib/majority_voting[...])
box::use(lib/df_to_list[...])
box::use(lib/are_vectors_identical[...])
box::use(lib/extract_tools_contributing_to_vote[...])
box::use(lib/extract_na_tools[...])
box::use(lib/weighted_vote[...])

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
  make_option(c("-w", "--weights"), type="character", default=NULL,
              help="path to weights - converts to rankings if using default majority vote", metavar="character"),
  make_option(c("-m", "--method"), type="character", default="majority",
              help="use default majority vote (majority) or weighted voting (weighted)"),
  make_option(c("-d", "--depth"), type="character", default=NULL,
              help="Path to mean depth file for HLA exons 2 and 3", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 7){
  print_help(opt_parser)
  stop("Seven arguments must be supplied", call.=FALSE)
}

# assign arguments to variables
opt$samplename -> samplename
opt$optitype -> optitype_in
opt$polysolver -> polysolver_in
opt$hlala -> hlala_in
opt$kourami -> kourami_in
opt$weights -> weights_in
opt$method -> method
opt$depth -> depth

message(paste("method: ", method, sep = ""))

# Read in files

# samplename <- "HG00097"
# optitype_in <- "../assets/test-outputs/test-outputs-1000genomes/HG00097/optitype_calls/"
# polysolver_in <- "../assets/test-outputs/test-outputs-1000genomes/HG00097/polysolver_calls/HG00097/polysolver_calls/"
# hlala_in <- "../assets/test-outputs/test-outputs-1000genomes/HG00097/hlala_calls/"
# kourami_in <- "../assets/test-outputs/test-outputs-1000genomes/HG00097/kourami_calls/"
# benchmark_in <- "../assets/benchmarking_results_claeys_cleaned.csv"
# depth <- "../assets/test-outputs/test-outputs-1000genomes/HG00097/mosdepth/HG00097-mosdepth-summary.tsv"

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
# read in depth
depth_df <- read.delim(depth, stringsAsFactors = FALSE)
print("combining...")
combined <- rbind(hlala, optitype, polysolver, kourami)

rownames(combined) <- c("hlala", "optitype", "polysolver", "kourami")
combined$tool <- rownames(combined)
combined$sample <- rep(samplename, nrow(combined))
colnames(combined) <- c("A1", "A2", "B1", "B2", "C1", "C2", "tool", "sample")
print(combined)

# convert to list format
A_list <- df_to_list(combined, cols = c("A1", "A2"))
B_list <- df_to_list(combined, cols = c("B1", "B2"))
C_list <- df_to_list(combined, cols = c("C1", "C2"))

# remove NAs
A_list_notna <- A_list[not_na(A_list)]
B_list_notna <- B_list[not_na(B_list)]
C_list_notna <- C_list[not_na(C_list)]

# Read in benchmarking
weights <- read.csv(weights_in)

if (opt$method == "majority") {
  print("using unweighted majority voting")
  weights$A <- rev(rank(weights$A))
  weights$B <- rev(rank(weights$B))
  weights$C <- rev(rank(weights$C))
  
  ### the commented out code tests out different different scenarios
  # calls_all_identical <- list(optitype = "03:01", polysolver = "03:01",  kourami = "03:01", hlala = "03:01")
  # calls_all_identical <- data.frame(A1 = rep("03:01",4), A2 = rep("05:01", 4), B1 = rep("01:01", 4), B2 = rep("02:02",4), C1 = rep("03:03",4), C2 = rep("06:01",4))
  # rownames(calls_all_identical) <- c("hlala", "kourami", "optitype", "polysolver")
  # A_list_test <- df_to_list(calls_all_identical, cols = c("A1", "A2"))
  # A_list_test_notna <- A_list_test[not_na(A_list_test)]
  # A_test_identical <- outer(A_list_test_notna, A_list_test_notna, FUN = are_vectors_identical_vectorised)
  # A_test_identical_vote <- majority_vote_comparison(A_test_identical, A_list_test_notna, benchmark, "A")
  
  
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
  
  # all_na <- list(optitype = c(NA,NA), polysolver = c(NA,NA),  kourami = c(NA,NA), hlala = c(NA, NA))
  # all_na <- all_na[not_na(all_na)]
  # all_na_identical <- outer(all_na, all_na, FUN = are_vectors_identical_vectorised)
  # outfile <- majority_vote_comparison(comparison = all_na_identical, calls = all_na, benchmark = benchmark, hla = "A")
  # matching_tools <- extract_tools_contributing_to_vote(genotype_list = all_na, genotype_call_vector = outfile)
  # n_tools_called_A <- length(all_na)
  # n_tools_support_A <- length(strsplit(matching_tools, ",")[[1]])
  # support_A <- n_tools_support_A / n_tools_called_A
  
  # Run majority voting for HLA-A
  print("A_list")
  print(A_list)
  print("A_list_notna")
  print(A_list_notna)
  A_identical <- outer(A_list_notna, A_list_notna, FUN = are_vectors_identical_vectorised)
  print("A_identical")
  print(A_identical)
  A_vote <- majority_vote_comparison(A_identical, A_list, weights, "A")
  print("A_vote")
  print(A_vote)
  matching_tools_A <- extract_tools_contributing_to_vote(genotype_list = A_list, genotype_call_vector = A_vote)
  n_tools_called_A <- length(A_list_notna)
  n_tools_support_A <- length(strsplit(matching_tools_A, ",")[[1]])
  support_A <- round(n_tools_support_A / n_tools_called_A, 2)
  
  A_df <- data.frame(
    sample = samplename,
    gene = "HLA-A",
    allele1 = A_vote[[1]],
    allele2 = A_vote[[2]],
    support = support_A,
    matching_tools = matching_tools_A,
    method = "majority_vote",
    weight_winner = n_tools_support_A,
    total_weight = n_tools_called_A,
    n_tools_support = n_tools_support_A,
    n_tools_called = n_tools_called_A
  )
  
  # Run majority voting for HLA-B
  B_identical <- outer(B_list_notna, B_list_notna, FUN = are_vectors_identical_vectorised)
  B_vote <- majority_vote_comparison(B_identical, B_list, weights, "B")
  matching_tools_B <- extract_tools_contributing_to_vote(genotype_list = B_list, genotype_call_vector = B_vote)
  n_tools_called_B <- length(B_list_notna)
  n_tools_support_B <- length(strsplit(matching_tools_B, ",")[[1]])
  support_B <- round(n_tools_support_B / n_tools_called_B, 2)
  
  B_df <- data.frame(
    sample = samplename,
    gene = "HLA-B",
    allele1 = B_vote[[1]],
    allele2 = B_vote[[2]],
    support = support_B,
    matching_tools = matching_tools_B,
    method = "majority_vote",
    weight_winner = n_tools_support_B,
    total_weight = n_tools_called_B,
    n_tools_support = n_tools_support_B,
    n_tools_called = n_tools_called_B
  )
  
  # Run majority voting for HLA-C
  C_identical <- outer(C_list_notna, C_list_notna, FUN = are_vectors_identical_vectorised)
  C_vote <- majority_vote_comparison(C_identical, C_list, weights, "C")
  matching_tools_C <- extract_tools_contributing_to_vote(genotype_list = C_list, genotype_call_vector = C_vote)
  n_tools_called_C <- length(C_list_notna)
  n_tools_support_C <- length(strsplit(matching_tools_C, ",")[[1]])
  support_C <- round(n_tools_support_C / n_tools_called_C, 2)
  
  C_df <- data.frame(
    sample = samplename,
    gene = "HLA-C",
    allele1 = C_vote[[1]],
    allele2 = C_vote[[2]],
    support = support_C,
    matching_tools = matching_tools_C,
    method = "majority_vote",
    weight_winner = n_tools_support_C,
    total_weight = n_tools_called_C,
    n_tools_support = n_tools_support_C,
    n_tools_called = n_tools_called_C
  )
  rownames(combined) <- NULL
  full_output <- combined %>% relocate(., sample, .before = A1) %>% relocate(., tool, .before = A1)
  write.table(full_output, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_all_calls_mhci.tsv", sep = ""))
  output <- rbind.data.frame(A_df, B_df, C_df) %>% 
    dplyr::mutate(method = "majority_vote") %>% 
    relocate(method, .before = weight_winner)
  output <- output %>%
    left_join(depth_df, by = c("sample", "gene"))
  output_clean <- output %>% dplyr::select(sample, gene, allele1, allele2, matching_tools, method, support, mean_depth_hla_exons_2_3_gene)
  write.table(output, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_votes_mhci_stats.tsv", sep = ""))
  write.table(output_clean, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_votes_mhci.tsv", sep = ""))
} else if (method == "weighted") {
  print("applying weighted voting")
  tools <- c("hlala", "kourami", "optitype", "polysolver")
  if (!(all(weights$tool %in% tools))){
    stop("names of tools in tool column of weights file are incorrect")
  }
  alleles <- c("A", "B", "C")
  if (!all(alleles %in% colnames(weights))){
    stop("Allele names in weights file are incorrect")
  }
  weights_selected <- weights %>% dplyr::select(tool, A, B, C)
  weights_A <- weights_selected$A
  names(weights_A) <- weights_selected$tool
  
  weights_B <- weights_selected$B
  names(weights_B) <- weights_selected$tool
  
  weights_C <- weights_selected$C
  names(weights_C) <- weights_selected$tool
  A_weighted_vote <- weighted_vote(A_list_notna, weights_A)
  A_vote <- A_weighted_vote$genotype 
  A_data <- data.frame(A_weighted_vote[2:length(A_weighted_vote)])
  matching_tools_A <- extract_tools_contributing_to_vote(A_list_notna, A_vote)
  A_df <- data.frame(
                    sample = samplename,
                    gene = "HLA-A",
                    allele1 = A_vote[1],
                    allele2 = A_vote[2]
                    )
  A_df <- cbind(A_df, matching_tools_A) %>% 
          dplyr::rename("matching_tools" = matching_tools_A) 
  A_df <- cbind(A_df, A_data) %>% 
    dplyr::relocate(support, .before = matching_tools)
  B_weighted_vote <- weighted_vote(B_list_notna, weights_B)
  B_vote <- B_weighted_vote$genotype
  B_data <- data.frame(B_weighted_vote[2:length(B_weighted_vote)])
  matching_tools_B <- extract_tools_contributing_to_vote(B_list_notna, B_vote)
  B_df <- data.frame(
    sample = samplename,
    gene = "HLA-B",
    allele1 = B_vote[1],
    allele2 = B_vote[2]
  )
  B_df <- cbind(B_df, matching_tools_B) %>% 
    dplyr::rename("matching_tools" = matching_tools_B) 
  B_df <- cbind(B_df, B_data) %>% 
    dplyr::relocate(support, .before = matching_tools)
  C_weighted_vote <- weighted_vote(C_list_notna, weights_C)
  C_vote <- C_weighted_vote$genotype
  C_data <- data.frame(C_weighted_vote[2:length(C_weighted_vote)])
  matching_tools_C <- extract_tools_contributing_to_vote(C_list_notna, C_vote)
  C_df <- data.frame(
    sample = samplename,
    gene = "HLA-C",
    allele1 = C_vote[1],
    allele2 = C_vote[2]
  )
  C_df <- cbind(C_df, matching_tools_C) %>% 
    dplyr::rename("matching_tools" = matching_tools_C) 
  C_df <- cbind(C_df, C_data) %>% 
    dplyr::relocate(support, .before = matching_tools)
  rownames(combined) <- NULL
  full_output <- combined %>% relocate(., sample, .before = A1) %>% relocate(., tool, .before = A1)
  write.table(full_output, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_all_calls_mhci.tsv", sep = ""))
  output <- rbind.data.frame(A_df, B_df, C_df) %>% 
            dplyr::mutate(method = "weighted_vote") %>% 
            relocate(method, .before = weight_winner)
  output <- output %>%
    left_join(depth_df, by = c("sample", "gene"))
  output_clean <- output %>% dplyr::select(sample, gene, allele1, allele2, matching_tools, method, support, mean_depth_hla_exons_2_3_gene)
  write.table(output, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_votes_mhci_stats.tsv", sep = ""))
  write.table(output_clean, quote = F, row.names = F, sep = "\t", file = paste(samplename, "_votes_mhci.tsv", sep = ""))
} else {
  stop("--method must be majority or weighted")
}
