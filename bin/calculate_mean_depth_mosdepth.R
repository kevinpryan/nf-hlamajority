#!/usr/bin/env Rscript
library(optparse)
library(stringr)

option_list = list(
  make_option(c("-s", "--samplename"), type="character", default=NULL,
              help="name of sample", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="path to input file", metavar="character"),
   make_option(c("-o", "--output"), type="character", default=NULL,
              help="name of output file", metavar="character")
   );

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 3){
  print_help(opt_parser)
  stop("3 arguments must be supplied", call.=FALSE)
}

opt$samplename -> samplename
opt$input -> input
opt$output -> output

#input <- "../assets/test-outputs/test-outputs-1000genomes/HG00097/mosdepth/HG00097.regions.bed.gz"

regions <- read.table(gzfile(input), sep = "\t", stringsAsFactors = FALSE)
colnames(regions) <- c("chr", "start", "end", "region", "depth")

# extract gene
regions$gene <- str_extract(regions$region, "HLA-[ABC]")

# per-gene mean (exons 2+3)
gene_means <- aggregate(
  depth ~ gene,
  data = regions,
  FUN = mean
)

# global mean
global_mean <- mean(regions$depth)

out <- cbind(
  data.frame(sample = samplename),
  gene_means,
  mean_depth_hla_exons_2_3_classI = round(global_mean, 2)
)

out$depth <- round(out$depth, 2)
colnames(out)[colnames(out) == "depth"] <- "mean_depth_hla_exons_2_3_gene"
#output <- "../assets/test-outputs/test-outputs-1000genomes/HG00097/mosdepth/HG00097-mosdepth-summary.tsv"
write.table(out, file = output, row.names = FALSE, quote = FALSE, sep = "\t")
