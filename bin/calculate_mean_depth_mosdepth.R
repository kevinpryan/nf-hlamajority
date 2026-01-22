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

regions <- read.table(gzfile(input), sep = "\t")
mean_depth <- round(mean(regions$V5),2)
df <- data.frame(sample = samplename, mean_depth_hla_exons_2_3 = mean_depth)
write.csv(df, file = output, row.names = F, quote = F)
