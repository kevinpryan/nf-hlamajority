library(dplyr)
benchmark <- read.csv("../assets/benchmarking_results_claeys.csv")
benchmark <- benchmark %>% dplyr::filter(tool %in% c("HLA*LA", "Kourami", "Optitype", "Polysolver") & seq_type == "WES")
benchmark <- benchmark %>% dplyr::select(c(tool, A, B, C))
stopifnot(identical(c("HLA*LA", "Kourami", "Optitype", "Polysolver"), benchmark$tool))
benchmark$tool <- c("hlala", "kourami", "optitype", "polysolver")
benchmark$A <- benchmark$A/100
benchmark$B <- benchmark$B/100
benchmark$C <- benchmark$C/100
write.csv(benchmark, file = "../assets/benchmarking_results_claeys_cleaned.csv", row.names = F, quote = F)