#BiocManager::install("rtracklayer")
library(rtracklayer)
library(stringr)
gtf <- import("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz")
gtf
gtf$transcript_id_trimmed <- str_split_fixed(gtf$transcript_id, pattern = "\\.", n = 2)[,1]
genes <- c("HLA-A", "HLA-B", "HLA-C")
# HLA-A, HLA-B and HLA-C MANE select transcripts from the ENSEMBL website, information not in GTF for HLA genes
mane_select_transcripts <- c("ENST00000376809", "ENST00000412585", "ENST00000376228")
subset_gtf <- gtf[gtf$type == "exon" & 
                  gtf$gene_name %in% genes & 
                  gtf$exon_number %in% c("2", "3") & 
                  gtf$transcript_id_trimmed %in% mane_select_transcripts]
subset_gtf$score <- 0
subset_gtf$name <- paste0(subset_gtf$gene_name, "_Exon", subset_gtf$exon_number)

export.bed(con = "hla-a-b-c-exons-2-3.bed", object = subset_gtf)
