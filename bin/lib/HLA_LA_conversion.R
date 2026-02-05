# script source: https://github.com/CCGGlab/mhc_genotyping/blob/main/scripts/functions/Conversion/HLA-LA_conversion.R
#' @export
toolOutputToR.HLA_LA <- function(outputFolder, mhci_only = F, trim = F) {
  
  # Get a list of specific txt files in the given output folder
  fileList <- list.files(path = outputFolder, pattern = "*_bestguess_G.txt$")
  # Define the expected loci
  if (mhci_only == T){
  loci <- c("A", "B", "C")
  } else {
    loci <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
  }
  if (length(fileList) == 0) {
    results <- data.frame(matrix(NA, nrow = 0, ncol = length(loci) * 2))
    return(results)
  }
  # Create a data frame to store the results in
  results <- data.frame(matrix(NA, nrow = length(fileList), ncol = length(loci) * 2))
  names(results) <- rep(loci, each = 2)
  IDs <- c()
  
  # For every file in the given folder, extract and store sample ID and genotype output
  for (i in 1:length(fileList)) {
    # Load in file as R object
    data <- utils::read.delim(file = paste(outputFolder, fileList[i], sep = "/"), sep = "\t", header = TRUE)
    
    # Extract sample ID
    IDs[i] <- strsplit(fileList[i], "_")[[1]][1]
        # Loop through alleles, store the type of gene and allele result
	for (row in 1:nrow(data)) {

	  # initialise allele as NA for this row
	  allele <- NA

	  # skip loci we don't care about
	  if (!data[row, 1] %in% loci) {
	    next
	  }

	  locus <- data[row, 1]

	  # extract allele if present and valid
	  if (!is.na(data[row, 3]) &&
	      grepl("[A-Z]{1,3}[0-9]*\\*[0-9]+", data[row, 3])) {
            message("extracting allele as it is present and valid...")
	    allele <- gsub("[A-Z]{1,3}[0-9]*\\*", "", data[row, 3])
	    allele <- gsub("\\s", "", allele)

	    if (trim == TRUE) {
	      message("trimming allele...")
	      allele <- paste(strsplit(allele, ":")[[1]][1:2], collapse = ":")
	    }
	  }

	  colIndex <- which(names(results) %in% locus)
          if (length(colIndex) == 0) next

          hap <- suppressWarnings(as.integer(data[row, 2]))
          if (!is.na(hap) && hap %in% c(1,2)) {
              results[i, colIndex[hap]] <- allele
          }
	  # store result (allele or NA)
	  #results[i, colIndex[data[row, 2]]] <- allele
	}	    
  }
  
  # Add the IDs as rownames to the data frame
  row.names(results) <- IDs
  
  return(results)
}
