extract_na_tools <- function(genotype_list){
  na_call <- c(NA, NA)
  result <- names(genotype_list[sapply(genotype_list, function(x) setequal(x, na_call))])
  result_ordered <- paste(result[order(result)], collapse = ",")
  if(result_ordered == c("")){
    result_ordered = "none"
  } 
  return(result_ordered)
}
