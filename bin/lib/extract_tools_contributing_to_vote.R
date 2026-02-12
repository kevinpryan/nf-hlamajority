#' @export
extract_tools_contributing_to_vote <- function(genotype_list, genotype_call_vector){
  if (all(is.na(genotype_call_vector))){
    result_ordered <- "" 
  } else {
  result <- names(genotype_list[sapply(genotype_list, function(x) setequal(x, genotype_call_vector))])
  result_ordered <- paste(result[order(result)], collapse = ",")
  }
  return(result_ordered)
}