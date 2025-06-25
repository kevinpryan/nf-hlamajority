#' @export
majority_vote <- function(comparison, calls) {
  stopifnot(nrow(comparison) == 4)
  stopifnot(length(calls) == 4)
  stopifnot(all(c("optitype", "polysolver", "kourami", "hlala") %in% rownames(comparison)))
  # all different, go with optitype (best tool for all MHC Class I alleles)
  if (max(colSums(comparison)) == 1) {
    output <- calls[["optitype"]]
    # all agree, doesn't matter which one we go with
  } else if (max(colSums(comparison)) == 4) {
    print("all identical and all not NA - go for optitype")
    output <- calls[["optitype"]]
  } else if (max(colSums(comparison)) == 3) {
    print("max colsums is 3")
    # get cols with that colsum
    colsum_3 <- which(colSums(comparison) == 3)
    tool_name <- colnames(colsum_3)[1]
    output <- calls[[tool_name]]
    # there is a tie, max colsum is 2
  } else {
    print("2 callers give same result")
    # check if total sum is 8 - if so, we have two equal ties - go for optitype
    if (sum(comparison) == 8) {
      output <- calls[["optitype"]]
    } else {
      # otherwise, get the name of one of the tools in the tie and go for the call
      colsum_2 <- which(colSums(comparison) == 2)
      tool_name <- colnames(comparison[,colsum_2])[1]
      print("tool name to get call for")
      print(tool_name)
      output <- calls[[tool_name]]
    }
  }
  return(output)
}

#' @export
not_na <- function(vector_list) {
  not.na.names <- c()
  for (vec_name in names(vector_list)) {
    na_indices <- which(is.na(vector_list[[vec_name]]))
    if (length(na_indices) == 0) {
      not.na.names <- c(not.na.names, vec_name)
    } #else {
      #cat("NA values in", vec_name, "\n")
    #}
  }
  return(not.na.names)
}

#' @export
majority_vote2 <- function(comparison, calls, benchmark, hla) {
  stopifnot(nrow(comparison) == 4)
  stopifnot(length(calls) == 4)
  stopifnot(all(c("optitype", "polysolver", "kourami", "hlala") %in% rownames(comparison)))
  # check if there are any NAs
  #print(calls)
  print(benchmark)
  ranks <- benchmark[[hla]]
  print(ranks)
  names(ranks) <- benchmark$tool
  print(ranks)
  #print(calls[which(any(is.na()))])
  if (!anyNA(calls, recursive = TRUE)){
	  print("no NAs")
	  # all different, go with optitype (best tool for all MHC Class I alleles)
	  if (max(colSums(comparison)) == 1) {
	    output <- calls[["optitype"]]
	    # all agree, doesn't matter which one we go with
	  } else if (max(colSums(comparison)) == 4) {
	    output <- calls[["optitype"]]
	  } else if (max(colSums(comparison)) == 3) {
	    print("max colsums is 3")
	    # get cols with that colsum
	    colsum_3 <- which(colSums(comparison) == 3)
	    tool_name <- colnames(colsum_3)[1]
	    output <- calls[[tool_name]]
	    # there is a tie, max colsum is 2
	  } else {
	    print("2 callers give same result")
	    # check if total sum is 8 - if so, we have two equal ties - go for optitype
	    if (sum(comparison) == 8) {
	      output <- calls[["optitype"]]
	    } else {
	      # otherwise, get the name of one of the tools in the tie and go for the call
	      colsum_2 <- which(colSums(comparison) == 2)
	      tool_name <- colnames(comparison[,colsum_2])[1]
	      print("tool name to get call for")
	      print(tool_name)
	      output <- calls[[tool_name]]
	    }
	 }
  } else {
    print("at least one element NA")
    output <- c(NA, NA)
    not.nas <- not_na(calls)
    #print(unlist(calls[[not.nas]])) 
    print(calls[not.nas])
  }
  return(output)
}

#' @export
majority_vote3 <- function(comparison, calls, benchmark, hla) {
  #stopifnot(nrow(comparison) == 4)
  #stopifnot(length(calls) == 4)
  #stopifnot(all(c("optitype", "polysolver", "kourami", "hlala") %in% rownames(comparison)))
  # check if there are any NAs
  #print(calls)
  #print(benchmark)
  tools <- rownames(comparison)
  benchmark <- benchmark[tools,]
  ranks <- benchmark[[hla]]
  #print(ranks)
  names(ranks) <- benchmark$tool
  print(ranks)
  
  top_tool <- max(ranks)
  top_tool <- names(top_tool)
  sum_cols <- colSums(comparison)
  #print(calls[which(any(is.na()))])
  # all the same, go with any call
  if (colSums(comparison) == nrow(comparison)**2) {
      output <- calls[1]
  # all different
  } else if (max(colSums(comparison)) == 1) {
      output <- calls[[top_tool]]
  # check for ties
  } else if (length(which(sum_cols == max(sum_cols))) > 1) {
      if (top_tool %in% asdf){
      # get cols with that colsum
      colsum_3 <- which(colSums(comparison) == 3)
      tool_name <- colnames(colsum_3)[1]
      output <- calls[[tool_name]]
      # there is a tie, max colsum is 2
    } else {
      print("2 callers give same result")
      # check if total sum is 8 - if so, we have two equal ties - go for optitype
      if (sum(comparison) == 8) {
        output <- calls[["optitype"]]
      } else {
        # otherwise, get the name of one of the tools in the tie and go for the call
        colsum_2 <- which(colSums(comparison) == 2)
        tool_name <- colnames(comparison[,colsum_2])[1]
        print("tool name to get call for")
        print(tool_name)
        output <- calls[[tool_name]]
      }
    }
  } else {
    print("at least one element NA")
    output <- c(NA, NA)
    not.nas <- not_na(calls)
    #print(unlist(calls[[not.nas]])) 
    print(calls[not.nas])
  }
  return(output)
}

#' @export
majority_vote_comparison <- function(comparison, calls, benchmark, hla) {
  # if all NA
  if (nrow(comparison) == 0){
    print("all NA")
    return(c(NA,NA))
  }
  tools <- rownames(comparison)
  benchmark <- benchmark[which(benchmark$tool %in% tools),]
  ranks <- benchmark[[hla]]
  names(ranks) <- benchmark$tool
  top_tool <- max(ranks)
  top_tool <- names(top_tool)
  # Remove diagonal elements as they represent self-comparisons
  comparison[lower.tri(comparison, diag = T)] <- 0
  # Check for ties
  tie_indices <- which(comparison == max(comparison), arr.ind = TRUE)
  tie_calls <- rownames(comparison)[tie_indices[,1]]
  # Get row and column indices of maximum values
  max_indices <- which(comparison == max(comparison), arr.ind = TRUE)
  print("max_indices")
  print(max_indices)
  # If there's only one maximum value, return the corresponding row name
  if (nrow(max_indices) == 1) {
    print("only 1 max value")
    tool_chosen <- rownames(comparison)[max_indices[1,1]]
    return(calls[[tool_chosen]])
  # If there's no tie, return the call with the highest frequency
  } else if (length(tie_calls) == 1) {
    print("no tie")
    return(calls[[tie_calls[1]]])
  } else if (max(colSums(comparison)) == 0) {
    tool_chosen = rownames(comparison)[which.min(ranks[rownames(comparison)])]
    return(calls[[tool_chosen]])
  } else {
    # If there's a tie, return the call from the tool with the highest rank
    #print("there's a tie")
    tool_chosen <- tie_calls[which.min(ranks[tie_calls])]
    return(calls[[tool_chosen]])
  }
}

