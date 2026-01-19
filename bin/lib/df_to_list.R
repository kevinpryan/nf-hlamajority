#' @export
df_to_list <- function(df, cols = c()) {
  num_cols <- length(cols)
  output <- list()
  for (i in 1:nrow(df)) {
    listdf <- df[i, which(colnames(df) %in% cols)]
    vec <- c(listdf[1, 1], listdf[1, 2])
    if (!all(is.na(vec))){
    vec_sorted <- sort(vec)
    names(vec_sorted) <- cols
    output[[rownames(df)[i]]] <- vec_sorted
    } else {
        output[[rownames(df)[i]]] <- c(NA,NA)
    }
  }
  return(output)
}
