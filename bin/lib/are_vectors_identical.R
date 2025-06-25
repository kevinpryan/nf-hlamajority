#' @export
are_vectors_identical <- function(vec1, vec2) {
  return(identical(sort(vec1), sort(vec2)))
}

#' @export
are_vectors_identical_vectorised <- Vectorize(are_vectors_identical)
