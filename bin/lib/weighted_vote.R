#' @export
weighted_vote <- function(calls, weights) {
  box::use(dplyr)
  if (length(calls) == 0) {
    return(list(
      genotype = c(NA, NA),
      weight_winner = NA_real_
    ))
  }
  
  tools <- names(calls)
  weights <- weights[tools]
  
  canonical <- function(x) paste(sort(x), collapse = "|")
  genotype_strings <- sapply(calls, canonical)
  df <- data.frame(
    tool = tools,
    genotype = genotype_strings,
    weight = weights,
    stringsAsFactors = FALSE
  )
  
  scores <- df |> 
    dplyr$group_by(genotype) |> 
    dplyr$summarise(
      total_weight = sum(weight),
      max_weight   = max(weight),
      n_tools      = length(weight),
      tools        = list(tool),
      .groups = "drop"
    )
  
  top <- scores |> dplyr$filter(total_weight == max(total_weight))
  
  if (nrow(top) == 1) {
    winner <- top$genotype
    weight_winner <- top$total_weight
  } else {
    top <- top |> dplyr$filter(max_weight == max(max_weight))
    if (nrow(top) == 1) {
      winner <- top$genotype
      weight_winner <- top$total_weight
    } else {
      top <- top |> dplyr$filter(n_tools == max(n_tools))
      if (nrow(top) == 1) {
        winner <- top$genotype
        weight_winner <- top$total_weight
      } else {
        best_tool <- names(which.max(weights))
        winner <- df$genotype[df$tool == best_tool][1]
        weight_winner <- df$weight[df$tool == best_tool][1]
      }
    }
  }
  n_tools_called  <- length(calls)
  n_tools_support <- top$n_tools[1]
  genotype <- strsplit(winner, "\\|")[[1]]
  genotype <- genotype[order(
    as.integer(sub(":.*", "", genotype)),      # first field
    as.integer(sub(".*:", "", genotype))        # second field
  )]
  
  return(list(
    genotype        = genotype,
    support         = weight_winner / sum(weights),
    weight_winner   = weight_winner,
    total_weight    = sum(weights),
    n_tools_support = n_tools_support,
    n_tools_called  = n_tools_called
    )
  )
  
}
