#' Select globally optimal non-overlapping TADs using dynamic programming
#'
#' @description
#' This function implements a dynamic programming algorithm to select a set of
#' non-overlapping TADs that maximize the total MoC score. The algorithm sorts TADs
#' by their end coordinates and builds an optimal solution by either including or
#' excluding each TAD based on which choice yields the highest total score.
#'
#' @param tad_all Data frame containing TAD information with columns: chr, start, end, moc_score, score_source
#'
#' @return Data frame containing the selected non-overlapping TADs that maximize total score
#'
#' @importFrom dplyr filter arrange select
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' # Prepare input data
#' tad_data <- data.frame(
#'   chr = rep("chr1", 4),
#'   start = c(10000, 20000, 50000, 70000),
#'   end = c(30000, 45000, 65000, 90000),
#'   moc_score = c(2.5, 3.2, 1.8, 4.1),
#'   score_source = c("tool1, tool2", "tool1, tool3", "tool2, tool3", "tool1, tool4")
#' )
#'
#' # Select optimal TADs
#' optimal_tads <- select_global_optimal_tads(tad_data)
#' }
#'
#' @export
select_global_optimal_tads <- function(tad_all) {
  if (NROW(tad_all) == 0) {
    return(tibble::tibble(
      chr = factor(character(0)),
      start = integer(0),
      end = integer(0),
      moc_score = double(0),
      score_source = character(0)
    ))
  }
  if (NROW(tad_all) == 1) {
    res_df <- tad_all %>%
      dplyr::select(chr, start, end, moc_score, score_source)
    return(res_df)
  }
  tb_input <- tad_all %>%
    dplyr::arrange(end)
  n <- NROW(tb_input)
  dp <- numeric(n)
  prev <- integer(n)
  prev[1] <- 0
  for (i in 2:n) {
    non_overlap_region <- which(tb_input$end[seq(1, i - 1, 1)] < tb_input$start[i])
    if (length(non_overlap_region) == 0) {
      prev[i] <- 0
    } else {
      prev[i] <- max(non_overlap_region)
    }
  }
  dp[1] <- tb_input$moc_score[1]
  for (i in 2:n) {
    if (prev[i] == 0) {
      dp[i] <- max(dp[i - 1], tb_input$moc_score[i])
    } else {
      dp[i] <- max(dp[i - 1], tb_input$moc_score[i] + dp[prev[i]])
    }
  }
  selected_intervals <- list()
  i <- n
  while (i > 0) {
    if (i == 1 || dp[i] > dp[i - 1]) {
      selected_intervals <- append(selected_intervals, list(i))
      i <- prev[i]
    } else {
      i <- i - 1
    }
  }
  res_df <- tb_input[unlist(selected_intervals),] %>%
    dplyr::select(chr, start, end, moc_score, score_source)
  res_df
}
