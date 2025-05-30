#' Generate consensus TADs from multiple tool predictions
#'
#' @description
#' This function generates consensus TADs from predictions made by multiple tools.
#' It applies an iterative threshold approach to select optimal non-overlapping TADs
#' that represent the consensus across different prediction methods.
#'
#' @param df_tools Data frame containing TAD information with columns: chr, start, end, meta.tool
#'   where meta.tool identifies the prediction tool source
#' @param threshold Numeric value, the minimum threshold for MoC filtering, default is 0
#' @param step Numeric vector, sequence of threshold values to use in the iterative
#'   selection process, default is -0.05
#'
#' @return Data frame containing both the original tool TADs and consensus TADs with additional
#'   columns: score_source (metadata about contributing tools), threshold (the MoC threshold
#'   at which each consensus TAD was selected)
#'
#' @importFrom dplyr bind_rows mutate full_join select
#' @importFrom purrr discard
#' @importFrom tidyr replace_na
#'
#' @examples
#' \dontrun{
#' # Prepare input data with predictions from multiple tools
#' tad_data <- data.frame(
#'   chr = rep("chr1", 6),
#'   start = c(10000, 20000, 50000, 12000, 22000, 48000),
#'   end = c(30000, 45000, 65000, 32000, 43000, 67000),
#'   meta.tool = c(rep("tool1", 3), rep("tool2", 3))
#' )
#'
#' # Generate consensus TADs with default parameters
#' consensus_results <- generate_tad_consensus(tad_data)
#'
#' # Generate consensus TADs with custom threshold values
#' consensus_results <- generate_tad_consensus(
#'   tad_data,
#'   threshold = 0.8,
#'   step = -0.05
#' )
#' }
#'
#' @export
generate_tad_consensus <- function(df_tools,
                                   threshold = 0,
                                   step = -0.05) {
  if (length(unique(df_tools$meta.tool)) == 1) {
    return(df_tools)
  }
  data_input <- df_tools %>%
    split(.$chr) %>%
    purrr::discard(function(x) NROW(x) == 0)
  moc_cut_c <- round(seq(1, threshold, step), 2)
  df_consensus <- data_input %>%
    lapply(select_tads_by_threshold_series,
           threshold_c = moc_cut_c,
           include_threshold = TRUE, considering_width = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::select(chr, start, end, score_source, threshold)
  df_consensus
}
