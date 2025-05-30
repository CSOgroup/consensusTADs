#' Select TADs using a series of decreasing thresholds
#'
#' @description
#' This function iteratively selects TADs using a series of decreasing MoC thresholds.
#' It starts with the highest threshold and gradually processes the remaining unselected
#' TADs with lower thresholds. For each iteration, it removes previously selected TADs
#' from consideration to avoid redundancy.
#'
#' @param tb_tool_sel Data frame containing TAD information with columns: chr, start, end, meta.tool
#' @param threshold_c Numeric vector, series of decreasing MoC thresholds
#' @param include_threshold Logical, whether to include TADs equal to the threshold, default is TRUE
#' @param considering_width Logical, whether to adjust scores by TAD width, default is TRUE
#'
#' @return Data frame containing the selected optimal non-overlapping TADs
#'
#' @importFrom dplyr bind_rows mutate select
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges overlapsAny
#' @importFrom tibble as_tibble
#'
#' @examples
#' \dontrun{
#' # Prepare input data
#' tad_data <- data.frame(
#'   chr = rep("chr1", 5),
#'   start = c(10000, 20000, 50000, 70000, 90000),
#'   end = c(30000, 45000, 65000, 85000, 110000),
#'   meta.tool = c("tool1", "tool2", "tool3", "tool1", "tool2")
#' )
#'
#' # Select TADs using a series of thresholds
#' selected_tads <- select_tads_by_threshold_series(
#'   tad_data,
#'   moc_cut_c = round(seq(1, 0.2, 0.05), 2)
#' )
#' }
#'
#' @export
select_tads_by_threshold_series <- function(tb_tool_sel, threshold_c,
                                            include_threshold = TRUE, considering_width = TRUE) {
  threshold_c <- sort(threshold_c, decreasing = TRUE)
  res_list <- vector("list", length(threshold_c))
  data_input_list <- vector("list", length(threshold_c))
  data_input_list[[1]] <- tb_tool_sel
  for (i in seq(threshold_c)) {
    res_list[[i]] <- select_tads_by_threshold(data_input_list[[i]],
                                              threshold_c[i], include_threshold, considering_width) %>%
      dplyr::mutate(threshold = threshold_c[i])
    if (NROW(res_list[[i]]) == 0) {
      data_input_list[[i + 1]] <- data_input_list[[i]]
      next
    }
    gr01 <- GenomicRanges::makeGRangesFromDataFrame(data_input_list[[i]], keep.extra.columns = TRUE)
    gr02 <- GenomicRanges::makeGRangesFromDataFrame(res_list[[i]], keep.extra.columns = TRUE)
    data_input_list[[i + 1]] <- gr01[!IRanges::overlapsAny(gr01, gr02)] %>%
      tibble::as_tibble() %>%
      dplyr::select(chr = seqnames, start, end, meta.tool)
    if (NROW(data_input_list[[i + 1]]) == 0) {
      break
    }
  }
  dplyr::bind_rows(res_list)
}

#' Select optimal non-overlapping TADs using a single threshold
#'
#' @description
#' This function selects a set of optimal non-overlapping TADs by first filtering
#' TADs based on the provided MoC threshold, then applying a global optimization
#' algorithm to select TADs that maximize the total score without overlaps.
#'
#' @param tb_tool_sel Data frame containing TAD information with columns: chr, start, end, meta.tool
#' @param threshold Numeric value, threshold for filtering TADs based on MoC score
#' @param include_threshold Logical, whether to include TADs equal to the threshold, default is TRUE
#' @param considering_width Logical, whether to adjust scores by TAD width, default is TRUE
#'
#' @return Data frame containing the selected optimal non-overlapping TADs
#'
#' @importFrom dplyr mutate
#'
#' @examples
#' \dontrun{
#' # Prepare input data
#' tad_data <- data.frame(
#'   chr = rep("chr1", 4),
#'   start = c(10000, 20000, 50000, 70000),
#'   end = c(30000, 45000, 65000, 90000),
#'   meta.tool = c("tool1", "tool2", "tool3", "tool4")
#' )
#'
#' # Select TADs with threshold 0.2
#' selected_tads <- select_tads_by_threshold(tad_data, threshold = 0.2)
#' }
#'
#' @export
select_tads_by_threshold <- function(tb_tool_sel, threshold, include_threshold = TRUE, considering_width = TRUE) {
  tad_all <- moc_score_filter(tb_tool_sel, threshold, include_threshold)
  if (considering_width) {
    tad_all <- tad_all %>%
      dplyr::mutate(moc_score = moc_score * (end + 1 - start) / 10000)
  }
  tad_selected <- select_global_optimal_tads(tad_all)
  tad_selected
}
