#' Generate hierarchical consensus TADs through iterative rounds
#'
#' @description
#' This function generates consensus TADs through multiple rounds of iteration.
#' In each round, it identifies consensus TADs and removes partially overlapping
#' regions from the input data for the next round. This hierarchical approach
#' helps identify TADs at different levels of consensus strength.
#'
#' Parallel processing is controlled by the future framework. Configure it before calling:
#' future::plan(future::multisession(workers = 4))
#'
#' @param df_tools Data frame containing TAD information with columns: chr, start, end, meta.tool
#' @param threshold Numeric value, the minimum threshold for MoC filtering, default is 0
#' @param step Numeric vector, sequence of threshold values to use in the iterative
#'   selection process, default is -0.05
#' @param split_vars Character vector, variables to split data by for parallel processing,
#'   default is c("chr")
#' @param max_round Integer, maximum number of rounds to perform. If NULL, continues
#'   until no more TADs remain in the input data. Default is 10
#' @param include_isolated Logical, whether to include isolated TADs (with no overlaps) when threshold is 0, default is FALSE
#' @param consider_level Logical, whether to consider meta.tool_level when filtering overlaps, default is FALSE
#'
#' @return Data frame containing all consensus TADs with round information
#'
#' @importFrom dplyr bind_rows mutate anti_join
#' @importFrom GenomicRanges GRanges findOverlaps pintersect width
#' @importFrom IRanges IRanges
#' @importFrom tibble as_tibble
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- generate_tad_consensus_hierarchy(tad_data, max_round = 5)
#'
#' # Parallel processing
#' options(future.globals.maxSize = 100 * 1024^3)
#' future::plan(future::multisession(workers = 4))
#' result <- generate_tad_consensus_hierarchy(tad_data, max_round = 5)
#' future::plan(future::sequential)
#'
#' # With tool levels
#' tad_data_with_level <- data.frame(
#'   chr = rep("chr1", 8),
#'   start = c(10000, 15000, 20000, 50000, 55000, 15000, 50000, 80000),
#'   end = c(30000, 35000, 45000, 70000, 75000, 35000, 70000, 100000),
#'   meta.tool = c("tool1", "tool1", "tool2", "tool3", "tool3", "tool2", "tool1", "tool4"),
#'   meta.tool_level = c("L1", "L2", NA, "L1", "L2", NA, "L2", NA)
#' )
#'
#' result_hierarchy <- generate_tad_consensus_hierarchy(
#'   tad_data_with_level,
#'   max_round = 3,
#'   consider_level = TRUE
#' )
#' }
#'
#' @export
generate_tad_consensus_hierarchy <- function(df_tools,
                                             threshold = 0,
                                             step = -0.05,
                                             split_vars = c("chr"),
                                             max_round = 10,
                                             include_isolated = FALSE,
                                             consider_level = FALSE) {

  consensus_list <- list()
  df_remaining <- df_tools
  round_num <- 1

  continue_condition <- if (is.null(max_round)) {
    function(round, nrow_remaining) nrow_remaining > 0
  } else {
    function(round, nrow_remaining) round <= max_round && nrow_remaining > 0
  }

  while (continue_condition(round_num, nrow(df_remaining))) {

    round_consensus <- generate_tad_consensus(
      df_remaining,
      threshold = threshold,
      step = step,
      split_vars = split_vars,
      include_isolated = include_isolated,
      consider_level = consider_level
    )

    if (nrow(round_consensus) == 0) {
      break
    }

    round_consensus$round <- round_num
    consensus_list[[round_num]] <- round_consensus

    df_remaining_diff <- dplyr::anti_join(
      df_remaining,
      round_consensus,
      by = c("chr", "start", "end")
    )

    if (nrow(df_remaining_diff) == 0) {
      break
    }

    gr_remaining <- GenomicRanges::GRanges(
      seqnames = df_remaining_diff$chr,
      ranges = IRanges::IRanges(start = df_remaining_diff$start, end = df_remaining_diff$end)
    )

    gr_consensus <- GenomicRanges::GRanges(
      seqnames = round_consensus$chr,
      ranges = IRanges::IRanges(start = round_consensus$start, end = round_consensus$end)
    )

    df_intersect <- GenomicRanges::findOverlaps(gr_remaining, gr_consensus) %>%
      tibble::as_tibble()

    if (nrow(df_intersect) > 0) {
      gr_remaining_long <- gr_remaining[df_intersect$queryHits]
      gr_consensus_long <- gr_consensus[df_intersect$subjectHits]

      intersect_width <- GenomicRanges::pintersect(gr_remaining_long, gr_consensus_long) %>%
        GenomicRanges::width() %>%
        as.double()
      remaining_width <- gr_remaining_long %>%
        GenomicRanges::width() %>%
        as.double()
      consensus_width <- gr_consensus_long %>%
        GenomicRanges::width() %>%
        as.double()

      is_partial_overlap <- (intersect_width != remaining_width) & (intersect_width != consensus_width)
      partial_overlap_query_idx <- unique(df_intersect$queryHits[is_partial_overlap])

      if (length(partial_overlap_query_idx) > 0) {
        df_remaining <- df_remaining_diff[-partial_overlap_query_idx, ]
      } else {
        df_remaining <- df_remaining_diff
      }
    } else {
      df_remaining <- df_remaining_diff
    }

    if (nrow(df_remaining) == 0) {
      break
    }

    round_num <- round_num + 1
  }

  consensus_combined <- dplyr::bind_rows(consensus_list)

  return(consensus_combined)
}
