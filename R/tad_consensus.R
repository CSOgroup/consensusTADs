#' Generate consensus TADs from multiple tool predictions
#'
#' @description
#' This function generates consensus TADs from predictions made by multiple tools.
#' It applies an iterative threshold approach to select optimal non-overlapping TADs
#' that represent the consensus across different prediction methods.
#'
#' Parallel processing is controlled by the future framework. Configure it before calling:
#' future::plan(future::multisession(workers = 4))
#'
#' @param df_tools Data frame containing TAD information with columns: chr, start, end, meta.tool
#'   where meta.tool identifies the prediction tool source
#' @param threshold Numeric value, the minimum threshold for MoC filtering, default is 0
#' @param step Numeric vector, sequence of threshold values to use in the iterative
#'   selection process, default is -0.05
#' @param split_vars Character vector, variables to split data by for parallel processing,
#'   default is c("chr")
#' @param include_isolated Logical, whether to include isolated TADs (with no overlaps) when threshold is 0, default is FALSE
#' @param consider_level Logical, whether to consider meta.tool_level when filtering overlaps, default is FALSE
#'
#' @return Data frame containing both the original tool TADs and consensus TADs with additional
#'   columns: score_source (metadata about contributing tools), threshold (the MoC threshold
#'   at which each consensus TAD was selected)
#'
#' @importFrom dplyr bind_rows mutate full_join select
#' @importFrom purrr discard
#' @importFrom tidyr replace_na
#' @importFrom furrr future_map2_dfr
#'
#' @examples
#' \dontrun{
#' tad_data <- data.frame(
#'   chr = rep("chr1", 6),
#'   start = c(10000, 20000, 50000, 12000, 22000, 48000),
#'   end = c(30000, 45000, 65000, 32000, 43000, 67000),
#'   meta.tool = c(rep("tool1", 3), rep("tool2", 3))
#' )
#'
#' # Sequential (default)
#' consensus_results <- generate_tad_consensus(tad_data)
#'
#' # Generate consensus TADs with custom threshold values
#' consensus_results <- generate_tad_consensus(
#'   tad_data,
#'   threshold = 0.8,
#'   step = -0.05
#' )
#'
#' # Parallel controlled by environment
#' options(future.globals.maxSize = 100 * 1024^3)
#' future::plan(future::multisession(workers = 4))
#' consensus_results <- generate_tad_consensus(tad_data)
#' future::plan(future::sequential)
#' }
#'
#' @export
generate_tad_consensus <- function(df_tools,
                                   threshold = 0,
                                   step = -0.05,
                                   split_vars = c("chr"),
                                   include_isolated = FALSE,
                                   consider_level = FALSE) {
  if (length(unique(df_tools$meta.tool)) == 1) {
    if (threshold == 0 && include_isolated) {
      has_level <- "meta.tool_level" %in% names(df_tools) &&
                   any(!is.na(df_tools$meta.tool_level))

      if (consider_level && has_level) {
        return(df_tools %>%
          dplyr::mutate(
            score_source = dplyr::case_when(
              !is.na(meta.tool_level) ~ paste0(meta.tool, "(", meta.tool_level, ")_1"),
              TRUE ~ paste0(meta.tool, "_1")
            ),
            threshold = 0
          ) %>%
          dplyr::select(dplyr::all_of(split_vars), start, end, score_source, threshold)
        )
      } else {
        return(df_tools %>%
          dplyr::mutate(
            score_source = paste0(meta.tool, "_1"),
            threshold = 0
          ) %>%
          dplyr::select(dplyr::all_of(split_vars), start, end, score_source, threshold)
        )
      }
    } else {
      return(tibble::tibble(
        chr = character(0),
        start = integer(0),
        end = integer(0),
        score_source = character(0),
        threshold = numeric(0)
      ))
    }
  }

  moc_cut_c <- round(seq(1, threshold, step), 2)

  df_grouped <- df_tools %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(split_vars)))

  data_input <- dplyr::group_split(df_grouped, .keep = TRUE)
  keys       <- dplyr::group_keys(df_grouped)

  df_consensus <- furrr::future_map2_dfr(
    data_input, seq_len(nrow(keys)),
    ~ {
      res <- select_tads_by_threshold_series(
        .x,
        threshold_c = moc_cut_c,
        include_threshold = TRUE,
        considering_width = TRUE,
        include_isolated = include_isolated,
        consider_level = consider_level
      )
      dplyr::bind_cols(dplyr::select(keys[.y, , drop = FALSE], -chr), res)
    }
  )

  df_consensus %>%
    dplyr::select(dplyr::all_of(split_vars),
                  start, end, score_source, threshold)
}
