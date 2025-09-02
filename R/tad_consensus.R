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
#' @param split_vars Character vector, variables to split data by for parallel processing,
#'   default is c("chr")
#' @param core_number Integer, number of CPU cores to use for parallel processing.
#'   If set to 1 (default), sequential processing is used. If greater than 1,
#'   parallel processing with the specified number of cores is enabled
#' @param include_isolated Logical, whether to include isolated TADs (with no overlaps) when threshold is 0, default is FALSE
#' @param .skip_plan_setup Logical, internal parameter to skip future plan setup.
#'   Used when called from hierarchy function. Default is FALSE
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
                                   step = -0.05,
                                   split_vars = c("chr"),
                                   core_number = 1,
                                   include_isolated = FALSE,
                                   .skip_plan_setup = FALSE) {
  if (length(unique(df_tools$meta.tool)) == 1) {
    return(df_tools)
  }

  moc_cut_c <- round(seq(1, threshold, step), 2)

  df_grouped <- df_tools %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(split_vars)))

  data_input <- dplyr::group_split(df_grouped, .keep = TRUE)
  keys       <- dplyr::group_keys(df_grouped)

  if (core_number > 1) {
    if (!.skip_plan_setup) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)

      options(future.globals.maxSize = 100 * 1024^3)
      future::plan(future::multisession(workers = core_number))
    }

    df_consensus <- furrr::future_map2_dfr(
      data_input, seq_len(nrow(keys)),
      ~ {
        res <- select_tads_by_threshold_series(
          .x,
          threshold_c = moc_cut_c,
          include_threshold = TRUE,
          considering_width = TRUE,
          include_isolated = include_isolated
        )
        dplyr::bind_cols(dplyr::select(keys[.y, , drop = FALSE], -chr), res)
      }
    )
  } else {
    df_consensus <- purrr::map2_dfr(
      data_input, seq_len(nrow(keys)),
      ~ {
        res <- select_tads_by_threshold_series(
          .x,
          threshold_c = moc_cut_c,
          include_threshold = TRUE,
          considering_width = TRUE,
          include_isolated = include_isolated
        )
        dplyr::bind_cols(dplyr::select(keys[.y, , drop = FALSE], -chr), res)
      }
    )
  }

  df_consensus %>%
    dplyr::select(dplyr::all_of(split_vars),
                  start, end, score_source, threshold)
}
