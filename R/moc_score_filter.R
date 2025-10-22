#' Calculate and filter Measure of Concordance (MoC) between TADs
#'
#' @description
#' This function calculates the Measure of Concordance (MoC) between TADs in the input data frame
#' and filters significant overlaps based on a threshold. The MoC is calculated as:
#' intersect.width^2 / (width1 * width2), where intersect.width is the overlap length between two regions,
#' and width1 and width2 are the lengths of the two regions.
#'
#' @param tb_tool_sel Data frame containing TAD coordinates. Must include columns: chr, start, end, meta.tool.
#'   Optionally can include meta.tool_level for finer tool classification
#' @param moc_cut Numeric value, threshold for MoC
#' @param include_moc_cut Logical, whether to include results equal to MoC threshold, default is TRUE
#' @param include_isolated Logical, whether to include isolated TADs (with no overlaps) when moc_cut is 0.
#'   These TADs will have moc_score = 0. Default is FALSE
#' @param consider_level Logical, whether to consider meta.tool_level when filtering overlaps.
#'   If TRUE and meta.tool_level exists, different levels of the same tool are treated as different tools.
#'   Default is FALSE
#'
#' @return Data frame containing merged TAD information with calculated MoC scores and the following columns:
#'   \item{chr}{Character, the chromosome name where the TAD is located}
#'   \item{start}{Integer, the start coordinate of the TAD}
#'   \item{end}{Integer, the end coordinate of the TAD}
#'   \item{moc_score}{Numeric, the Measure of Concordance (MoC) score calculated for the TAD,
#'                    representing the level of agreement between different TADs}
#'   \item{score_source}{Character, a string containing information about the tools that
#'                      contributed to this TAD and their individual MoC scores}
#'
#' @importFrom GenomicRanges GRanges findOverlaps pintersect width
#' @importFrom tibble as_tibble
#' @importFrom dplyr bind_cols select mutate filter arrange group_by summarise ungroup distinct bind_rows
#' @importFrom tidyr unite
#' @importFrom stringr str_split_1
#'
#' @examples
#' \dontrun{
#' # Prepare input data
#' tad_data <- data.frame(
#'   chr = rep("chr1", 3),
#'   start = c(10000, 20000, 50000),
#'   end = c(30000, 45000, 65000),
#'   meta.tool = c("tool1", "tool2", "tool3")
#' )
#'
#' # Calculate MoC
#' results <- moc_score_filter(tad_data, moc_cut = 0.1)
#'
#' # Include isolated TADs when moc_cut is 0
#' results_with_isolated <- moc_score_filter(tad_data, moc_cut = 0, include_isolated = TRUE)
#'
#' # With tool levels
#' tad_data_with_level <- data.frame(
#'   chr = rep("chr1", 8),
#'   start = c(10000, 15000, 20000, 50000, 55000, 15000, 50000, 50000),
#'   end = c(30000, 35000, 45000, 70000, 75000, 35000, 70000, 70000),
#'   meta.tool = c("tool1", "tool1", "tool2", "tool3", "tool3", "tool2", "tool1", "tool4"),
#'   meta.tool_level = c("L1", "L2", NA, "L1", "L2", NA, "L2", NA)
#' )
#'
#' # Without considering levels - tool1(L1) and tool1(L2) are treated as same tool
#' results_no_level <- moc_score_filter(tad_data_with_level, moc_cut = 0.1, consider_level = FALSE)
#' # Output shows overlaps between tool1, tool2, tool3
#'
#' # With considering levels - tool1(L1) and tool1(L2) are treated as different tools
#' results_with_level <- moc_score_filter(tad_data_with_level, moc_cut = 0.1, consider_level = TRUE)
#' # Output shows overlaps between tool1(L1), tool1(L2), tool2, tool3(L1), tool3(L2)
#' # score_source will show format like: tool1(L1)_0.5; tool2_0.3
#' }
#'
#' @export
moc_score_filter <- function(tb_tool_sel, moc_cut, include_moc_cut = TRUE, include_isolated = FALSE, consider_level = FALSE) {
  # Check if meta.tool_level exists and has non-NA values
  has_level <- "meta.tool_level" %in% names(tb_tool_sel) &&
               any(!is.na(tb_tool_sel$meta.tool_level))

  # If no level info or all NA, set consider_level to FALSE
  if (!has_level) {
    consider_level <- FALSE
  }

  gr_whole <- GenomicRanges::GRanges(
    tb_tool_sel
  )
  set1DT <- gr_whole
  set2DT <- gr_whole
  df_intersect <- GenomicRanges::findOverlaps(set1DT, set2DT) %>%
    tibble::as_tibble() %>%
    dplyr::filter(queryHits != subjectHits)
  set1DT.long <- set1DT[df_intersect$queryHits]
  set2DT.long <- set2DT[df_intersect$subjectHits]
  intersect.width <- GenomicRanges::pintersect(set1DT.long, set2DT.long) %>%
    GenomicRanges::width() %>%
    as.double()
  set1DT.long.width <- set1DT.long %>%
    GenomicRanges::width() %>%
    as.double()
  set2DT.long.width <- set2DT.long %>%
    GenomicRanges::width() %>%
    as.double()
  moc_score <- intersect.width^2 / (set1DT.long.width * set2DT.long.width)

  consider_only_max_contribution <- FALSE
  if (consider_level && has_level) {
    tb_moc_good01 <- dplyr::bind_cols(
      set1DT.long %>%
        tibble::as_tibble() %>%
        dplyr::select(chr = seqnames, start = start, end = end, meta.tool, meta.tool_level),
      set2DT.long %>%
        tibble::as_tibble() %>%
        dplyr::select(start02 = start, end02 = end, meta.tool02 = meta.tool, meta.tool_level02 = meta.tool_level)
    ) %>%
      dplyr::mutate(moc_score = moc_score) %>%
      dplyr::filter(meta.tool != meta.tool02)

    if (consider_only_max_contribution) {
      tb_moc_good01 <- tb_moc_good01 %>%
        dplyr::group_by(chr, start, end, meta.tool, meta.tool_level, meta.tool02) %>%
        dplyr::slice_max(moc_score, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()
    }
  } else {
    tb_moc_good01 <- dplyr::bind_cols(
      set1DT.long %>%
        tibble::as_tibble() %>%
        dplyr::select(chr = seqnames, start = start, end = end, meta.tool),
      set2DT.long %>%
        tibble::as_tibble() %>%
        dplyr::select(start02 = start, end02 = end, meta.tool02 = meta.tool)
    ) %>%
      dplyr::mutate(moc_score = moc_score) %>%
      dplyr::filter(meta.tool != meta.tool02)

    if (consider_only_max_contribution) {
        tb_moc_good01 <- tb_moc_good01 %>%
          dplyr::group_by(chr, start, end, meta.tool, meta.tool02) %>%
          dplyr::slice_max(moc_score, n = 1, with_ties = FALSE) %>%
          dplyr::ungroup()
    }
  }

  process_row <- function(input_string) {
    paste(sort(stringr::str_split_1(input_string, ", ")), collapse = "; ")
  }

  if (include_moc_cut) {
    tb_moc_good01_input <- tb_moc_good01 %>%
      dplyr::filter(moc_score >= moc_cut)
  } else {
    tb_moc_good01_input <- tb_moc_good01 %>%
      dplyr::filter(moc_score > moc_cut)
  }

  if (consider_level && has_level) {
    tb_moc <- tb_moc_good01_input %>%
      dplyr::mutate(
        score_source = dplyr::case_when(
          !is.na(meta.tool_level02) ~ paste0(meta.tool02, "(", meta.tool_level02, ")_", moc_score),
          TRUE ~ paste0(meta.tool02, "_", moc_score)
        )
      ) %>%
      dplyr::group_by(chr, start, end, meta.tool, meta.tool_level) %>%
      dplyr::summarise(
        moc_score = sum(moc_score),
        MoC_source = paste(sort(score_source), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        meta.tool = dplyr::case_when(
          !is.na(meta.tool_level) ~ paste0(meta.tool, "(", meta.tool_level, ")_1"),
          TRUE ~ paste0(meta.tool, "_1")
        )
      ) %>%
      tidyr::unite("meta.tool", MoC_source, meta.tool, sep = ", ") %>%
      dplyr::mutate(meta.tool = sapply(meta.tool, process_row)) %>%
      dplyr::select(chr, start, end, moc_score, score_source = meta.tool) %>%
      dplyr::distinct() %>%
      dplyr::arrange(-moc_score)
  } else {
    tb_moc <- tb_moc_good01_input %>%
      tidyr::unite("score_source", meta.tool02, moc_score, sep = "_", remove = FALSE) %>%
      dplyr::group_by(chr, start, end, meta.tool) %>%
      dplyr::summarise(
        moc_score = sum(moc_score),
        MoC_source = paste(sort(score_source), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(meta.tool = paste0(meta.tool, "_1")) %>%
      tidyr::unite("meta.tool", MoC_source, meta.tool, sep = ", ") %>%
      dplyr::mutate(meta.tool = sapply(meta.tool, process_row)) %>%
      dplyr::select(chr, start, end, moc_score, score_source = meta.tool) %>%
      dplyr::distinct() %>%
      dplyr::arrange(-moc_score)
  }

  if (nrow(tb_moc) == 0) {
    tb_moc <- tb_moc %>%
      dplyr::mutate(score_source = character(0))
  }

  if (moc_cut == 0 && include_isolated) {
    if (consider_level && has_level) {
      tb_isolated <- tb_tool_sel %>%
        dplyr::anti_join(
          tb_moc %>% dplyr::select(start, end),
          by = c("start", "end")
        ) %>%
        dplyr::mutate(
          moc_score = 0,
          score_source = dplyr::case_when(
            !is.na(meta.tool_level) ~ paste0(meta.tool, "(", meta.tool_level, ")_1"),
            TRUE ~ paste0(meta.tool, "_1")
          )
        ) %>%
        dplyr::select(chr, start, end, moc_score, score_source)
    } else {
      tb_isolated <- tb_tool_sel %>%
        dplyr::anti_join(
          tb_moc %>% dplyr::select(start, end),
          by = c("start", "end")
        ) %>%
        dplyr::mutate(
          moc_score = 0,
          score_source = paste0(meta.tool, "_1")
        ) %>%
        dplyr::select(chr, start, end, moc_score, score_source)
    }

    tb_moc <- dplyr::bind_rows(tb_moc, tb_isolated)
  }

  return(tb_moc)
}
