#' Calculate and filter Measure of Concordance (MoC) between TADs
#'
#' @description
#' This function calculates the Measure of Concordance (MoC) between TADs in the input data frame
#' and filters significant overlaps based on a threshold. The MoC is calculated as:
#' intersect.width^2 / (width1 * width2), where intersect.width is the overlap length between two regions,
#' and width1 and width2 are the lengths of the two regions.
#'
#' @param tb_tool_sel Data frame containing TAD coordinates. Must include columns: chr, start, end, meta.tool
#' @param moc_cut Numeric value, threshold for MoC
#' @param include_moc_cut Logical, whether to include results equal to MoC threshold, default is TRUE
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
#' @importFrom dplyr bind_cols select mutate filter arrange group_by summarise ungroup distinct
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
#' }
#'
#' @export
moc_score_filter <- function(tb_tool_sel, moc_cut, include_moc_cut = TRUE) {
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
  tb_moc_good01 <- dplyr::bind_cols(
    set1DT.long %>%
      tibble::as_tibble() %>%
      dplyr::select(start = start, end = end, meta.tool),
    set2DT.long %>%
      tibble::as_tibble() %>%
      dplyr::select(start02 = start, end02 = end, meta.tool02 = meta.tool)
  ) %>%
    dplyr::mutate(moc_score = moc_score)

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

  tb_moc <- tb_moc_good01_input %>%
    tidyr::unite("score_source", meta.tool02, moc_score, sep = "_", remove = FALSE) %>%
    dplyr::group_by(start, end, meta.tool) %>%
    dplyr::summarise(
      moc_score = sum(moc_score),
      MoC_source = paste(sort(score_source), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::mutate(meta.tool = paste0(meta.tool, "_1")) %>%
    tidyr::unite("meta.tool", MoC_source, meta.tool, sep = ", ") %>%
    dplyr::mutate(meta.tool = sapply(meta.tool, process_row)) %>%
    dplyr::select(start, end, moc_score, meta.tool) %>%
    dplyr::distinct() %>%
    dplyr::arrange(-moc_score)
  tad_all <- tb_moc %>%
    dplyr::mutate(
      chr = unique(tb_tool_sel$chr)
    ) %>%
    dplyr::rename(score_source = meta.tool)
  return(tad_all)
}
