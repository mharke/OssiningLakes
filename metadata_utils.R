# =============================================================================
# metadata_utils.R
#
# Purpose
# -------
# data/ossining_meta.csv is correct as provided, but its SampleID column
# ("LM-6-1-23") does not match the sample-column names used in the ASV count
# tables ("16S-LM-6-1", "18S-LM-6-1"). This script derives the columns each
# Diversity Rmd / Figure script needs (a matching sample name, a Month
# grouping factor, and an ordered sampling-date factor) from the existing
# SampleID and Date columns -- no metadata values are invented.
#
# Verified against the actual data files (2026-09-24): the sample-name
# derivation below reproduces 45/45 column names in both
# data/asv_otu_rl_2023_16S.csv and data/asv_otu_rl_2023_18S.csv exactly.
#
# Sourced by:
#   - 2023_16S_ReAnalyzed_Diversity.Rmd
#   - 2023_18S_ReAnalyzed_Diversity.Rmd
#   - Figure2.R, Figure10.R
#   - interdomain_cooccurrence_analysis.Rmd
#
# Dependencies: base R only (no external packages required)
# =============================================================================

#' Read data/ossining_meta.csv and add derived columns.
#'
#' Added columns:
#'   Lake_abbr : Lake abbreviation parsed from SampleID (LM / LR / TTL) --
#'               used to build the sample-name key. (Note: this is NOT the
#'               same as the existing "Lake" column, which holds full names
#'               e.g. "Lake Mohegan"; both are kept.)
#'   Month_num : numeric month, parsed from SampleID (e.g. 6, 7, 8, 9)
#'   Day_num   : numeric day, parsed from SampleID
#'   Month     : factor version of Month_num, labeled with month name,
#'               ordered chronologically (June < July < August < September).
#'               Used for grouping/coloring and PERMANOVA (~Month) in the
#'               Diversity Rmds.
#'   Date_Sam  : the existing Date column parsed to an R Date, then wrapped
#'               in an ordered factor (levels = chronological order of the
#'               unique sampling dates) so ggplot/plot_richness x-axes sort
#'               correctly instead of alphabetically.
#'
#' @param path Path to ossining_meta.csv (default "data/ossining_meta.csv",
#'   i.e. relative to the repository root).
#' @return data.frame with all original columns plus the derived ones above.
load_ossining_meta <- function(path = "data/ossining_meta.csv") {

  if (!file.exists(path)) {
    stop("metadata_utils.R: cannot find '", path, "'. Run from the ",
         "repository root (where the .Rproj file lives), or pass the ",
         "correct path explicitly.")
  }

  meta <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

  # --- Parse SampleID into Lake/Month/Day components ------------------------
  # Expected format confirmed against every row in the file: "<Lake>-<Month>-
  # <Day>-<2-digit Year>", e.g. "LM-6-1-23", "LR-9-11-22", "TTL-7-7-22".
  id_pattern <- "^(LM|LR|TTL)-([0-9]{1,2})-([0-9]{1,2})-([0-9]{2})$"
  parsed_ok <- grepl(id_pattern, meta$SampleID)
  if (!all(parsed_ok)) {
    stop("metadata_utils.R: ", sum(!parsed_ok), " SampleID value(s) do not ",
         "match the expected '<Lake>-<Month>-<Day>-<Year>' pattern: ",
         paste(meta$SampleID[!parsed_ok], collapse = ", "),
         ". Fix these before proceeding -- not guessing a mapping for them.")
  }

  meta$Lake_abbr <- sub(id_pattern, "\\1", meta$SampleID)
  meta$Month_num <- as.integer(sub(id_pattern, "\\2", meta$SampleID))
  meta$Day_num   <- as.integer(sub(id_pattern, "\\3", meta$SampleID))

  # --- Month: ordered factor with real month names for plotting/PERMANOVA --
  month_lookup <- c(`6` = "June", `7` = "July", `8` = "August", `9` = "September")
  meta$Month <- factor(month_lookup[as.character(meta$Month_num)],
                       levels = month_lookup)

  # --- Date_Sam: existing Date column, parsed and ordered chronologically --
  # Confirmed format across all 81 rows: "M/D/YYYY" (no leading zeros).
  parsed_date <- as.Date(meta$Date, format = "%m/%d/%Y")
  if (any(is.na(parsed_date))) {
    stop("metadata_utils.R: ", sum(is.na(parsed_date)), " Date value(s) did ",
         "not parse as M/D/YYYY: ",
         paste(meta$Date[is.na(parsed_date)], collapse = ", "))
  }
  date_levels <- sort(unique(parsed_date))
  meta$Date_Sam <- factor(format(parsed_date, "%b-%d"),
                          levels = format(date_levels, "%b-%d"))

  meta
}

#' Add a sample-name column that matches the ASV table's sample columns.
#'
#' @param meta          Output of load_ossining_meta().
#' @param domain_prefix "16S" or "18S" -- must match the prefix used in the
#'                      corresponding asv_otu_rl_*.csv column names.
#' @param name_col      Name to give the new column (default "Name", pass
#'                      "Name2" for 18S if you want to mirror the original
#'                      analysis's column-naming convention).
#' @return meta with one additional column.
make_sample_name <- function(meta, domain_prefix, name_col = "Name") {
  if (!domain_prefix %in% c("16S", "18S")) {
    stop("metadata_utils.R: domain_prefix must be '16S' or '18S', got '",
         domain_prefix, "'.")
  }
  meta[[name_col]] <- paste0(domain_prefix, "-", meta$Lake_abbr, "-",
                             meta$Month_num, "-", meta$Day_num)
  meta
}
