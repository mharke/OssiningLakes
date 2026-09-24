# ------- load required packages -------#
# (Assumes ggplot2, ggpmisc, and patchwork are already installed,
#  as in the original script. Install with install.packages() if needed.)
library(ggplot2)
library(ggpmisc)
library(patchwork)

# ------- file paths (relative to repository root) -------#
path_16S_2022 <- "results/2022_16S.csv"
path_16S_2023 <- "results/2023_16S.csv"
path_meta     <- "data/ossining_meta.csv"

# NOT YET RESOLVED: results/2022_16S.csv does not exist in this repository
# and there is no data/asv_otu_rl_2022_16S.csv / asv_tax_rl_2022_16S.csv to
# curate it from either. This script cannot produce the 2022 panels of
# Figure 2 until that raw 2022 ASV+taxonomy data is added to data/ and run
# through an asv_curation_16S-style script (or an already-curated
# results/2022_16S.csv is supplied directly). Stopping here rather than
# silently proceeding without the 2022 data.
if (!file.exists(path_16S_2022)) {
  stop("Figure2.R: missing '", path_16S_2022, "'. This script needs 2022 ",
       "16S ASV+taxonomy data (curated the same way asv_curation_16S.R ",
       "curates the 2023 data) before it can run. See the comment above ",
       "this check for details.")
}

# ------- target genera list -------#
target_genera <- c("Aliinostoc", "Dolichospermum", "Geitlerinema",
                   "Leptolyngbya", "Microcystis", "Planktothrix",
                   "Pseudanabaena", "Radiocystis")

# non-sample columns present in the raw 16S tables
tax_cols <- c("asvID", "Kingdom", "Phylum", "Class", "Order",
              "Family", "Genus", "Species", "Taxa")

# ============================================================
# FUNCTION: calc_relabund
#   Reads a raw 16S ASV table and returns, for every sample column,
#   the relative abundance (%) of the target genera actually present
#   in that table's taxonomy, out of the TOTAL 16S reads in that
#   sample (see limitation #2 above).
# ============================================================
calc_relabund <- function(file) {
  
  dta <- read.csv(file, header = TRUE, check.names = FALSE)
  
  sample_cols <- setdiff(colnames(dta), tax_cols)
  
  # only sum genera that actually exist in this file's taxonomy
  genera_present <- intersect(target_genera, unique(dta$Genus))
  genera_missing <- setdiff(target_genera, genera_present)
  
  if (length(genera_missing) > 0) {
    message("NOTE (", file, "): target genera absent from taxonomy ",
            "and excluded from this year's sum: ",
            paste(genera_missing, collapse = ", "))
  }
  
  target_rows <- dta$Genus %in% genera_present
  
  total_reads  <- colSums(dta[, sample_cols, drop = FALSE])
  target_reads <- colSums(dta[target_rows, sample_cols, drop = FALSE])
  
  relabund <- (target_reads / total_reads) * 100
  
  data.frame(
    RawColumn     = sample_cols,
    cyanorelabund = as.numeric(relabund),
    stringsAsFactors = FALSE
  )
}

# ------- calculate relative abundance for each year -------#
relabund_2022 <- calc_relabund(path_16S_2022)
relabund_2023 <- calc_relabund(path_16S_2023)

# ------- convert raw ASV-table column names to metadata SampleIDs -------#
# 2022 columns (e.g., "LM-6-24") -> SampleID "LM-6-24-22"
relabund_2022$SampleID <- paste0(relabund_2022$RawColumn, "-22")

# 2023 columns (e.g., "16S-LM-6-1") -> strip "16S-" prefix,
# then append "-23" -> SampleID "LM-6-1-23"
relabund_2023$SampleID <- paste0(sub("^16S-", "", relabund_2023$RawColumn), "-23")

# ------- load metadata (microcystin + lake/year) -------#
meta <- read.csv(path_meta, header = TRUE, stringsAsFactors = FALSE)

# ------- merge relative abundance with metadata (inner join) -------#
merged_2022 <- merge(relabund_2022[, c("SampleID", "cyanorelabund")],
                     meta, by = "SampleID")
merged_2023 <- merge(relabund_2023[, c("SampleID", "cyanorelabund")],
                     meta, by = "SampleID")

# report any ASV-table samples that did not find a metadata match
unmatched_2022 <- setdiff(relabund_2022$SampleID, meta$SampleID)
unmatched_2023 <- setdiff(relabund_2023$SampleID, meta$SampleID)
if (length(unmatched_2022) > 0) {
  message("NOTE: 2022 samples with no metadata match (dropped): ",
          paste(unmatched_2022, collapse = ", "))
}
if (length(unmatched_2023) > 0) {
  message("NOTE: 2023 samples with no metadata match (dropped): ",
          paste(unmatched_2023, collapse = ", "))
}

# ------- combine years -------#
dta <- rbind(merged_2022, merged_2023)

# ------- exclude samples with missing microcystin (matches original "_nobdl" filtering) -------#
n_before <- nrow(dta)
dta <- subset(dta, !is.na(MC))
message("Excluded ", n_before - nrow(dta),
        " samples with missing microcystin (MC) values.")

#------- subset data by lake and year (mirrors original script) -------#
dta_lr   <- subset(dta, Lake == "Lake Rippowam")
dta_lr22 <- subset(dta_lr, Year == 2022)
dta_lr23 <- subset(dta_lr, Year == 2023)
dta_lm   <- subset(dta, Lake == "Lake Mohegan")
dta_lm22 <- subset(dta_lm, Year == 2022)
dta_lm23 <- subset(dta_lm, Year == 2023)
dta_ttl  <- subset(dta, Lake == "Teatown Lake")
dta_ttl22 <- subset(dta_ttl, Year == 2022)
dta_ttl23 <- subset(dta_ttl, Year == 2023)

#------- plot (structure/styling matches original script) -------#
lm22 <- ggplot(dta_lm22, aes(x = MC, y = cyanorelabund)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("n", "f", "p")), size = 2.5, label.y = 0.1) +
  geom_point() +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  labs(title = "Lake Mohegan 2022", x = "Microcystin [ug/l]", y = "Relative Abundance")

lr22 <- ggplot(dta_lr22, aes(x = MC, y = cyanorelabund)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("n", "f", "p")), size = 2.5, label.y = 0.1) +
  geom_point() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Lake Rippowam 2022", x = "Microcystin [ug/l]", y = "Relative Abundance")

ttl22 <- ggplot(dta_ttl22, aes(x = MC, y = cyanorelabund)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("n", "f", "p")), size = 2.5, label.y = 0.1) +
  geom_point() +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  labs(title = "Teatown Lake 2022", x = "Microcystin [ug/l]", y = "Relative Abundance")

lm23 <- ggplot(dta_lm23, aes(x = MC, y = cyanorelabund)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("n", "f", "p")), size = 2.5, label.y = 0.1) +
  geom_point() +
  theme_bw() +
  labs(title = "Lake Mohegan 2023", x = "Microcystin [ug/l]", y = "Relative Abundance")

lr23 <- ggplot(dta_lr23, aes(x = MC, y = cyanorelabund)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("n", "f", "p")), size = 2.5, label.y = 0.1) +
  geom_point() +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "Lake Rippowam 2023", x = "Microcystin [ug/l]", y = "Relative Abundance")

ttl23 <- ggplot(dta_ttl23, aes(x = MC, y = cyanorelabund)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("n", "f", "p")), size = 2.5, label.y = 0.1) +
  geom_point() +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  labs(title = "Teatown Lake 2023", x = "Microcystin [ug/l]", y = "Relative Abundance")

#------- merge plots (same layout as original: 2022 row on top, 2023 row below) -------#
lm22 + lr22 + ttl22 + lm23 + lr23 + ttl23 +
  plot_layout(ncol = 3)

# ------- save figure -------#
dir.create("figs", showWarnings = FALSE, recursive = TRUE)
ggsave("figs/Figure2.tif", width = 8, height = 5, dpi = 600, compression = "lzw")
