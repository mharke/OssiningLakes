# ============================================================
# toxinvscyano_planktothrix_LR2023.R
#
# PURPOSE:
#   Single-panel version of the Figure 2 analysis, restricted to:
#     - Lake:  Lake Rippowam
#     - Year:  2023
#     - Taxon: Planktothrix only (not the full 8-genus MC-producer
#              sum used in the multi-panel script)
#
# Relative abundance = (Planktothrix reads in that sample) /
#                       (TOTAL 16S reads in that sample) * 100
#   Denominator kept consistent with the multi-panel script
#   (% of whole 16S community), per your prior decision.
#
# Microcystin (MC) source: ossining_meta.csv. Samples with missing
# MC are excluded (same "_nobdl" convention as before).
#
# NOTE: This uses only the 2023 16S table, so the Leptolyngbya/
# Geitlerinema genus-availability limitation from the 2022 table
# does not apply here.
# ============================================================

library(ggplot2)
library(ggpmisc)

# ------- file paths (edit to match your project structure) -------#
path_16S_2023 <- "results/2023_16S.csv"
path_meta     <- "data/ossining_meta.csv"

tax_cols <- c("asvID", "Kingdom", "Phylum", "Class", "Order",
              "Family", "Genus", "Species", "Taxa")

# ------- load 2023 16S ASV table -------#
dta_16S <- read.csv(path_16S_2023, header = TRUE, check.names = FALSE)
sample_cols <- setdiff(colnames(dta_16S), tax_cols)

# confirm Planktothrix is present before proceeding
if (!"Planktothrix" %in% dta_16S$Genus) {
  stop("Planktothrix not found in Genus column of ", path_16S_2023)
}

target_rows <- dta_16S$Taxa == "Planktothrix"

total_reads     <- colSums(dta_16S[, sample_cols, drop = FALSE])
planktothrix_reads <- colSums(dta_16S[target_rows, sample_cols, drop = FALSE])

relabund <- (planktothrix_reads / total_reads) * 100

relabund_df <- data.frame(
  RawColumn         = sample_cols,
  planktothrixabund = as.numeric(relabund),
  stringsAsFactors  = FALSE
)

# 2023 columns (e.g., "16S-LR-6-1") -> strip "16S-" prefix,
# append "-23" -> SampleID "LR-6-1-23"
relabund_df$SampleID <- paste0(sub("^16S-", "", relabund_df$RawColumn), "-23")

# ------- load metadata and merge -------#
meta <- read.csv(path_meta, header = TRUE, stringsAsFactors = FALSE)

merged <- merge(relabund_df[, c("SampleID", "planktothrixabund")],
                meta, by = "SampleID")

# ------- subset: Lake Rippowam, 2023, non-missing MC -------#
dta_lr23 <- subset(merged, Lake == "Lake Rippowam" & Year == 2023)

n_before <- nrow(dta_lr23)
dta_lr23 <- subset(dta_lr23, !is.na(MC))
message("Excluded ", n_before - nrow(dta_lr23),
        " Lake Rippowam 2023 samples with missing microcystin (MC) values.")
message("Final n for plot = ", nrow(dta_lr23))

# ------- plot -------#
lr23_plankto <- ggplot(dta_lr23, aes(x = MC, y = planktothrixabund)) +
  stat_poly_line() +
  stat_poly_eq(use_label(c("n", "f", "p")), size = 5, label.y = 0.7) +
  geom_point() +
  theme_bw() +
  labs(
    title = "Lake Rippowam 2023",
    x = "Microcystin [ug/l]",
    y = "Planktothrix Relative Abundance"
  )

lr23_plankto

# ------- save figure -------#
dir.create("figs", showWarnings = FALSE, recursive = TRUE)
ggsave("figs/Figure10.tif",
       plot = lr23_plankto, width = 6, height = 5, dpi = 600, compression = "lzw")
