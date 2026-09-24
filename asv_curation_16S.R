# =============================================================================
# fill_taxa_and_merge.R
#
# Purpose
# -------
# 1. Fill the "Taxa" column of an ASV taxonomy table:
#      - If Genus is present, Taxa = Genus.
#      - If Genus is NA, walk up the ranks (Family -> Order -> Class ->
#        Phylum -> Kingdom) and use the first non-NA rank found, formatted
#        as "<Rank>_unclassified" (e.g., "Campylobacterales_unclassified").
# 2. Merge the filled taxonomy table with an ASV x sample read-count table,
#    producing a single table ordered as: asvID, count columns, taxonomy
#    columns.
#
# Data-integrity note (confirmed with user before running)
# ----------------------------------------------------------
# Taxonomy file:  8,837 ASVs
# Count file:     8,874 ASVs
# 37 ASVs appear in the count file but have no corresponding taxonomy entry.
# Per user instruction, these 37 ASVs are DROPPED from the merged output
# (inner join on asvID). This script reports the dropped IDs so the
# exclusion is documented and reproducible.
#
# Inputs (edit paths below if needed)
# ------------------------------------
# - asv_tax_rl_2023_16S.csv : columns asvID,Kingdom,Phylum,Class,Order,
#   Family,Genus,Species,Taxa
# - asv_otu_rl_2023_16S.csv : first (unnamed) column = ASV ID, remaining
#   columns = per-sample read counts
#
# Outputs
# -------
# - results/2023_16S.csv          : merged asvID + counts + taxonomy (Taxa
#   column filled). Used by interdomain_cooccurrence_analysis.Rmd.
# - results/tax_filled_2023_16S.csv : taxonomy-only table (asvID + Kingdom..
#   Species + Taxa, no count columns). Used by
#   2023_16S_ReAnalyzed_Diversity.Rmd to build the phyloseq tax_table (which
#   requires taxonomy and counts as separate objects, not merged).
#
# Run this script from the repository root (paths below are relative to it).
#
# Dependencies: base R only (no external packages required)
# =============================================================================

tax_path     <- "data/asv_tax_rl_2023_16S.csv"
otu_path     <- "data/asv_otu_rl_2023_16S.csv"
out_path     <- "results/2023_16S.csv"
tax_out_path <- "results/tax_filled_2023_16S.csv"

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------
# 1. Read input files
#    na.strings = c("NA", "") ensures both the literal string "NA" and
#    blank cells are treated as missing data (matches how the CSVs are
#    coded).
# -----------------------------------------------------------------------
tax <- read.csv(tax_path, na.strings = c("NA", ""), stringsAsFactors = FALSE)

otu <- read.csv(otu_path, na.strings = c("NA", ""), stringsAsFactors = FALSE,
                check.names = FALSE)
# name first column explicitly so it can be used as the join key.
colnames(otu)[1] <- "asvID"
colnames(tax)[1] <- "asvID"

# -----------------------------------------------------------------------
# 2. Fill the Taxa column
#    Rank order from highest to lowest (Genus is the lowest rank used
#    here, consistent with the columns present in the taxonomy file).
# -----------------------------------------------------------------------
rank_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

fill_taxa <- function(kingdom, phylum, class, order, family, genus) {
  if (!is.na(genus)) {
    return(genus)
  }
  # Walk from the rank just above Genus (Family) up to Kingdom and return
  # the first non-missing value found, appended with "_unclassified".
  ranks_above_genus <- c(Family = family, Order = order, Class = class,
                         Phylum = phylum, Kingdom = kingdom)
  for (val in ranks_above_genus) {
    if (!is.na(val)) {
      return(paste0(val, "_unclassified"))
    }
  }
  # Only reached if every rank column is missing (not expected in this
  # dataset, but handled for completeness/reproducibility).
  return("Unclassified")
}

tax$Taxa <- mapply(fill_taxa,
                   tax$Kingdom, tax$Phylum, tax$Class,
                   tax$Order, tax$Family, tax$Genus)

# -----------------------------------------------------------------------
# 3. Report and drop count-only ASVs (per user-confirmed decision)
# -----------------------------------------------------------------------
tax_ids <- tax$asvID
otu_ids <- otu$asvID
otu_only <- sort(setdiff(otu_ids, tax_ids))

cat("Taxonomy ASVs:", length(tax_ids), "\n")
cat("Count-table ASVs:", length(otu_ids), "\n")
cat("ASVs in count table with no taxonomy match (dropped):", length(otu_only), "\n")
if (length(otu_only) > 0) {
  cat("Dropped ASV IDs:\n")
  print(otu_only)
}

# -----------------------------------------------------------------------
# 4. Merge: inner join on asvID -> asvID, counts, taxonomy
# -----------------------------------------------------------------------
count_cols <- setdiff(colnames(otu), "asvID")
tax_cols   <- setdiff(colnames(tax), "asvID")

merged <- merge(otu, tax, by = "asvID")           # inner join (default)
merged <- merged[, c("asvID", count_cols, tax_cols)]

cat("Final merged table:", nrow(merged), "ASVs x", ncol(merged), "columns\n")

# -----------------------------------------------------------------------
# 5. Write outputs
# -----------------------------------------------------------------------
write.csv(merged, out_path, row.names = FALSE)
cat("Written to:", out_path, "\n")

# Taxonomy-only table (asvID + Kingdom..Species + Taxa), no count columns --
# this is the shape 2023_16S_ReAnalyzed_Diversity.Rmd expects for tax_table().
write.csv(tax, tax_out_path, row.names = FALSE)
cat("Written to:", tax_out_path, "\n")
