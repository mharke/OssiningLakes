# Ossining Lakes Microbial Community Analysis

## Overview

Associated code for data analysis and figure generation for the study of microbial (16S and 18S rRNA gene) communities in three Ossining-area lakes: Lake Mohegan (LM), Lake Rippowam (LR), and Teatown Lake (TTL).

**Authors:** Taylor F. Gibson, Kayley Pugh, Matthew J. Harke

**Citation:** Gibson, T.F., Pugh, K., & Harke, M.J. (*In Press*) Spatiotemporal Dynamics of Whole Lake Microbiomes and potential Inter-Domain Associations with Toxic Cyanobacterial Bloom. *Harmful Algae*.

## Description of contents

### 1. Study background

Three lakes in the Ossining, NY area were sampled from June to September 2023 to characterize prokaryotic (16S rRNA gene) and eukaryotic (18S rRNA gene) microbial communities using amplicon sequencing. Environmental metadata including nutrients (nitrate, phosphate), temperature, pH, dissolved oxygen, and microcystin were measured at each sampling event. This analysis investigates the co-occurrence of microbial taxa across domains and their relationships with environmental conditions, with a focus on harmful algal bloom (HAB)-associated taxa such as *Planktothrix*.

### 2. Working R environment

Using R version 4.5.2 (2025-10-31 ucrt) with RMarkdown. Key packages:

- [vegan](https://cran.r-project.org/package=vegan) — Community ecology (Hellinger transformation, PERMANOVA, Mantel tests)
- [phyloseq](https://joey711.github.io/phyloseq/) — Microbiome data handling and visualization
- [microbiome](https://microbiome.github.io/tutorials/) — Alpha diversity metrics
- [circlize](https://cran.r-project.org/package=circlize) — Circos-style network plots
- [pheatmap](https://cran.r-project.org/package=pheatmap) — Heatmaps of taxon co-occurrence and environmental correlations
- [tidyverse](https://www.tidyverse.org/) — Data wrangling and plotting
- [patchwork](https://patchwork.data-imaginist.com/) — Combining ggplot panels
- [metagMisc](https://github.com/vmikk/metagMisc) — Pairwise PERMANOVA (GitHub only: `remotes::install_github("vmikk/metagMisc")`)
- [SpiecEasi](https://github.com/zdk123/SpiecEasi) — Cross-domain compositional association estimation, used in `interdomain_cooccurrence_analysis.Rmd` (GitHub only: `remotes::install_github("zdk123/SpiecEasi")`)
- [NetCoMi](https://github.com/stefpeschel/NetCoMi) — Network construction/comparison, used in `interdomain_cooccurrence_analysis.Rmd` (GitHub only: `remotes::install_github("stefpeschel/NetCoMi", repos = c("https://cloud.r-project.org/", BiocManager::repositories()))`)

### 3. Running the pipeline in order

```r
# 1. Curate raw ASV/taxonomy tables (writes results/2023_16S.csv,
#    results/2023_18S.csv, and the tax-only tables the Diversity Rmds need)
source("asv_curation_16S.R")
source("asv_curation_18S.R")

# 2. Knit the per-domain diversity analyses
rmarkdown::render("2023_16S_ReAnalyzed_Diversity.Rmd")
rmarkdown::render("2023_18S_ReAnalyzed_Diversity.Rmd")

# 3. Knit the interdomain co-occurrence network analysis
#    (requires SpiecEasi/NetCoMi above; chunks default to eval = FALSE --
#    see the note at the top of that Rmd before flipping it on)
rmarkdown::render("interdomain_cooccurrence_analysis.Rmd")
```

All scripts/Rmds above assume the working directory is the repository root
(where `ossining-lakes-2023.Rproj` lives).

### 4. 16S and 18S Diversity Analysis (2023)

Analysis of prokaryotic and eukaryotic community diversity, composition, and structure across three lakes over the 2023 sampling season.

- Import DADA2-derived ASV count tables and taxonomy
- Alpha diversity (observed richness)
- Community composition at phylum (16S) and division (18S) levels
- Cyanobacteria, chytrid, and metazoa community composition
- PCoA ordination (Bray-Curtis dissimilarity) per lake and overall
- PERMANOVA and pairwise comparisons (Month, Lake, Lake × Month)

> **Note:** The DADA2-derived ASV/taxonomy tables for 2023 (`data/asv_otu_rl_2023_16S.csv`, `data/asv_tax_rl_2023_16S.csv`, and the 18S equivalents) **are included** in this repository, so the 2023 diversity sections run as committed. The three `*_ReAnalyzed_Dada2.Rmd` files that generate those tables from raw sequences are reference-only (`eval = FALSE`): they use absolute HPC cluster paths and reference databases (SILVA, CyanoSeq, PR2) that are not part of this repository. Raw sequences are available under BioProject PRJNA_XXXXXXX.

### 5. Cross-Domain (16S × 18S) Co-occurrence Network Analysis

Everything below runs in `interdomain_cooccurrence_analysis.Rmd` (chunks default to `eval = FALSE`; see the note at the top of that file). It builds one SpiecEasi/NetCoMi co-occurrence network per lake from the combined 16S + 18S count data, then characterizes cyanobacteria/HAB-taxon associations within and across those networks.

**5a. Data preparation**
- Load the curated 2023 16S/18S count+taxonomy tables and metadata
- Remove non-target 16S sequences (chloroplast/mitochondrial/organelle reads), with a full record of what was removed and why
- Aggregate ASVs to a working "Taxa" level; identify Cyanobacteria and user-specified HAB (harmful-algal-bloom) genera, flagging any HAB genus not classified as Cyanobacteriota in this taxonomy
- Per-lake sample matching and prevalence filtering before network construction

**5b. Network construction and cyanobacteria/HAB focus**
- Cross-domain network per lake (`SpiecEasi::multi.spiec.easi` → `NetCoMi::netConstruct`/`netAnalyze`), saved as an SVG per lake with a numbered taxon lookup table (taxon, domain, HAB status, degree, centrality measures, hub status)
- Cyanobacteria-centered edge extraction (biotic co-occurrence involving cyanobacterial taxa)
- HAB-focused ego-networks (subnetworks centered on HAB genera) with matching lookup tables
- Sorted diverging ("tornado") bar plots and a faceted heatmap of each HAB taxon's top co-occurrence partners, colored by Phylum (16S non-cyanobacteria), Taxa/genus (cyanobacteria), or Taxa2 (18S) to match the diversity-analysis figures

**5c. Cross-lake comparison**
- Cross-lake network comparison via `NetCoMi::netCompare`: global network properties, Graphlet Correlation Distance (GCD), and Jaccard index of central nodes across lake pairs, exported per-pair and combined

**5d. Environmental associations**
- Mantel tests (community dissimilarity vs. environmental distance) per lake
- HAB-genus × environmental-variable Spearman correlation heatmap, BH-adjusted, per lake
- Individual taxon-vs-environmental-variable scatter plots with OLS fit (e.g., *Planktothrix* relative abundance vs. microcystin in Lake Rippowam)

**5e. Diagnostics and exports**
- Per-lake edge counts, full edge tables (association/dissimilarity/adjacency), edge-weight quantiles, and hub-taxon tables exported to `results/network_diagnostics_export/`
- HAB-genus edges at a data-driven 75th-percentile association cutoff, combined across lakes

## Repository structure

```
ossining-lakes-2023/
├── README.md
├── ossining-lakes-2023.Rproj
├── .gitignore
├── .gitattributes
├── metadata_utils.R                        # shared metadata-massaging helpers (sourced by files below)
├── asv_curation_16S.R                      # curates raw 16S ASV+taxonomy -> results/
├── asv_curation_18S.R                      # curates raw 18S ASV+taxonomy -> results/
├── 2023_16S_ReAnalyzed_Dada2.Rmd            # reference only, eval = FALSE (HPC paths)
├── 2023_18S_ReAnalyzed_Dada2.Rmd            # reference only, eval = FALSE (HPC paths)
├── 2022_16S_ReAnalyzed_Dada2.Rmd            # reference only, eval = FALSE (HPC paths)
├── 2023_16S_ReAnalyzed_Diversity.Rmd
├── 2023_18S_ReAnalyzed_Diversity.Rmd
├── 2022_16S_ReAnalyzed_Diversity.Rmd         # BLOCKED: 2022 ASV/taxonomy data not yet in data/
├── interdomain_cooccurrence_analysis.Rmd     # SpiecEasi/NetCoMi cross-domain network analysis
├── Figure1.R                                 # sampling-station map (needs a shapefile not yet in data/)
├── Figure2.R                                 # multi-lake MC vs. cyanobacteria relabund (2023 only; 2022 panel blocked)
├── Figure10.R                                # Lake Rippowam 2023 Planktothrix panel
├── data/                                     # Input data files
│   ├── asv_otu_rl_2023_16S.csv
│   ├── asv_tax_rl_2023_16S.csv
│   ├── asv_otu_rl_2023_18S.csv
│   ├── asv_tax_rl_2023_18S.csv
│   └── ossining_meta.csv
├── results/                          # Generated outputs (curated tables, network diagnostics, netCompare CSVs)
├── figs/                             # Saved figures
└── docs/                             # Rendered HTML for GitHub Pages
    └── index.html
```

**Known open items (not yet resolved as of this README update):**
- 2022 16S ASV/taxonomy data and the `Ossining_Lake_Metadata_2022.xlsx`-derived fields are not present in `data/`. `2022_16S_ReAnalyzed_Diversity.Rmd`, `2022_16S_ReAnalyzed_Dada2.Rmd`, and the 2022 panel of `Figure2.R` are blocked on this.
- The NY-state shoreline shapefile `data/NYS_Civil_Boundaries.shp/State_Shoreline.shp` used by `Figure1.R` is not present in `data/`.

## Data availability

Raw amplicon sequences are publicly available under BioProject [PRJNA1443331](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1443331).

## Last updated

September 24, 2026
