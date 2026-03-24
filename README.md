# Ossining Lakes Microbial Community Analysis (2023)

## Overview

Associated code for data analysis and figure generation for the study of microbial (16S and 18S rRNA gene) communities in three Ossining-area lakes: Lake Mohegan (LM), Lake Rippowam (LR), and Teatown Lake (TTL).

**Authors:** Taylor F. Gibson, Kayley Pugh, Matthew J. Harke

**Citation:** Gibson, T.F., Pugh, K., & Harke, M.J. (*Submitted*) [Manuscript title]. *Harmful Algae*.

### [All code and analyses presented here.](https://mharke.github.io/OssiningLakes/)

## Description of contents

### 1. Study background

Three lakes in the Ossining, NY area were sampled from June to September 2023 to characterize prokaryotic (16S rRNA gene) and eukaryotic (18S rRNA gene) microbial communities using amplicon sequencing. Environmental metadata including nutrients (nitrate, phosphate), temperature, pH, dissolved oxygen, and microcystin were measured at each sampling event. This analysis investigates the co-occurrence of microbial taxa across domains and their relationships with environmental conditions, with a focus on harmful algal bloom (HAB)-associated taxa such as *Planktothrix*.

### 2. Working R environment

Using R version 4.5.2 (2025-10-31 ucrt) with RMarkdown. Key packages:

- [WGCNA](https://cran.r-project.org/package=WGCNA) — Weighted correlation network analysis
- [vegan](https://cran.r-project.org/package=vegan) — Community ecology (Hellinger transformation, PERMANOVA)
- [phyloseq](https://joey711.github.io/phyloseq/) — Microbiome data handling and visualization
- [microbiome](https://microbiome.github.io/tutorials/) — Alpha diversity metrics
- [igraph](https://cran.r-project.org/package=igraph) / [ggraph](https://cran.r-project.org/package=ggraph) — Network construction and visualization
- [tidyverse](https://www.tidyverse.org/) — Data wrangling and plotting
- [patchwork](https://patchwork.data-imaginist.com/) — Combining ggplot panels
- [metagMisc](https://github.com/vmikk/metagMisc) — Pairwise PERMANOVA

### 3. 16S and 18S Diversity Analysis (2023)

Analysis of prokaryotic and eukaryotic community diversity, composition, and structure across three lakes over the 2023 sampling season.

- Import DADA2-derived ASV count tables and taxonomy
- Alpha diversity (observed richness)
- Community composition at phylum (16S) and division (18S) levels
- Cyanobacteria, chytrid, and metazoa community composition
- PCoA ordination (Bray-Curtis dissimilarity) per lake and overall
- PERMANOVA and pairwise comparisons (Month, Lake, Lake × Month)

> **Note:** The diversity sections require DADA2-derived ASV tables (`asv_otu_rl.csv`, `asv_tax_rl.csv`) that are generated upstream and are not included in this repository. These sections are included with `eval = FALSE` to document the methods. Raw sequences are available under BioProject PRJNA_XXXXXXX.

### 4. WGCNA — Weighted Gene Co-expression Network Analysis

Combined 16S and 18S ASV data analyzed using WGCNA to identify co-occurring modules of microbial taxa.

- Data loading, Hellinger transformation, and low-abundance ASV filtering
- Sample clustering and outlier detection
- Soft-threshold power selection (signed network)
- Network construction and module detection (`blockwiseModules`)
- Module eigengene dendrogram (Figure S11)

### 5. Module-Trait Correlation

Correlate module eigengenes with environmental parameters to identify modules responsive to specific conditions.

- Pearson correlation of module eigengenes vs. environmental traits
- Labeled heatmap of module-trait relationships (Figure 8)
- Gene Significance (GS) and Module Membership (MM) scatter plots

### 6. Intramodular Connectivity and Hub Taxa

Identify hub ASVs within key modules based on intramodular connectivity.

- Adjacency matrix calculation
- Intramodular connectivity (kWithin)
- Top hub taxa identification per module

### 7. Inter-Domain (16S vs 18S) Correlation Analysis

Cross-domain correlations within each WGCNA module to identify significant prokaryote-eukaryote associations.

- Inter-domain Pearson correlation within each module
- Significance filtering (p < 0.05, then stricter p < 0.01 and |r| > 0.4)
- Module membership and taxonomy annotation of significant pairs

### 8. Network Visualization

Network graph visualization of filtered inter-domain correlations.

- igraph/ggraph network plots per module
- Node coloring by taxonomy, shape by domain (16S/18S)
- Edge coloring by correlation sign, width by correlation strength

### 9. Module-Trait-Taxa Summary

Bubble chart summarizing the top taxon per module-trait combination (Figure 9).

- Most abundant taxon per module-trait pair
- Average Gene Significance direction and strength

## Repository structure

```
ossining-lakes-2023/
├── README.md
├── ossining-lakes-2023.Rproj
├── .gitignore
├── .gitattributes
├── ossining_lakes_analysis.Rmd      # R Markdown source (renders to docs/)
├── data/                            # Input data files
│   ├── Table_SX_16S_2023_ASV_Count_and_Taxonomy.csv
│   ├── Table_SX_18S_2023_ASV_Count_and_Taxonomy.csv
│   ├── Ossining_Lake_Metadata_2023.csv
│   ├── lakesites.csv
│   ├── 16S_2022_Abundance_FIXED.xlsx    (supplementary reference)
│   ├── 16S_2023_Abundance_FIXED.xlsx    (supplementary reference)
│   ├── 18S_2023_Abundance_FIXED.xlsx    (supplementary reference)
│   ├── Ossining_Lake_Metadata_2023.xlsx (supplementary reference)
│   └── Ossining_ANOVA_Stats.xlsx        (supplementary reference)
├── results/                         # Generated outputs
│   └── wgcna/                       # WGCNA result CSVs
├── figs/                            # Saved figures
└── docs/                            # Rendered HTML for GitHub Pages
    └── index.html
```

## Data availability

Raw amplicon sequences are publicly available under BioProject [PRJNA_XXXXXXX](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA_XXXXXXX).

## Rendering the analysis

To reproduce the HTML report locally:

```r
# From the project root directory in RStudio, click Knit or run:
rmarkdown::render("ossining_lakes_analysis.Rmd",
                  output_file = "docs/index.html")
```

The `knit:` header in the Rmd file will automatically render to `docs/index.html`.

## Hosting with GitHub Pages

1. Push this repository to GitHub.
2. Go to **Settings > Pages**.
3. Under **Source**, select **Deploy from a branch**.
4. Set branch to `main` and folder to `/docs`.
5. Save — your site will be live at `https://mharke.github.io/OssiningLakes/`.

## Last updated

[Date]
