# Ossining Lakes Microbial Community Analysis (2023)

## Overview

Code for data analysis and figure generation for the study of microbial (16S and 18S rRNA gene) communities in three Ossining-area lakes: Lake Mohegan, Lake Rippowam, and Teatown Lake.

**Authors:** Taylor F. Gibson, Kayley Pugh, Matthew J. Harke

**Citation:** Submitted to Harmful Algae

### [All code and analyses presented here.](https://mharke.github.io/OssiningLakes/)

## Contents

Analysis of 16S (prokaryote) and 18S (eukaryote) amplicon sequence variant (ASV) data from three lakes sampled June–September 2023. Includes:

- Weighted Gene Co-expression Network Analysis (WGCNA) to identify co-occurring microbial modules
- Module–trait (environmental parameter) correlation analysis
- Inter-domain (16S vs 18S) correlation analysis within modules
- Network visualization of inter-domain associations
- Core taxa–environment time series plots

## Repository structure

```
ossining-lakes-2023/
├── README.md
├── ossining-lakes-2023.Rproj
├── .gitignore
├── ossining_lakes_analysis.Rmd      # R Markdown source (renders to docs/)
├── data/                            # Raw input data
│   ├── Table_SX_16S_2023_ASV_Count_and_Taxonomy.csv
│   ├── Table_SX_18S_2023_ASV_Count_and_Taxonomy.csv
│   ├── Ossining_Lake_Metadata_2023.csv
│   ├── lakesites.csv
│   ├── correlations_yellow_sdo.csv
│   └── correlations_blue_sdo.csv
├── results/                         # Generated outputs
│   └── wgcna/                       # WGCNA result CSVs and filtered correlations
├── figs/                            # Saved figures
└── docs/                            # Rendered HTML for GitHub Pages
    └── index.html
```

## Data availability

Raw amplicon sequences are publicly available under BioProject [PRJNA_XXXXXXX](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA_XXXXXXX).

## Working R environment

Analysis was performed in **R version 4.5.2 (2025-10-31 ucrt)**. Key packages:

- [WGCNA](https://cran.r-project.org/package=WGCNA) — Weighted correlation network analysis
- [vegan](https://cran.r-project.org/package=vegan) — Community ecology (Hellinger transformation)
- [igraph](https://cran.r-project.org/package=igraph) / [ggraph](https://cran.r-project.org/package=ggraph) — Network construction and visualization
- [tidyverse](https://www.tidyverse.org/) — Data wrangling and plotting
- [patchwork](https://patchwork.data-imaginist.com/) — Combining ggplot panels

## Rendering the analysis

To reproduce the HTML report locally:

```r
# From the project root directory in RStudio:
rmarkdown::render("ossining_lakes_analysis.Rmd",
                  output_file = "docs/index.html")
```

Or use the **Knit** button in RStudio.

## Hosting with GitHub Pages

1. Push this repository to GitHub.
2. Go to **Settings > Pages**.
3. Under **Source**, select **Deploy from a branch**.
4. Set branch to `main` and folder to `/docs`.
5. Save — your site will be live at `https://mharke.github.io/ossining-lakes-2023/`.

## Last updated

[Date]
