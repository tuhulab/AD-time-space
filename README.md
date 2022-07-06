# AD in time and space assessed by mini-biopsy

This [bookdown](https://github.com/rstudio/bookdown) documents the reproducible data analysis pipelines for Tu Hu's PhD project.

## Data analysis pipeline

[![DOI](https://zenodo.org/badge/378928145.svg)](https://zenodo.org/badge/latestdoi/378928145)

[Assessment of Spatial and Temporal Variation in the Skin Transcriptome of Atopic Dermatitis by Use of Minimally Invasive Punch Biopsies](02-AD-transcriptomics-time-space.Rmd)

| rmd                               | description                                                                         |
|-----------------------------------|-------------------------------------------------------------------------------------|
| 1-transcriptome-data-cleaning.Rmd | Data cleaning and curation                                                          |
| 2-geo-upload.Rmd                  | Upload data to GEO                                                                  |
| 3-table1-tableS1.Rmd              | Table 1 (Baseline characteristics) and Table S1 (sample metadata, subject metadata) |
| 4-tableS2.Rmd                     | Table S2 (differential gene expression analysis)                                    |
| 5-figure1.Rmd                     | Figure 1 (PCA, heatmap, Venn)                                                       |
| 6-figure2.Rmd                     | Figure 2 (Across-study functional enrichment analysis)                              |
| 7-figure3.Rmd                     | Figure 3 (Variance partition analysis)                                              |
| 8-figure4.Rmd                     | Figure 4 (Space variation)                                                          |
| 9-figureS1.Rmd                    | Figure S1 (Transcriptome heatmap, cosine distance)                                  |
| 10-figureS2.Rmd                   | Figure S2 (Genome regulatory elements)                                              |
| 11-figureS3.Rmd                   | Figure S3 (Time variation)                                                          |
| 12-figureS4.Rmd                   | Figure S4 (Intraindividual variation)                                               |
| 13-figureS5.Rmd                   | Figure S5 (Time fluctuation of disease severity)                                    |
| 14-figureS6.Rmd                   | Figure S6 (Correlation heatmap of IL34, IL37, UGT3A2 and inflammatory biomarkers)   |

## Shiny application for data exploration and downloading
