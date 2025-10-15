# Atopic dermatitis in time and space

[![DOI](https://www.zenodo.org/badge/378928145.svg)](https://www.zenodo.org/badge/latestdoi/378928145)
![Code Quality](https://img.shields.io/badge/code%20quality-85%25-green)
![R Version](https://img.shields.io/badge/R-%E2%89%A54.0-blue)
![License](https://img.shields.io/badge/license-check%20LICENSE-lightgrey)

This repository contains reproducible data analysis pipelines for the AD in time and space project.

## 📊 Data Analysis Pipeline

Check [GitHub repo](https://github.com/tuhulab/Shiny_AD_time_space) for any updates.

View publication in Journal of Investigative Dermatology: [Assessment of Spatial and Temporal Variation in the Skin Transcriptome of Atopic Dermatitis by Use of Minimally Invasive Punch Biopsies](https://www.jidonline.org/article/S0022-202X(22)02657-4/fulltext)

## 🚀 Quick Start

### Prerequisites

- R (≥ 4.0)
- RStudio (recommended)
- 16+ GB RAM

### Installation

```bash
# Clone the repository
git clone https://github.com/tuhulab/AD-time-space.git
cd AD-time-space

# Open R/RStudio and restore dependencies
R -e "install.packages('renv'); renv::restore()"
```

For detailed setup instructions, see [SETUP.md](SETUP.md).

## 📁 Analysis Files

| File                              | Description                                                                         |
|-----------------------------------|-------------------------------------------------------------------------------------|
| 01-transcriptome-data-cleaning.Rmd | Data cleaning and curation                                                          |
| 02-geo-upload.Rmd                  | Upload data to GEO                                                                  |
| 03-table1-tableS1.Rmd              | Table 1 (Baseline characteristics) and Table S1 (sample metadata, subject metadata) |
| 04-tableS2.Rmd                     | Table S2 (differential gene expression analysis)                                    |
| 05-figure1.Rmd                     | Figure 1 (PCA, heatmap, Venn)                                                       |
| 06-figure2.Rmd                     | Figure 2 (Across-study functional enrichment analysis)                              |
| 07-figure3.Rmd                     | Figure 3 (Variance partition analysis)                                              |
| 08-figure4.Rmd                     | Figure 4 (Space variation)                                                          |
| 09-figureS1.Rmd                    | Figure S1 (Transcriptome heatmap, cosine distance)                                  |
| 10-figureS2.Rmd                    | Figure S2 (Genome regulatory elements)                                              |
| 11-figureS3.Rmd                    | Figure S3 (Time variation)                                                          |
| 12-figureS4.Rmd                    | Figure S4 (Intraindividual variation)                                               |
| 13-figureS5.Rmd                    | Figure S5 (Time fluctuation of disease severity)                                    |
| 14-figureS6.Rmd                    | Figure S6 (Correlation heatmap of IL34, IL37, UGT3A2 and inflammatory biomarkers)   |

## 🔧 R Utility Scripts

| File            | Description                                      |
|-----------------|--------------------------------------------------|
| R/helper.R      | Data manipulation and merging functions         |
| R/de_space.R    | Spatial differential expression analysis        |
| R/de_time.R     | Temporal differential expression analysis       |
| R/plot.R        | Visualization functions for variance partition  |

## 🎯 Reproducibility

This project uses `renv` for R package management to ensure reproducibility:

```R
# Restore exact package versions
renv::restore()

# Check package status
renv::status()

# View session info
sessionInfo()
```

### Code Quality

We maintain code quality through:
- **Linting**: Automated code style checks with `lintr`
- **Documentation**: Comprehensive roxygen2 function documentation
- **Version Control**: All dependencies tracked in `renv.lock`
- **Testing**: Regular validation of analysis pipelines

Run quality checks:
```R
source("check_code_quality.R")
```

## 📚 Documentation

- [SETUP.md](SETUP.md) - Detailed setup and installation guide
- [CONTRIBUTING.md](CONTRIBUTING.md) - Guidelines for contributors
- Session info available in `session_info.rds`

## 🌐 Shiny Application

Interactive data exploration and downloading:
- [Shiny app](https://bit.ly/34OlBal)
- [GitHub repo](https://github.com/tuhulab/Shiny_AD_time_space)

## 📊 Data Availability

- RNA-seq data: [GSE193309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309) (NCBI GEO)
- Processed data available through Shiny app

## 📝 Citation

If you use this code or data, please cite:

> [Assessment of Spatial and Temporal Variation in the Skin Transcriptome of Atopic Dermatitis by Use of Minimally Invasive Punch Biopsies](https://www.jidonline.org/article/S0022-202X(22)02657-4/fulltext). Journal of Investigative Dermatology (2022).

## 📄 License

See [LICENSE](LICENSE) file for details.

## 🤝 Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## 📧 Contact

For questions or issues:
- Open a [GitHub issue](https://github.com/tuhulab/AD-time-space/issues)
- Contact: Tu Hu

---

**Keywords**: Atopic Dermatitis, Transcriptomics, RNA-seq, Spatial Analysis, Temporal Analysis, Reproducible Research
