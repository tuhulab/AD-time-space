# Quick Reference Guide

## Common Tasks

### Initial Setup

```bash
# Clone repository
git clone https://github.com/tuhulab/AD-time-space.git
cd AD-time-space

# Start R
R
```

```R
# Install renv
install.packages("renv")

# Restore dependencies
renv::restore()

# Initialize reproducibility
source("setup_reproducibility.R")
```

### Running Analyses

#### Data Cleaning
```R
rmarkdown::render("01-transcriptome-data-cleaning.Rmd")
```

#### Differential Expression

**Spatial Analysis**:
```bash
Rscript R/de_space.R
```

**Temporal Analysis**:
```bash
Rscript R/de_time.R
```

#### Generate Figures
```R
# Figure 1
rmarkdown::render("05-figure1.Rmd")

# Figure 2
rmarkdown::render("06-figure2.Rmd")

# Continue for other figures...
```

#### All Figures at Once
```R
figure_files <- list.files(pattern = "^\\d{2}-figure.*\\.Rmd$")
lapply(figure_files, rmarkdown::render)
```

### Code Quality Checks

#### Run Lintr
```R
library(lintr)

# Lint all R files
lint_dir("R")

# Lint specific file
lint("R/helper.R")

# Run full quality check
source("check_code_quality.R")
```

#### Generate Quality Badge
```bash
Rscript generate_quality_badge.R
```

### Package Management

#### Check Package Status
```R
renv::status()
```

#### Update Packages
```R
# Update all packages
renv::update()

# Update specific package
renv::update("dplyr")
```

#### Install New Package
```R
# Install package
install.packages("new_package")

# Update renv.lock
renv::snapshot()
```

#### Restore Original Packages
```R
renv::restore()
```

### Docker Usage

#### Build Image
```bash
docker build -t ad-time-space .
```

#### Run Container
```bash
docker run -v $(pwd):/usr/src ad-time-space
```

#### Interactive Container
```bash
docker run -it -v $(pwd):/usr/src ad-time-space /bin/sh
```

### Git Workflow

#### Create Feature Branch
```bash
git checkout -b feature/my-feature
```

#### Check Status
```bash
git status
```

#### Commit Changes
```bash
git add .
git commit -m "Description of changes"
```

#### Push Changes
```bash
git push origin feature/my-feature
```

### Troubleshooting

#### View Session Info
```R
sessionInfo()
```

#### Check R Version
```bash
R --version
```

#### Clear renv Cache
```R
renv::purge()
renv::restore()
```

#### Check Data Files
```bash
ls -lh data/
md5sum data/*.gz
```

## File Locations

```
AD-time-space/
├── R/                          # R utility scripts
│   ├── helper.R               # Data manipulation functions
│   ├── de_space.R             # Spatial DE analysis
│   ├── de_time.R              # Temporal DE analysis
│   └── plot.R                 # Plotting functions
├── data/                       # Data files (download from GEO)
│   ├── counts.txt.gz          # Raw counts
│   ├── se.rds                 # SummarizedExperiment
│   └── metadata/              # Clinical metadata
├── *.Rmd                       # Analysis notebooks
├── DESCRIPTION                 # Package metadata
├── renv.lock                   # Dependency versions
├── README.md                   # Project overview
├── SETUP.md                    # Setup instructions
├── CONTRIBUTING.md             # Contribution guidelines
└── REPRODUCIBILITY.md          # Reproducibility statement
```

## Key Functions

### From R/helper.R

```R
# Load and merge data
full_data <- pull_full_data("data/processed.rds")

# Merge technical replicates
merged_counts <- counttable_merge_library_fun(
  counttable_data = counts,
  lib_to_merge_vector = c("lib01", "lib02")
)
```

### From R/plot.R

```R
# Plot variance partition results
plotVarPart_internal(
  obj = variance_results,
  main = "Variance Partition",
  ylab = "Fraction of Variance Explained (%)"
)
```

## Quick Checks

### Verify Installation
```R
library(dplyr)
library(ggplot2)
library(DESeq2)
library(SummarizedExperiment)
```

### Test Data Loading
```R
se <- readRDS("data/se.rds")
dim(se)
colData(se)
```

### Check Outputs
```bash
ls -lh data/dge_*.rds
```

## Resources

- **GitHub**: https://github.com/tuhulab/AD-time-space
- **Paper**: [JID Article](https://www.jidonline.org/article/S0022-202X(22)02657-4/fulltext)
- **Shiny App**: https://bit.ly/34OlBal
- **GEO Data**: [GSE193309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309)

## Getting Help

1. Check documentation: README.md, SETUP.md, CONTRIBUTING.md
2. Search issues: https://github.com/tuhulab/AD-time-space/issues
3. Open new issue with:
   - Your session info
   - Steps to reproduce
   - Error messages

## Keyboard Shortcuts (RStudio)

- `Ctrl/Cmd + Enter`: Run current line/selection
- `Ctrl/Cmd + Shift + Enter`: Run current chunk (Rmd)
- `Ctrl/Cmd + Shift + K`: Knit document
- `Ctrl/Cmd + Shift + C`: Comment/uncomment lines
- `Ctrl/Cmd + Shift + M`: Insert pipe operator ` %>%`

---

For detailed instructions, see full documentation files.
