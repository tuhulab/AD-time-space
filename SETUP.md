# Setup and Reproducibility Guide

This guide provides detailed instructions for setting up the AD-time-space analysis environment to ensure full reproducibility of results.

## System Requirements

### Software Requirements

- **R**: Version 4.0.0 or higher
- **RStudio**: Version 1.4 or higher (recommended)
- **Git**: For version control
- **Minimum RAM**: 16 GB (32 GB recommended for DESeq2 analyses)
- **Storage**: At least 50 GB free space

### Operating System

Tested on:
- macOS (10.15 or higher)
- Linux (Ubuntu 18.04 or higher, CentOS 7 or higher)
- Windows 10 with WSL2 (Windows Subsystem for Linux)

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/tuhulab/AD-time-space.git
cd AD-time-space
```

### 2. Install R and Dependencies

#### On macOS

```bash
# Install R using Homebrew
brew install r

# Or download from CRAN: https://cran.r-project.org/bin/macosx/
```

#### On Ubuntu/Debian

```bash
# Add CRAN repository
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Install R
sudo apt-get update
sudo apt-get install r-base r-base-dev

# Install system dependencies for R packages
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev
```

### 3. Set Up R Environment

Open R or RStudio and run:

```R
# Install renv for dependency management
install.packages("renv")

# Restore project dependencies
renv::restore()

# This will install all required packages with exact versions specified in renv.lock
```

### 4. Verify Installation

```R
# Check R version
R.version.string

# Verify key packages are installed
library(dplyr)
library(ggplot2)
library(DESeq2)
library(SummarizedExperiment)

# Run setup script
source("setup_reproducibility.R")
```

## Data Setup

### Download Data from GEO

The RNA-seq data is available at GEO under accession [GSE193309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309).

```bash
# Create data directory
mkdir -p data

# Download processed data
# (Specific download commands depend on your preferred method)
```

### Expected Data Structure

```
data/
├── counts.txt.gz              # Raw count matrix
├── se.rds                     # SummarizedExperiment object
├── metadata/
│   └── LEO-AD-metadata.csv   # Clinical metadata
├── dge_space.rds             # Spatial DE results (generated)
└── dge_time.rds              # Temporal DE results (generated)
```

## Running the Analysis

### Step-by-Step Execution

1. **Data Cleaning** (requires raw data):
   ```R
   rmarkdown::render("01-transcriptome-data-cleaning.Rmd")
   ```

2. **Generate Tables**:
   ```R
   rmarkdown::render("03-table1-tableS1.Rmd")
   rmarkdown::render("04-tableS2.Rmd")
   ```

3. **Generate Figures**:
   ```R
   rmarkdown::render("05-figure1.Rmd")
   rmarkdown::render("06-figure2.Rmd")
   # ... and so on
   ```

4. **Run Differential Expression** (requires computational resources):
   ```bash
   # For spatial analysis
   Rscript R/de_space.R
   
   # For temporal analysis
   Rscript R/de_time.R
   ```

### Using the Bookdown Site

To generate the complete bookdown site:

```R
bookdown::render_book("index.Rmd", "bookdown::gitbook")
```

The output will be in the `_book/` directory.

## Docker Setup (Alternative)

For maximum reproducibility, use Docker:

```bash
# Build Docker image
docker build -t ad-time-space .

# Run container
docker run -v $(pwd):/usr/src ad-time-space
```

## Reproducibility Checklist

- [ ] R version matches (4.0+)
- [ ] All packages installed via `renv::restore()`
- [ ] Data files downloaded and placed in correct directories
- [ ] Working directory set correctly in scripts
- [ ] Sufficient computational resources available
- [ ] Session info captured: `sessionInfo()` or `session_info.rds`

## Environment Verification

Create and check your session info:

```R
# Capture current session
session_info <- sessionInfo()
print(session_info)
saveRDS(session_info, "my_session_info.rds")

# Compare with original
original_session <- readRDS("session_info.rds")
# Check for major version differences
```

## Code Quality Verification

Run code quality checks:

```R
# Run lintr on all R files
source("check_code_quality.R")

# Check specific file
lintr::lint("R/helper.R")
```

## Troubleshooting

### Common Issues

#### Issue: Package installation fails

**Solution**: Install system dependencies first
```bash
# Ubuntu/Debian
sudo apt-get install build-essential gfortran

# macOS
xcode-select --install
```

#### Issue: Memory errors during DESeq2

**Solution**: Increase available memory or use parallel processing
```R
# Reduce number of cores
BiocParallel::register(MulticoreParam(workers = 2))

# Or increase RAM allocation
options(java.parameters = "-Xmx8g")
```

#### Issue: File paths not found

**Solution**: Check and update working directory
```R
# Check current directory
getwd()

# Set to project root
setwd("/path/to/AD-time-space")
```

#### Issue: renv::restore() fails

**Solution**: Clear renv cache and retry
```R
renv::purge()
renv::restore()
```

## Performance Optimization

### For Large Datasets

1. **Use parallel processing**:
   ```R
   library(BiocParallel)
   register(MulticoreParam(workers = 8))
   ```

2. **Increase memory**:
   ```R
   # In bash before starting R
   ulimit -s unlimited
   ```

3. **Use faster linear algebra libraries**:
   - Install OpenBLAS or MKL
   - Link R to optimized BLAS

## Getting Help

- **Issues**: Open a GitHub issue
- **Questions**: See [CONTRIBUTING.md](CONTRIBUTING.md)
- **Data access**: Contact authors or check GEO
- **Citations**: See [published paper](https://www.jidonline.org/article/S0022-202X(22)02657-4/fulltext)

## Session Info Template

Always include session info when reporting issues:

```R
sessionInfo()
# Or more detailed:
devtools::session_info()
```

## Updates and Maintenance

To update the environment:

```R
# Update all packages
renv::update()

# Create new snapshot
renv::snapshot()

# Check for updates
renv::status()
```

## Additional Resources

- [renv documentation](https://rstudio.github.io/renv/)
- [Bioconductor installation guide](https://www.bioconductor.org/install/)
- [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

---

Last updated: 2025
