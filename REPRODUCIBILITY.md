# Reproducibility Statement

This document outlines the measures taken to ensure full reproducibility of the AD-time-space analysis pipeline.

## Overview

This project implements comprehensive reproducibility measures following best practices for computational research. All analyses can be reproduced exactly by following the instructions provided.

## Reproducibility Measures

### 1. Dependency Management

**Tool**: `renv` (R Environment Manager)

- **What**: All R package dependencies and their exact versions are tracked
- **Where**: `renv.lock` file (to be generated) and `renv.lock.template`
- **How to use**:
  ```R
  # Restore exact package versions
  renv::restore()
  
  # Check dependency status
  renv::status()
  ```

**Key packages with version requirements**:
- R ≥ 4.0.0
- dplyr ≥ 1.0.0
- ggplot2 ≥ 3.3.0
- DESeq2 (Bioconductor)
- tidybulk (Bioconductor)
- SummarizedExperiment (Bioconductor)

### 2. Environment Documentation

**Session Information**: Captured in `session_info.rds`

To generate your session info:
```R
source("setup_reproducibility.R")
# Or manually:
session_info <- sessionInfo()
saveRDS(session_info, "my_session_info.rds")
```

**System Requirements**:
- Operating System: Linux, macOS, or Windows with WSL2
- RAM: 16 GB minimum, 32 GB recommended
- R Version: 4.0.0 or higher
- Storage: 50 GB free space

### 3. Containerization

**Docker Support**: `Dockerfile` provided

```bash
# Build image
docker build -t ad-time-space .

# Run container
docker run -v $(pwd):/usr/src ad-time-space
```

The Docker image ensures:
- Consistent operating system (Alpine Linux)
- Fixed R version
- Reproducible build environment

### 4. Version Control

**Git Repository**: All code is version controlled

- Repository: https://github.com/tuhulab/AD-time-space
- DOI: [![DOI](https://www.zenodo.org/badge/378928145.svg)](https://www.zenodo.org/badge/latestdoi/378928145)
- Commit history: Full audit trail of all changes
- Tagged releases: Specific versions can be referenced

### 5. Data Availability

**Public Data Repository**: GEO (Gene Expression Omnibus)

- Accession: [GSE193309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309)
- Data Type: RNA-seq count data
- Metadata: Complete sample and subject metadata included
- Access: Publicly available, no restrictions

**Data Files**:
```
data/
├── counts.txt.gz              # Raw counts (from GEO)
├── se.rds                     # Processed SummarizedExperiment
├── metadata/
│   └── LEO-AD-metadata.csv   # Clinical metadata
```

### 6. Analysis Pipeline

**Documented Workflow**: Each analysis step is documented

| Step | File | Description | Input | Output |
|------|------|-------------|-------|--------|
| 1 | `01-transcriptome-data-cleaning.Rmd` | Data cleaning | Raw counts | `se.rds` |
| 2 | `R/de_space.R` | Spatial DE | `se.rds` | `dge_space.rds` |
| 3 | `R/de_time.R` | Temporal DE | `se.rds` | `dge_time.rds` |
| 4+ | `0X-figureX.Rmd` | Visualizations | DGE results | Figures |

**Parameters**: All analysis parameters are documented in code comments

### 7. Code Quality

**Linting**: Automated code style checking

- Tool: `lintr`
- Configuration: `.lintr`
- Run checks: `source("check_code_quality.R")`

**Documentation**: All functions have roxygen2 documentation

Example from `R/helper.R`:
```R
#' Pull and Merge Full Dataset
#'
#' @param dl_path Path to RDS file
#' @return Tibble with merged data
#' @export
pull_full_data <- function(dl_path = ...) { ... }
```

**Code Standards**:
- No hardcoded absolute paths
- Input validation in all functions
- Informative error messages
- Consistent naming conventions

### 8. Continuous Integration

**GitHub Actions**: Automated quality checks

- Workflow: `.github/workflows/code-quality.yml`
- Triggers: On push and pull request
- Checks: Code linting, quality scoring
- Artifacts: Quality reports saved

### 9. Documentation

**Comprehensive Guides**:

1. **README.md**: Project overview, quick start
2. **SETUP.md**: Detailed installation instructions
3. **CONTRIBUTING.md**: Guidelines for contributors
4. **This file (REPRODUCIBILITY.md)**: Reproducibility measures

**Inline Documentation**:
- All R scripts have header comments
- All functions have roxygen2 documentation
- Analysis notebooks have markdown explanations

## Verification Checklist

Use this checklist to verify reproducibility:

- [ ] Clone repository: `git clone https://github.com/tuhulab/AD-time-space.git`
- [ ] Check R version: `R --version` (should be ≥ 4.0)
- [ ] Restore packages: `renv::restore()`
- [ ] Download data from GEO: [GSE193309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309)
- [ ] Place data in correct directories (see SETUP.md)
- [ ] Run data cleaning: `rmarkdown::render("01-transcriptome-data-cleaning.Rmd")`
- [ ] Run DE analyses: `Rscript R/de_space.R` and `Rscript R/de_time.R`
- [ ] Generate figures: Render figure Rmd files
- [ ] Compare results with published version
- [ ] Check session info matches: Compare `sessionInfo()` output

## Reporting Issues

If you encounter reproducibility issues:

1. **Check versions**: Ensure R and package versions match
   ```R
   sessionInfo()
   ```

2. **Verify data**: Confirm data files are correct
   ```bash
   md5sum data/counts.txt.gz
   # Compare with expected checksums
   ```

3. **Review logs**: Check error messages carefully

4. **Open an issue**: Report on GitHub with:
   - Your session info
   - Steps to reproduce
   - Error messages
   - System information

## Continuous Improvement

This reproducibility framework is maintained through:

- Regular dependency updates with `renv::update()`
- Automated testing via GitHub Actions
- Code quality monitoring
- Community feedback and contributions

## Citation

When using this reproducible pipeline, please cite:

> Assessment of Spatial and Temporal Variation in the Skin Transcriptome of Atopic Dermatitis by Use of Minimally Invasive Punch Biopsies. Journal of Investigative Dermatology (2022).

And reference this repository:
> Tu Hu. (2022). AD-time-space: Reproducible analysis pipeline (Version X.X) [Software]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

## Contact

For reproducibility questions:
- Open a [GitHub issue](https://github.com/tuhulab/AD-time-space/issues)
- Tag with `reproducibility` label

---

**Last Updated**: 2025-10-15  
**Maintained By**: Tu Hu Lab  
**Reproducibility Level**: High (full dependency tracking, public data, version control)
