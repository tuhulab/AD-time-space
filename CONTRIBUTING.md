# Contributing to AD-time-space

Thank you for your interest in contributing to this project! This document provides guidelines for contributing to the AD-time-space reproducible research repository.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Code Quality Standards](#code-quality-standards)
- [Submitting Changes](#submitting-changes)
- [Reporting Issues](#reporting-issues)

## Code of Conduct

This project adheres to standard scientific research ethics and collaborative development practices. Please be respectful and constructive in all interactions.

## Getting Started

### Prerequisites

- R (version 4.0 or higher recommended)
- RStudio (optional but recommended)
- Basic understanding of transcriptomics data analysis
- Familiarity with R packages: dplyr, ggplot2, DESeq2, tidybulk

### Initial Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/tuhulab/AD-time-space.git
   cd AD-time-space
   ```

2. **Install dependencies**
   ```R
   # Install renv for package management
   install.packages("renv")
   
   # Restore project dependencies
   renv::restore()
   ```

3. **Configure your environment**
   - Update the working directory in analysis scripts to your local path
   - Ensure data files are available in the `data/` directory

## Development Setup

### Package Management with renv

This project uses `renv` to ensure reproducibility:

- **Initialize renv** (first time only):
  ```R
  source("setup_reproducibility.R")
  ```

- **Install new packages**:
  ```R
  # Install package as usual
  install.packages("package_name")
  
  # Update renv.lock
  renv::snapshot()
  ```

- **Update existing packages**:
  ```R
  renv::update()
  ```

### Project Structure

```
AD-time-space/
├── R/                    # Reusable R functions
│   ├── helper.R         # Data manipulation helpers
│   ├── de_space.R       # Spatial DE analysis
│   ├── de_time.R        # Temporal DE analysis
│   └── plot.R           # Visualization functions
├── data/                # Data files (not tracked in git)
├── *.Rmd                # Analysis notebooks
├── DESCRIPTION          # Package metadata
└── renv.lock            # Dependency specifications
```

## Code Quality Standards

### R Code Style

We follow the [tidyverse style guide](https://style.tidyverse.org/) with these key points:

1. **Naming conventions**:
   - Use snake_case for function and variable names
   - Use meaningful, descriptive names
   - Avoid abbreviations unless widely recognized

2. **Function documentation**:
   - Include roxygen2-style documentation for all functions
   - Document parameters, return values, and examples
   - Use `@param`, `@return`, `@examples` tags

3. **Code organization**:
   - Group related functions together
   - Use section headers (# Section ----) for navigation
   - Keep functions focused and modular

4. **Error handling**:
   - Validate function inputs
   - Use informative error messages
   - Handle edge cases gracefully

### Code Quality Checks

Before submitting code:

1. **Run lintr**:
   ```R
   source("check_code_quality.R")
   ```

2. **Check for common issues**:
   - No hardcoded absolute paths
   - No sensitive information in code
   - All dependencies declared
   - Code is well-commented

3. **Test your changes**:
   - Ensure scripts run without errors
   - Verify output is as expected
   - Check that dependencies are correctly specified

### Commit Messages

Use clear, descriptive commit messages:

```
Add variance partition plotting function

- Implement plotVarPart_internal() for variance visualization
- Add comprehensive roxygen2 documentation
- Include input validation and error handling
```

Format:
- First line: Brief summary (50 chars or less)
- Blank line
- Detailed description with bullet points if needed

## Submitting Changes

### Pull Request Process

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**:
   - Write clean, documented code
   - Follow the style guide
   - Update relevant documentation

3. **Test your changes**:
   ```R
   # Run code quality checks
   source("check_code_quality.R")
   
   # Test affected scripts
   rmarkdown::render("affected_script.Rmd")
   ```

4. **Update renv.lock if needed**:
   ```R
   renv::snapshot()
   ```

5. **Commit and push**:
   ```bash
   git add .
   git commit -m "Your descriptive message"
   git push origin feature/your-feature-name
   ```

6. **Create Pull Request**:
   - Provide clear description of changes
   - Reference any related issues
   - Include screenshots if relevant

### Review Process

- Maintainers will review your PR
- Address any feedback or requested changes
- Once approved, changes will be merged

## Reporting Issues

### Bug Reports

Include:
- Clear description of the problem
- Steps to reproduce
- Expected vs. actual behavior
- R session info (`sessionInfo()`)
- Relevant error messages

### Feature Requests

Include:
- Use case description
- Proposed solution
- Potential alternatives considered
- Impact on existing functionality

### Questions

For questions about:
- **Usage**: Open a discussion or issue
- **Data access**: See README for GEO accession
- **Methods**: Refer to the published paper

## Additional Resources

- [Project README](README.md)
- [Published Paper](https://www.jidonline.org/article/S0022-202X(22)02657-4/fulltext)
- [Shiny Application](https://bit.ly/34OlBal)
- [GEO Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309)

## License

By contributing, you agree that your contributions will be licensed under the same license as the project.

## Contact

For questions or concerns, please open an issue or contact the repository maintainers.

Thank you for contributing to reproducible research in dermatology!
