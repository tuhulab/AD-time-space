# Differential Gene Expression Analysis - Spatial Variation
# 
# This script performs differential gene expression analysis to identify genes
# that vary across different spatial locations (biopsy areas) in atopic dermatitis samples.
#
# Dependencies: tibble, dplyr, tidyr, readr, stringr, purrr, tidybulk, 
#               tidySummarizedExperiment, BiocParallel, Biobase, DESeq2
#
# Author: Tu Hu
# Date: 2022

# Setup environment ----
# Note: Set working directory to project root before running
# setwd("path/to/project")  # Uncomment and modify as needed

# Load required libraries
pacman::p_load(
  tibble, dplyr, tidyr, readr, stringr, purrr,
  tidybulk, tidySummarizedExperiment,
  BiocParallel, Biobase, DESeq2
)

# Configure knitr options
knitr::opts_chunk$set(message = FALSE, warning = FALSE, echo = FALSE, fig.align = "center")

# Configure parallel processing
core_n <- future::availableCores()
register(MulticoreParam(ifelse(core_n <= 8, core_n - 2, core_n - 6)))

# Load data ----
if (!file.exists("data/se.rds")) {
  stop("Data file not found: data/se.rds\n",
       "Please ensure the SummarizedExperiment object is available.")
}

se <- readr::read_rds("data/se.rds")

# Perform differential expression analysis ----
# Define contrasts for spatial comparisons
DGE_space <-
  tibble(
    C1 = c("LS", "NN"),  # Lesional skin, Non-lesional (neighbor)
    C2 = c("NL", "NN")   # Non-lesional, Non-lesional (neighbor)
  ) %>%
  mutate(
    contrast = paste0(C1, "vs", C2),
    se_exp = map2(C1, C2, ~ se[, se$skin_type %in% c(.x, .y)]),
    # Design formula adjusts for subject, visit, and biopsy area
    design_f = ifelse(
      C1 == "NN" & C2 == "NN",
      "~ subject + visit + biopsy_area",
      "~ subject + skin_type + visit + biopsy_area"
    ),
    deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)),
    deseq = map(deseq, ~ DESeq(.x, parallel = TRUE)),
    res_name = map(deseq, ~ resultsNames(.x) %>% grep(pattern = "biopsy_area", ., value = TRUE))
  ) %>%
  unnest(cols = res_name) %>%
  mutate(
    # Apply shrinkage estimation for better log2 fold change estimates
    res_shrink = map2(
      deseq, res_name, 
      ~ lfcShrink(.x, coef = .y, type = "apeglm", parallel = TRUE)
    ),
    res_shrink = map(res_shrink, ~ as_tibble(.x, rownames = "gene_name"))
  )

# Save results ----
output_dir <- "data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

saveRDS(
  DGE_space %>% select(contrast, res_name, res_shrink), 
  file.path(output_dir, "dge_space.rds")
)

cat("Spatial differential gene expression analysis completed successfully.\n")
cat("Results saved to:", file.path(output_dir, "dge_space.rds"), "\n")

