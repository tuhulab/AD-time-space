# Differential Gene Expression Analysis - Temporal Variation
# 
# This script performs differential gene expression analysis to identify genes
# that vary over time in atopic dermatitis samples, including condition-specific
# and non-specific temporal changes.
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

# Configure parallel processing
core_n <- future::availableCores()
register(MulticoreParam(ifelse(core_n <= 8, core_n - 2, 20)))

# Load data ----
if (!file.exists("data/se.rds")) {
  stop("Data file not found: data/se.rds\n",
       "Please ensure the SummarizedExperiment object is available.")
}

se <- readr::read_rds("data/se.rds")

# Prepare data for longitudinal analysis ----
# Filter samples with visit quarter information
se <- se[, !is.na(se$visit_quarter)]

# Identify subjects with multiple time points (longitudinal data)
longitudinal_subject <-
  colData(se) %>% 
  as_tibble() %>%
  select(subject, visit) %>% 
  distinct() %>%
  group_by(subject) %>%
  summarise(n = n()) %>% 
  filter(n > 1) %>% 
  pull(subject)

# Keep only longitudinal subjects
se <- se %>% filter(subject %in% longitudinal_subject)

# Perform differential expression analysis ----
# Define contrasts and design formulas for temporal comparisons
DGE_time_d <-
  tibble(
    C1 = rep(c("LS", "LS", "NL"), 2),  # Lesional, Lesional, Non-lesional
    C2 = rep(c("NL", "NN", "NN"), 2),  # Non-lesional, Neighbor, Neighbor
    t = c(rep("visit", 3), rep("visit_quarter", 3))  # Time variables
  ) %>%
  mutate(
    contrast = paste0(C1, "vs", C2),
    # Full model with interaction between skin type and time
    design_f = ifelse(
      C1 == "LS" & C2 == "NL",
      paste("~", "subject", "+", t, "+", "skin_type"),
      paste("~", "gender", "+", t, "+", "skin_type")
    ),
    design_f = paste(design_f, "+", "skin_type", "*", t),
    # Reduced model for likelihood ratio test
    reduce_f = case_when(
      (C1 == "LS") & (C2 == "NL") ~ paste("~", "subject", "+", t, "+", "skin_type"),
      (C2 == "NN") ~ paste("~", "gender", "+", t, "+", "skin_type"),
      TRUE ~ "NULL"
    ),
    se_exp = map2(C1, C2, ~ se[, se$skin_type %in% c(.x, .y)]),
    deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)),
    # Perform likelihood ratio test
    deseq = map2(
      reduce_f, deseq, 
      ~ DESeq(.y, test = "LRT", reduced = .x %>% as.formula, parallel = TRUE)
    ),
    result_n = map2(
      deseq, t, 
      ~ resultsNames(.x) %>% grep(paste0(.y, "_"), ., value = TRUE)
    )
  )

# Extract condition-specific results (LRT test)
DGE_condition_specific <-
  DGE_time_d %>%
  mutate(result = map(deseq, ~ results(.x) %>% as_tibble(rownames = "gene_name"))) %>%
  select(t, contrast, result)

# Extract condition-specific non-interaction results (with shrinkage)
DGE_condition_specific_non <-
  DGE_time_d %>%
  select(t, contrast, deseq, result_n) %>%
  unnest(cols = "result_n") %>%
  mutate(result = map2(deseq, result_n, ~ lfcShrink(.x, coef = .y, type = "apeglm"))) %>%
  mutate(result = map(result, ~ .x %>% as_tibble(rownames = "gene_name"))) %>%
  select(-deseq)

# Combine results ----
DGE_time <- list(
  condition_specific = DGE_condition_specific,
  condition_specific_non = DGE_condition_specific_non,
  description = paste(
    "condition_specific: DGE on time condition specific (interaction term);",
    "condition_specific_non: DGE on time condition non-specific (main effects)"
  )
)

# Save results ----
output_dir <- "data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

saveRDS(DGE_time, file.path(output_dir, "dge_time.rds"))

cat("Temporal differential gene expression analysis completed successfully.\n")
cat("Results saved to:", file.path(output_dir, "dge_time.rds"), "\n")

