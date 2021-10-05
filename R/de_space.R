# SPACE R JOB
setwd("/home/projects/ku_00015/people/tuhu/multiomics-ad-phd") # In computerome
# load libraries
pacman::p_load(tibble, dplyr, tidyr, readr, stringr, purrr,
               tidybulk, tidySummarizedExperiment,
               BiocParallel, Biobase,
               DESeq2)
knitr::opts_chunk$set(message = FALSE, warning = FALSE, echo = FALSE, fig.align = 'center')
core_n <- future::availableCores()
register(MulticoreParam(ifelse(core_n <= 8, core_n - 2, core_n - 6)))
se <- readr::read_rds("data/se.rds")
# execuate jobs
DGE_space <-
  tibble(
    C1 = c("LS", "NN"),
    C2 = c("NL", "NN")
  ) %>%
  mutate(
    contrast = paste0(C1, "vs", C2),
    se_exp = map2(C1, C2, ~ se[,se$skin_type %in% c(.x, .y)]),
    design_f = ifelse(C1 == "NN" & C2 == "NN",
                      "~ subject + visit + biopsy_area",
                      "~ subject + skin_type + visit + biopsy_area"),
    deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)),
    deseq = map(deseq, ~ DESeq(.x, parallel = T)),
    res_name = map(deseq, ~ resultsNames(.x) %>% grep(pattern = "biopsy_area", ., value = T))) %>%
  unnest(cols = res_name) %>%
  mutate(
    res_shrink = map2(deseq, res_name, ~ lfcShrink(.x, coef = .y, type = "apeglm", parallel = T)),
    res_shrink = map(res_shrink, ~ as_tibble(.x, rownames = "gene_name"))
  )
saveRDS(DGE_space %>% select(contrast, res_name, res_shrink), "data/dge_space.rds")
