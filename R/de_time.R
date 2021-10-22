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

se <- se[,!is.na(se$visit_quarter)]

# execuate jobs
DGE_time <-
  tibble(
    C1 = rep(c("LS", "LS", "NL"), 2),
    C2 = rep(c("NL", "NN", "NN"), 2),
    t = c(rep("visit", 3), rep("visit_quarter", 3))
  ) %>%
  mutate(
    contrast = paste0(C1, "vs", C2),
    se_exp = map2(C1, C2, ~ se[,se$skin_type %in% c(.x, .y)]),
    design_f = ifelse(C1 == "LS" & C2 == "NL",
                      paste("~", "subject", "+", t, "+", "skin_type", "+", "skin_type", "*", t),
                      paste("~", "gender", "+", t, "+", "skin_type", "+", "skin_type", "*", t)),
    reduce_f = ifelse(C1 == "LS" & C2 == "NL",
                      paste("~", "subject", "+", t, "+", "skin_type"),
                      paste("~", "gender", "+", t, "+", "skin_type")),
    deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)),
    deseq = map2(deseq, reduce_f, ~ DESeq(.x, test = "LRT", reduced = .y %>% as.formula, parallel = T)),
    res = map(deseq, ~ DESeq2::results(.x)),
    res = map(res, ~ .x %>% as.data.frame %>% as_tibble(rownames = "gene_name"))
  )
saveRDS(DGE_time %>% select(C1:contrast, res), "data/dge_time.rds")
