# SPACE R JOB
setwd("/home/projects/ku_00015/people/tuhu/multiomics-ad-phd") # In computerome
# load libraries
pacman::p_load(tibble, dplyr, tidyr, readr, stringr, purrr,
               tidybulk, tidySummarizedExperiment,
               BiocParallel, Biobase,
               DESeq2)
core_n <- future::availableCores()
register(MulticoreParam(ifelse(core_n <= 8, core_n - 2, 20)))
se <- readr::read_rds("data/se.rds")

se <- se[,!is.na(se$visit_quarter)]

# filter out those who participated only in longitudinal studies
longitudinal_subject <-
  colData(se) %>% as_tibble() %>%
  select(subject, visit) %>% distinct() %>%
  group_by(subject) %>%
  summarise(n = n()) %>% filter(n>1) %>% pull(subject)

se <- se %>% filter(subject %in% longitudinal_subject)

DGE_time_d <-
  tibble(C1 = rep(c("LS", "LS", "NL"), 2),
         C2 = rep(c("NL", "NN", "NN"), 2),
         t = c(rep("visit", 3), rep("visit_quarter", 3))) %>%
  mutate(contrast = paste0(C1, "vs", C2),
         design_f = ifelse(C1 == "LS" & C2 == "NL",
                           paste("~", "subject", "+", t, "+", "skin_type"),
                           paste("~", "gender", "+", t, "+", "skin_type")),
         design_f = paste(design_f, "+", "skin_type", "*", t),
         reduce_f = case_when(
           (C1 == "LS") & (C2 == "NL")  ~ paste("~", "subject", "+", t, "+", "skin_type"),
           (C2 == "NN") ~ paste("~", "gender", "+", t, "+", "skin_type"),
           TRUE ~ "NULL"),
         se_exp = map2(C1, C2, ~ se[,se$skin_type %in% c(.x, .y)]),
         deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)),
         deseq = map2(reduce_f, deseq, ~ DESeq(.y, test = "LRT", reduced = .x %>% as.formula, parallel = T)),
         result_n = map2(deseq, t, ~ resultsNames(.x) %>% grep(paste0(.y, "_"), ., value = T)))

DGE_condition_specific <-
  DGE_time_d %>%
  mutate(result = map(deseq, ~ results(.x) %>% as_tibble(rownames = "gene_name"))) %>%
  select(t, contrast, result)

DGE_condition_specific_non <-
  DGE_time_d %>%
  select(t, contrast, deseq, result_n) %>%
  unnest(cols = "result_n") %>%
  mutate(result = map2(deseq, result_n, ~ lfcShrink(.x, coef = .y, type = "apeglm"))) %>%
  mutate(result = map(result, ~ .x %>% as_tibble(rownames = "gene_name"))) %>%
  select(-deseq)


# execuate jobs
DGE_time <-
  list(DGE_condition_specific,
       DGE_condition_specific_non)
DGE_time$des <- "DGE_condition_specific: dge on time condition specific; DGE_condition_specific_non: dge on time condition non-specific"
saveRDS(DGE_time, "data/dge_time.rds")
