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
         deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)))

deseq_r <- lapply(seq_along(rownames(DGE_time_d)),
                function(n){
                  deseq <- DESeq(DGE_time_d$deseq[[n]],
                                 test = "LRT",
                                 reduced = as.formula(DGE_time_d$reduce_f[[n]]),
                                 parallel = T)
                  return(deseq)
                })

result_n <- lapply(seq_along(deseq_r),
                   function(n){
                     rsname <- resultsNames(deseq_r[[n]]) %>% grep(paste0(DGE_time_d$t[[n]], "_"), ., value = T)
                     return(rsname)
                   })

DGE_condition_specific <-
  lapply(seq_along(deseq_r),
         function(n){
           rs <- results(deseq_r[[n]]) %>% as_tibble(rownames = "gene_name")
         })

DGE_condition_rs_specific <-
  DGE_time_d %>% select(t, contrast) %>% mutate(DGE_condition_specific)

DGE_condition_rs_specific_non <-
  DGE_time_d %>% mutate(result_n, deseq_r) %>%
  unnest(cols = "result_n") %>% select(t, contrast, deseq_r, result_n)

rs_specific_non <-
  lapply(seq_along(rownames(DGE_condition_rs_specific_non)),
       function(n){
         shrink_rs <-
           lfcShrink(DGE_condition_rs_specific_non$deseq_r[[n]],
                     coef = DGE_condition_rs_specific_non$result_n[[n]],
                     type = "apeglm")
         shrink_rs <-
           as_tibble(shrink_rs, rownames = "gene_name")
         return(shrink_rs)
       })
DGE_condition_rs_specific_non <-
  DGE_condition_rs_specific_non %>%
  mutate(rs_specific_non) %>%
  select(-deseq_r)

# execuate jobs
DGE_time <-
  list(DGE_condition_rs_specific,
       DGE_condition_rs_specific_non)
DGE_time$des <- "DGE_condition_rs_specific: dge on time condition specific; DGE_condition_rs_specific_non: dge on time condition non-specific; non-maps"
saveRDS(DGE_time, "data/dge_time_nmap.rds")
