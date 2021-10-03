# SPACE R JOB
setwd("/home/projects/ku_00015/people/tuhu/multiomics-ad-phd") # In computerome
# load libraries
pacman::p_load(tibble, dplyr, tidyr, readr, stringr, purrr,
               tidybulk, tidySummarizedExperiment,
               BiocParallel, Biobase,
               DESeq2)
knitr::opts_chunk$set(message = FALSE, warning = FALSE, echo = FALSE, fig.align = 'center')
core_n <- future::availableCores()
register(MulticoreParam(ifelse(core_n <= 8, core_n - 2, core_n -6)))
se <- readr::read_rds("data/se.rds")
# execuate jobs
DGE_space <-
  tibble(
    C1 = c("LS", "LS", "NL"),
    C2 = c("NL", "NN", "NN")
  ) %>%
  mutate(
    contrast = paste0(C1, "vs", C2),
    subject_matched = ifelse(C1 == "LS" & C2 =="NL", TRUE, FALSE),
    se_exp = map2(C1, C2, ~ se[,se$skin_type %in% c(.x, .y)]),
    design_f = ifelse(subject_matched,
                      "~ subject + skin_type + visit + skin_type * visit + biopsy_area",
                      "~ gender + skin_type + visit + skin_type * visit + biopsy_area"),
    deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)))
DGE_space_res <-
  lapply(DGE_space$deseq, DESeq, parallel = T)
res_name <-
  lapply(DGE_space_res,
         function(x){
           resultsNames(x) %>% grep(pattern = "biopsy_area", ., value = T)
           })
DGE_space_res_t <-
  DGE_space %>%
  select(contrast) %>%
  mutate(res_name, deseq_res_id = 1:nrow(.)) %>%
  unnest(res_name)
DGE_space_shrink <-
  lapply(
    1:nrow(DGE_space_res_t),
    function(x){
      res_n <- DGE_space_res_t$res_name[x]
      deseq_res <- DGE_space_res[[DGE_space_res_t$deseq_res_id[x]]]
      lfc_shrink <- lfcShrink(deseq_res, coef = res_n, type = "apeglm", parallel = T)
      lfc_shrink_t <- as_tibble(lfc_shrink, rownames = "gene_name")
    })
DGE_space_shrink_t <-
  DGE_space_res_t %>%
  mutate(DGE_space_shrink) %>%
  select(-deseq_res_id)
saveRDS(DGE_space_shrink_t, "data/dge_space_terminal.rds")
