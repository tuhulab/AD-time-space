# SPACE R JOB
setwd("/home/projects/ku_00015/people/tuhu/multiomics-ad-phd") # In computerome
# load libraries
pacman::p_load(tibble, dplyr, tidyr, readr, stringr, purrr,
               tidybulk, tidySummarizedExperiment,
               BiocParallel, Biobase,
               DESeq2)
core_n <- future::availableCores()
register(MulticoreParam(ifelse(core_n <= 8, core_n - 2, core_n / 2)))
se <- readr::read_rds("data/se.rds")

se <- se[,!is.na(se$visit_quarter)]

# filter out those who participated only in longitudinal studies
# longitudinal_subject <-
#   colData(se) %>% as_tibble() %>%
#   select(subject, visit) %>% distinct() %>%
#   group_by(subject) %>%
#   summarise(n = n()) %>% filter(n>1) %>% pull(subject)

# se <- se %>% filter(subject %in% longitudinal_subject)

DGE_time_d <-
  map_dfr(seq_len(2), ~ tibble(
          C1 = rep(c("LS", "LS", "NL"), 2),
          C2 = rep(c("NL", "NN", "NN"), 2),
          t = c(rep("visit", 3), rep("visit_quarter", 3))
          )) %>%
  mutate(condition_specific = c(rep(TRUE, 6), rep(FALSE, 6)),
         contrast = paste0(C1, "vs", C2),
         # se_exp = map2(C1, C2, ~ se[,se$skin_type %in% c(.x, .y)]),
         design_f = ifelse(C1 == "LS" & C2 == "NL",
                           paste("~", "subject", "+", t, "+", "skin_type"),
                           paste("~", "gender", "+", t, "+", "skin_type")),
         design_f = ifelse(condition_specific == TRUE,
                           paste(design_f, "+", "skin_type", "*", t),
                           design_f),
         reduce_f = case_when(
           (C1 == "LS") & (C2 == "NL") & (condition_specific == TRUE) ~ paste("~", "subject", "+", t, "+", "skin_type"),
           (C2 == "NN") & (condition_specific == TRUE)  ~ paste("~", "gender", "+", t, "+", "skin_type"),
           TRUE ~ "NULL"),
         se_exp = map2(C1, C2, ~ se[,se$skin_type %in% c(.x, .y)]),
         deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)))

DGE_time_d_t <-
  DGE_time_d %>%
  mutate(
    deseq = pmap(list(condition_specific, reduce_f, deseq),
                 function(x, y, z){
                   if(x == TRUE){
                     DESeq_res <- DESeq(z, test = "LRT", reduced = y %>% as.formula(), parallel = TRUE)
                   }
                   else{
                     DESeq_res <- DESeq(z)
                   }
                   return(DESeq_res)
                 }))





# execuate jobs
DGE_time <-
  DGE_time_d %>%
  mutate(
    contrast = paste0(C1, "vs", C2),
    se_exp = map2(C1, C2, ~ se[,se$skin_type %in% c(.x, .y)]),
    design_f = ifelse(C1 == "LS" & C2 == "NL",
                      paste("~", "subject", "+", t, "+", "skin_type", "+", "skin_type", "*", t),
                      paste("~", "gender", "+", t, "+", "skin_type", "+", "skin_type", "*", t)),
    reduce_f = ifelse(C1 == "LS" & C2 == "NL",
                      paste("~", "subject", "+", t, "+", "skin_type"),
                      paste("~", "gender", "+", t, "+", "skin_type")),
    deseq = map2(se_exp, design_f, ~ DESeqDataSet(.x, .y %>% as.formula)))

deseq <-
  lapply(1:nrow(DGE_time), function(x){
    deseq <- DESeq(DGE_time$deseq[[x]],
                   test = "LRT",
                   reduced = DGE_time$reduce_f[[x]] %>% as.formula,
                   parallel = T)
    return(deseq)
         })

res <-
  lapply(1:length(deseq), function(x){
    res <- results(deseq[[x]])
    })

res <-
  lapply(1:length(res), function(x){
    res <- res[[x]] %>% as.data.frame() %>% as_tibble(rownames = "gene_name")
  })

DGE_time <-
  DGE_time %>% select(C1:contrast) %>%
  mutate(res) %>% unnest() %>%
  filter(padj < .05, abs(log2FoldChange) > 1)

saveRDS(DGE_time, "data/dge_time.rds")
