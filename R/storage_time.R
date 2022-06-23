library(SummarizedExperiment)
library(tidybulk)
library(tidySummarizedExperiment)
library(lmerTest)

se <- readr::read_rds("data/se.rds")

coldata <-
  se %>% colData() %>% as_tibble() %>%
  mutate(date_visit = date_visit %>% as.Date("%d.%m.%Y"),
         time_storage = as.Date("2020-04-01") - date_visit)

se$time_storage <- coldata$time_storage

helen_emma_list <- data.frame(immu_func = "general_inflammation", gene = "MMP12") %>%
  add_row(immu_func = "Th1", gene = c("IFNG", "CXCL9", "CXCL10", "STAT1")) %>%
  add_row(immu_func = "Th2", gene = c("IL4R", "IL13", "IL31", "CCL17", "CCL18", "CCL22","CCL26", "OSM")) %>%
  add_row(immu_func = "Th17", gene = c("IL17A", "IL20", "PI3", "IL36G","CCL20", "CXCL1")) %>%
  add_row(immu_func = "Th17_Th22", gene = c("S100A7", "S100A8", "S100A9", "S100A12")) %>%
  add_row(immu_func = "Th22", gene = c("IL22", "IL32")) %>%
  add_row(immu_func = "Treg", gene = c("IL10", "FOXP3", "CTLA4")) %>%
  add_row(immu_func = "negative_regulators", gene=c("IL34", "IL37", "LORICRIN"))

se_panel <-
  se[intersect(rownames(se), helen_emma_list$gene), ] %>% as_tibble()


lme_m <-
  lapply(se_panel$.feature %>% unique,
       function(g){
         data <- se_panel %>% filter(.feature == g)
         model <- lmer(log2(counts_scaled + 1) ~ skin_type + time_storage + (1|subject),
                       data = data) %>%
           broom.mixed::tidy()
       })

lme_m %>% purrr::reduce(bind_rows) %>%
  filter(term == "time_storage") %>%
  mutate(padj = p.adjust(p.value)) %>%
  summarise(median_effect = median(estimate),
            median_p = median(p.value))
