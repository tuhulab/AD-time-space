# Output data for CytoSig

library(dplyr)
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193309/suppl/GSE193309%5Fcount%5Fmatrix%5Fnormalized%2Ecsv%2Egz",
              "_tmp/count_matrix_geo.csv.gz")

count_matrix <- gzfile("_tmp/count_matrix_geo.csv.gz")
count_matrix_csv <- readr::read_csv(count_matrix)

count_matrix_csv_log2 <-
  count_matrix_csv %>% mutate(across(-1, ~ log2(.x + .5) %>% round(2) ))

count_matrix_csv_log2[,1:40] %>% write.csv("_tmp/count_matrix_log2_1_40.csv")


