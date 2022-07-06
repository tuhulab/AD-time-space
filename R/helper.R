pull_full_data <- function(dl_path = ...) {
    dl <- readr::read_rds(dl_path)
    full_data <- bind_cols(dl$gene_annotation %>%
        select(Geneid), dl$countable) %>%
        tidyr::pivot_longer(!Geneid, names_to = "BAM_ID", values_to = "VSN_expr") %>%
        left_join(dl$metadata)
    return(full_data)
}

counttable_merge_library_fun <- function(counttable_data = ..., lib_to_merge_vector = ...) {
    lib_id <- counttable_data %>%
        colnames() %>%
        str_extract("lib\\d{1,}")
    merged_counttable <- sapply(lib_to_merge_vector, function(one_lib_id_to_merge) {
        merged_counts <- counttable_data %>%
            select((lib_id == one_lib_id_to_merge) %>%
                which()) %>%
            rowSums()
        merged_counts_df <- tibble(one_lib_id_to_merge = merged_counts)
        return(merged_counts_df)
    })
    # The function to merge libs for counttable ----------------- list tidy
    MergedLib <- do.call(rbind.data.frame, merged_counttable) %>%
        t() %>%
        as.data.frame() %>%
        tibble()
    rownames(MergedLib) <- c()
    # new colnames
    colnames(MergedLib) <- paste0(colnames(counttable_data)[which(lib_id %in% lib_to_merge_vector)] %>%
        str_extract("NG[:graph:]{1,}_lib\\d{1,}") %>%
        unique(), "_merged")
    # remove the old libs
    index_of_oldLibs <- which(lib_id %in% lib_to_merge_vector)
    counttable_data_rmOldLib <- counttable_data %>%
        select(-all_of(index_of_oldLibs))
    # merge with new libs
    counttable_data_rmOldLib_addMergedLib <- counttable_data_rmOldLib %>%
        bind_cols(MergedLib)
    return(counttable_data_rmOldLib_addMergedLib)
}
