#' Pull and Merge Full Dataset
#'
#' This function reads an RDS file containing gene annotation, count data, and metadata,
#' then merges them into a single long-format data frame.
#'
#' @param dl_path Character string specifying the path to the RDS file containing
#'   the data list with gene_annotation, countable, and metadata components.
#'
#' @return A tibble in long format containing gene IDs, BAM IDs, VSN expression values,
#'   and associated metadata for each sample.
#'
#' @examples
#' \dontrun{
#' full_data <- pull_full_data("data/processed_data.rds")
#' }
#'
#' @export
pull_full_data <- function(dl_path = ...) {
    if (!file.exists(dl_path)) {
        stop("File does not exist: ", dl_path)
    }
    
    dl <- readr::read_rds(dl_path)
    
    # Validate input structure
    required_cols <- c("gene_annotation", "countable", "metadata")
    if (!all(required_cols %in% names(dl))) {
        stop("Input RDS must contain: ", paste(required_cols, collapse = ", "))
    }
    
    full_data <- bind_cols(
        dl$gene_annotation %>% select(Geneid), 
        dl$countable
    ) %>%
        tidyr::pivot_longer(!Geneid, names_to = "BAM_ID", values_to = "VSN_expr") %>%
        left_join(dl$metadata, by = "BAM_ID")
    
    return(full_data)
}

#' Merge Technical Replicates in Count Table
#'
#' This function merges technical replicates (libraries with the same ID) by summing
#' their counts and creates a new count table with merged libraries.
#'
#' @param counttable_data Data frame containing count data with library IDs in column names
#' @param lib_to_merge_vector Character vector of library IDs to merge
#'
#' @return A data frame with merged library counts. Original libraries are removed and
#'   replaced with merged versions. Column names for merged libraries end with "_merged".
#'
#' @details
#' The function identifies all columns with matching library IDs, sums their counts row-wise,
#' and creates new columns with merged data. The original columns are removed.
#'
#' @examples
#' \dontrun{
#' # Merge libraries lib01 and lib02
#' merged_counts <- counttable_merge_library_fun(
#'   counttable_data = raw_counts,
#'   lib_to_merge_vector = c("lib01", "lib02")
#' )
#' }
#'
#' @export
counttable_merge_library_fun <- function(counttable_data = ..., lib_to_merge_vector = ...) {
    # Input validation
    if (!is.data.frame(counttable_data)) {
        stop("counttable_data must be a data frame")
    }
    
    if (length(lib_to_merge_vector) == 0) {
        warning("No libraries to merge, returning original data")
        return(counttable_data)
    }
    
    # Extract library IDs from column names
    lib_id <- counttable_data %>%
        colnames() %>%
        str_extract("lib\\d{1,}")
    
    # Merge counts for each library ID
    merged_counttable <- sapply(lib_to_merge_vector, function(one_lib_id_to_merge) {
        # Find columns matching this library ID
        cols_to_merge <- which(lib_id == one_lib_id_to_merge)
        
        if (length(cols_to_merge) == 0) {
            warning("Library ID ", one_lib_id_to_merge, " not found in data")
            return(NULL)
        }
        
        # Sum counts across technical replicates
        merged_counts <- counttable_data %>%
            select(all_of(cols_to_merge)) %>%
            rowSums()
        
        merged_counts_df <- tibble(one_lib_id_to_merge = merged_counts)
        return(merged_counts_df)
    })
    
    # Convert list to data frame
    MergedLib <- do.call(rbind.data.frame, merged_counttable) %>%
        t() %>%
        as.data.frame() %>%
        tibble()
    rownames(MergedLib) <- c()
    
    # Create new column names for merged libraries
    colnames(MergedLib) <- paste0(
        colnames(counttable_data)[which(lib_id %in% lib_to_merge_vector)] %>%
            str_extract("NG[:graph:]{1,}_lib\\d{1,}") %>%
            unique(), 
        "_merged"
    )
    
    # Remove original libraries and add merged ones
    index_of_oldLibs <- which(lib_id %in% lib_to_merge_vector)
    counttable_data_rmOldLib <- counttable_data %>%
        select(-all_of(index_of_oldLibs))
    
    counttable_data_rmOldLib_addMergedLib <- counttable_data_rmOldLib %>%
        bind_cols(MergedLib)
    
    return(counttable_data_rmOldLib_addMergedLib)
}
