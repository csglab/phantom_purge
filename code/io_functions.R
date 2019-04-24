
get_h5_filenames <- function(input_dir) {
  metadata <- list.files(path = input_dir,
                         pattern = "h5",
                         recursive = TRUE) %>%
    enframe(name = NULL, value ="sample_name_ext") %>%
    mutate(sample_name=  tools::file_path_sans_ext(sample_name_ext))  %>%
    mutate(filepath = file.path(input_dir, sample_name_ext)) 
  
  molecule_info_filepaths <- metadata$filepath
  names(molecule_info_filepaths) <- metadata$sample_name
  
  return(molecule_info_filepaths)
  
}

read10xMolInfo_3 <- function(sample)
  # Modified function from dropUtils
{
  data <- list()
  
  all.barcodes <- as.vector(h5read(sample, "/barcodes"))
  all.barcodes <-
    sub("-[0-9]+", "", all.barcodes) # removing GEM group.
  data$cell <-
    all.barcodes[as.vector(h5read(sample, "/barcode_idx")) + 1L]
  data$umi <- as.vector(h5read(sample, "/umi"))
  data$gene <- as.vector(h5read(sample, "/feature_idx")) + 1L
  data$reads <- as.vector(h5read(sample, "/count"))
  gene_name <- h5read(sample, "/features/name")
  gene_ids <- h5read(sample, "/features/id")
  
  nonunique <- 
    duplicated(gene_name) | 
    duplicated(gene_name, fromLast = TRUE)
  
  gene_name[nonunique] <- str_c(gene_name[nonunique], 
                                gene_ids[nonunique],
                                sep=":")
  
  data <- do.call(tibble, data)
  
  
  data <-
    data %>%
    filter(gene <= length(gene_name)) %>%
    mutate(gene= gene_name[gene] )
  
  # Don't define the total cell pool here, as higher level functions may want to use gem_group.
  return(data = data)
}


read10xMolInfo_2 <- function(sample)
  # Modified function from dropUtils
{
  data <- list()
  
  data$cell <-
    .Call(DropletUtils:::cxx_get_cell_barcodes,
          sample,
          "barcode",
          NULL)
  data$umi <- as.vector(h5read(sample, "/umi"))
  data$gene <- as.vector(h5read(sample, "/gene")) + 1L
  data$reads <- as.vector(h5read(sample, "/reads"))
  gene_name <- h5read(sample, "/gene_names")
  gene_ids <- h5read(sample, "/gene_ids")
  
  nonunique <- 
    duplicated(gene_name) | 
    duplicated(gene_name, fromLast = TRUE)
  
  gene_name[nonunique] <- str_c(gene_name[nonunique], 
                                gene_ids[nonunique],
                                sep=":")
  
  data <- do.call(tibble, data)
  
  data <-
    data %>%
    filter(gene <= length(gene_name)) %>%
    mutate(gene= gene_name[gene] )
  
  return(data)
}

load_molecule_info_data <- function(input_dir) {
  
  molecule_info_filepaths <- get_h5_filenames(input_dir)
  
  
  available <- h5ls(molecule_info_filepaths[[1]], recursive = FALSE)
  version <- if ("barcode_idx" %in% available$name)
    "3"
  else
    "2"
  
  if (version == "3") {
    read_counts <- future_map(molecule_info_filepaths,
                              read10xMolInfo_3)
    
  }
  else {
    read_counts <- map(molecule_info_filepaths,
                       read10xMolInfo_2)
    
  }
  
  return(read_counts)
}


save_phantom_purged_data <- function(read_counts,
                                     phantom_purged_data_filepath,
                                     sample_names,
                                     S) {
  
  read_counts <-
    read_counts %>% 
    filter(keep==TRUE)%>%
    select(c("cell", "umi", "gene", paste0("t", 1:S))) %>%
    setNames(c("cell", "umi", "gene", sample_names))
  
  
  
  saveRDS(read_counts, phantom_purged_data_filepath)
  
}


save_purged_umi_counts <- function(umi_purged,
                                   sample_name,
                                   output_dir
                                   ) {
  
  
  purged_data_filepath <- 
    file.path(output_dir,
              sprintf("%s_umi_count_matrix_purged.rds",
                      sample_name))

  saveRDS(umi_purged, purged_data_filepath)
  
}
# save_reassigned_read_counts <- function(read_counts,
#                                         data_outcome,
#                                         S,
#                                         reassigned_counts_filepath,
#                                         sample_names) {
#   data_outcome <-
#     data_outcome  %>%
#     select(c("outcome", paste0("t", 1:S))) %>%
#     setNames(c("outcome", sample_names))
#   
#   reassigned_read_counts <-
#     join_data(read_counts %>%
#                 select(outcome, cell, umi, gene),
#               data_outcome)
#   
#   saveRDS(reassigned_read_counts, reassigned_counts_filepath)
#   
# }