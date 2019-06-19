get_h5_filenames <- function(input_dir) {
  metadata <-
    list.files(
      path = input_dir,
      pattern = "h5",
      recursive = TRUE
    ) %>%
    enframe(name = NULL, value = "sample_name_ext") %>%
    mutate(sample_name = tools::file_path_sans_ext(sample_name_ext)) %>%
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

  gene_name[nonunique] <-
    str_c(gene_name[nonunique],
          gene_ids[nonunique],
          sep = ":")

  data <- do.call(tibble, data)

  data <-
    data %>%
    filter(gene <= length(gene_name)) %>%
    mutate(gene = gene_name[gene])

  return(data = data)
}


read10xMolInfo_2 <- function(sample)
# Modified function from dropUtils
{
  data <- list()

  data$cell <-
    .Call(
      DropletUtils:::cxx_get_cell_barcodes,
      sample,
      "barcode",
      NULL
    )
  data$umi <- as.vector(h5read(sample, "/umi"))
  data$gene <- as.vector(h5read(sample, "/gene")) + 1L
  data$reads <- as.vector(h5read(sample, "/reads"))
  gene_name <- h5read(sample, "/gene_names")
  gene_ids <- h5read(sample, "/gene_ids")

  nonunique <-
    duplicated(gene_name) |
      duplicated(gene_name, fromLast = TRUE)

  gene_name[nonunique] <-
    str_c(gene_name[nonunique],
      gene_ids[nonunique],
      sep = ":"
    )

  data <- do.call(tibble, data)

  data <-
    data %>%
    filter(gene <= length(gene_name)) %>%
    mutate(gene = gene_name[gene])

  return(data)
}

load_molecule_info_data <- function(input_dir) {
  molecule_info_filepaths <- get_h5_filenames(input_dir)


  available <-
    h5ls(molecule_info_filepaths[[1]],
      recursive = FALSE
    )
  version <- if ("barcode_idx" %in% available$name) {
    "3"
  } else {
    "2"
  }

  if (version == "3") {
    read_counts <- future_map(
      molecule_info_filepaths,
      read10xMolInfo_3
    )
  }
  else {
    read_counts <- map(
      molecule_info_filepaths,
      read10xMolInfo_2
    )
  }

  return(read_counts)
}


rename_var_data_list <- function(read_counts) {
  sample_names <- names(read_counts)

  name_list <- function(x, y) {
    names(x[[1]]) <- y
    return(x[[1]])
  }

  reads_vars <-
    rep("reads", length(sample_names)) %>%
    set_names(sample_names) %>%
    map(as.list) %>%
    map2(., sample_names, name_list)

  rename_var <- function(x, y) {
    x <-
      x %>%
      rename(!!!y)
    return(x)
  }

  read_counts <-
    read_counts %>%
    map2(., reads_vars, rename_var)
}

###########################################################

join_read_counts <- function(read_counts, keys = c("cell", "gene", "umi")) {
  # Join transciptome mapped data by cell, umi, and gene
  read_counts <-
    read_counts %>%
    rename_var_data_list() %>%
    map(setDT) %>%
    reduce(merge,
      all = TRUE,
      sort = FALSE,
      no.dups = TRUE,
      by = keys
    ) %>%
    replace(is.na(.), 0)

  return(read_counts)
}


add_outcome_variable <- function(read_counts, sample_names) {
  setDT(read_counts)

  read_counts[,
    outcome := do.call(paste, c(.SD, sep = ",")),
    .SDcols = sample_names
  ]

  read_counts[order(outcome)]

  return(read_counts)
}

########################################################

create_joined_counts <- function(input_dir, read_counts_filepath=NULL) {
  
  if (!is.null(read_counts_filepath) & file.exists(read_counts_filepath)) {
    read_counts <- readRDS(read_counts_filepath)
  } else if (is.null(read_counts_filepath)) {
    read_counts <- load_molecule_info_data(input_dir)
    read_counts <- join_read_counts(read_counts)
    S <- ncol(read_counts) - 3
    sample_names <- colnames(read_counts)[4:(S + 3)]
    read_counts <- add_outcome_variable(read_counts, sample_names)}
  
  if(!is.null(read_counts_filepath) & !file.exists(read_counts_filepath)){
    saveRDS(read_counts, read_counts_filepath)
  }
  
  return(read_counts)
}