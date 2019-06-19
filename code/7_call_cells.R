call_emptydrops <- function(sample_umi_count, lower = 200, fdr_thresh = 0.005) {
  sample_emptydrops <-
    DropletUtils::emptyDrops(sample_umi_count,
      lower = lower
    ) %>%
    as.data.frame() %>%
    rownames_to_column(var = "cell") %>%
    select(cell, FDR)

  sample_emptydrops[is.na(sample_emptydrops)] <- 1

  sample_emptydrops <-
    sample_emptydrops %>%
    mutate(
      is_cell = FDR <= fdr_thresh,
      FDR = NULL
    )

  return(sample_emptydrops)
}

call_defaultdrops <- function(sample_umi_count) {
  sample_defaultdrops <- DropletUtils::defaultDrops(sample_umi_count)
  sample_defaultdrops <- enframe(sample_defaultdrops,
    name = "cell",
    value = "is_cell"
  )

  return(sample_defaultdrops)
}

call_cells <- function(umi_count_matrix, sample_name) {
  called_cells <- future_map(
    umi_count_matrix,
    plyr::failwith(
      NULL,
      call_emptydrops
    )
  )


  if (is.null(called_cells$unpurged) | is.null(called_cells$purged)) {
    message(paste0(
      "Using defaultDrops cell calling instead.
      EmptyDrops function call failed for purged sample ",
      sample_name
    ))
    called_cells <- map(umi_count_matrix, call_defaultdrops)
  }


  called_cells <- left_join(called_cells$unpurged,
    called_cells$purged,
    by = "cell",
    suffix = c("_unpurged", "_purged")
  )


  return(called_cells)
}

make_umi_count_matrices <- function(umi_counts_cell_gene) {
  genes <- unique(umi_counts_cell_gene$gene)
  
  sample_data <-
    umi_counts_cell_gene %>%
    filter(m > 0)
  
  
  unpurged <- DropletUtils::makeCountMatrix(sample_data$gene,
                                            sample_data$cell,
                                            all.genes = genes,
                                            value = sample_data$m
  )
  
  sample_data <-
    sample_data %>%
    filter(rm_ret > 0)
  
  
  purged <- DropletUtils::makeCountMatrix(sample_data$gene,
                                          sample_data$cell,
                                          all.genes = genes,
                                          value = sample_data$rm_ret
  )
  
  return(list(
    unpurged = unpurged,
    purged = purged
  ))
}

get_cells_tally <- function(called_cells, sample_name) {

  # keys <- c(`FALSE`="background", `TRUE`="cell")

  cells_tally <-
    called_cells %>%
    group_by(
      is_cell_unpurged,
      is_cell_purged
    ) %>%
    tally() %>%
    ungroup() %>%
    complete(is_cell_unpurged, is_cell_purged, fill = list(n = 0)) %>%
    mutate(
      # is_cell_purged= recode(as.character(is_cell_purged), !!!keys),
      # is_cell_unpurged= recode(as.character(is_cell_unpurged), !!!keys),
      barcode = c(
        "consensus_background",
        "transition_cell",
        "phantom_background",
        "transition_background",
        "consensus_cell",
        "phantom_cell"
      )
    ) %>%
    select(barcode, n) %>%
    group_by(barcode) %>%
    set_names("barcode", sample_name)

  return(cells_tally)
}

call_cells_all_samples <- function(umi_counts_cell_gene, output_dir) {
  
  
  umi_count_matrices <- map(
    umi_counts_cell_gene,
    make_umi_count_matrices
  )
  
  called_cells <- imap(
    umi_count_matrices,
    call_cells
  )
  
  # Save purged umi counts
  
  save_purged_umi_counts_all_samples(
    called_cells,
    umi_count_matrices,
    output_dir
  )


  called_cells_tally <-
    imap_dfc(called_cells, get_cells_tally) %>%
    select(c("barcode", names(called_cells)))


  return(list(called_cells = called_cells, called_cells_tally = called_cells_tally))
}