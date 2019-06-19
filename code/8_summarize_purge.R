

get_purged_umi_counts <- function(called_cells_sample,
                                  umi_count_matrix_sample) {
  umi_purged <- umi_count_matrix_sample$purged
  keep <- drop_na(called_cells_sample)$is_cell_purged
  umi_purged <- umi_purged[, keep]
  
  return(umi_purged)
}


save_purged_umi_counts <- function(umi_purged,
                                   sample_name,
                                   output_dir) {
  purged_data_filepath <-
    file.path(
      output_dir,
      sprintf(
        "%s_umi_count_matrix_purged.rds",
        sample_name
      )
    )
  
  saveRDS(umi_purged, purged_data_filepath)
}



save_purged_umi_counts_all_samples <- function(called_cells, umi_count_matrices, output_dir) {
  purged_umi_counts <- map2(
    called_cells,
    umi_count_matrices,
    get_purged_umi_counts
  )
  
  imap(purged_umi_counts,
       save_purged_umi_counts,
       output_dir = output_dir
  )
}



get_umi_counts_cell <- function(sample_called_cells, sample_umi_counts_cell_gene) {
  sample_umi_counts_cell_gene <- left_join(sample_called_cells %>%
    select(-is_cell_purged),
  sample_umi_counts_cell_gene,
  by = "cell"
  )

  sample_umi_counts_cell <-
    sample_umi_counts_cell_gene %>%
    select(-c("gene")) %>%
    group_by(is_cell_unpurged, cell) %>%
    summarize(
      m = sum(m),
      rm_ret = sum(rm_ret),
      pm = sum(pm),
      rm_disc = sum(rm_disc)
    ) %>%
    arrange(-pm) %>%
    ungroup()

  sample_umi_counts_cell <-
    split(
      sample_umi_counts_cell %>%
        select(-is_cell_unpurged),
      sample_umi_counts_cell$is_cell_unpurged
    ) %>%
    set_names(c("background_cells", "called_cells"))


  return(sample_umi_counts_cell)
}



get_umi_counts_sample <- function(umi_counts_cell) {
  umi_counts_sample <-
    umi_counts_cell %>%
    ungroup() %>%
    select(c("m", "rm_ret", "rm_disc", "pm")) %>%
    summarize(
      m = sum(m),
      rm_ret = sum(rm_ret),
      rm_disc = sum(rm_disc),
      pm = sum(pm)
    )
}

update_summary_stats <- function(summary_stats, umi_counts_sample) {
  umi_counts_sample_all <-
    umi_counts_sample %>%
    select(-split) %>%
    group_by(sample) %>%
    summarize_all(list(~ sum(.))) %>%
    mutate(
      prm_ret = rm_ret / m,
      prm_disc = rm_disc / m,
      ppm = pm / m
    ) %>%
    ungroup() %>%
    mutate(prop_m = m / sum(m))

  umi_counts_sample_called_cells <-
    umi_counts_sample %>%
    filter(split == "called_cells") %>%
    select(-split) %>%
    group_by(sample) %>%
    summarize_all(list(~ sum(.))) %>%
    mutate(
      prm_ret = rm_ret / m,
      prm_disc = rm_disc / m,
      ppm = pm / m
    ) %>%
    ungroup() %>%
    mutate(prop_m = m / sum(m))

  summary_stats$marginal_called_cells <-
    left_join(umi_counts_sample_called_cells,
      summary_stats$marginal,
      by = "sample"
    ) %>%
    mutate(
      FRM = n_reads / m,
      m = m / 1000000
    ) %>%
    select(sample, m, prm_ret, prm_disc, ppm, prop_m, prop_reads, FRM, rm_ret, rm_disc, pm)

  summary_stats$marginal <-
    left_join(umi_counts_sample_all,
      summary_stats$marginal,
      by = "sample"
    ) %>%
    mutate(
      FRM = n_reads / m,
      m = m / 1000000
    ) %>%
    select(sample, m, prm_ret, prm_disc, ppm, prop_m, prop_reads, FRM, rm_ret, rm_disc, pm)

  return(summary_stats)
}



save_phantom_purged_data <- function(read_counts,
                                     phantom_purged_data_filepath,
                                     sample_names,
                                     S) {
  read_counts <-
    read_counts %>%
    filter(keep == TRUE) %>%
    select(c("cell", "umi", "gene", paste0("t", 1:S))) %>%
    setNames(c("cell", "umi", "gene", sample_names))



  saveRDS(read_counts, phantom_purged_data_filepath)
}