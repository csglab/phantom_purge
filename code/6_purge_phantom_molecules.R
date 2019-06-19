



dedup_reads <- function(read_counts, sample_names) {
  read_counts <-
    read_counts %>%
    mutate_at(
      c(sample_names, paste0(sample_names, "_hat")),
      list(~ as.integer(. > 0))
    )
  
  return(read_counts)
}


create_umi_counts <- function(read_counts,
                              sample_names) {
  
  S <- length(sample_names)
  
  read_counts <- dedup_reads(read_counts, sample_names)
  
  read_counts <- 
    read_counts %>%
    rename_at(sample_names, list(~ paste0("m", 1:S)))%>%
    rename_at(paste0(sample_names, "_hat"), list(~ paste0("t", 1:S)))
    

  # tallying molecule counts by cell-barcode and gene ID
  setDT(read_counts)
  umi_counts_cell_gene <-
    read_counts[, lapply(.SD, sum),
      keyby = "cell,gene,retain",
      .SDcols = -c("umi")
    ]


  # tranforming cell-gene molecule tally table into long format

  umi_counts_cell_gene <-
    melt(
      umi_counts_cell_gene,
      id = 1:3,
      variable.name = "sample",
      measure = patterns(m = "^m", t = "^t"),
      variable.factor = FALSE
    )

  sample_key <-
    tibble(
      sample = as.character(1:S),
      sample_name = sample_names
    )

  umi_counts_cell_gene <-
    left_join(umi_counts_cell_gene,
      sample_key,
      by = "sample"
    )

  umi_counts_cell_gene <-
    umi_counts_cell_gene %>%
    mutate(
      sample = sample_name,
      pm = m - t,
      rm_disc = t * (1 - retain), # discarded real
      rm_ret = t - rm_disc, # retained
      t = NULL,
      sample_name = NULL,
      retain = NULL
    )



  umi_counts_cell_gene <- split(
    umi_counts_cell_gene,
    umi_counts_cell_gene$sample
  )


  return(umi_counts_cell_gene)
}