run_workflow <- function(dataset_name, project_dir, max_fpr) {
  
  tic("Running workflow")
  
  
  data_dir <- file.path(project_dir, "data", dataset_name)
  output_dir <- file.path(data_dir, "output")
  input_dir <- file.path(data_dir, "input")
  
  read_counts_filepath <- file.path(output_dir,
                                    paste0(dataset_name,
                                           "_read_counts.rds"))
  
  
  
  # reassigned_counts_filepath <-
  #   file.path(output_dir,
  #             sprintf("%s_reassigned_read_counts.rds",
  #                     dataset_name))
  

  
  results_filepath <-
    file.path(output_dir,
              sprintf("%s_results.rds",
                      dataset_name))

    
  
  read_counts <- create_joined_read_counts_dt(input_dir,
                                              read_counts_filepath)
 
  
  S <- ncol(read_counts) -4
  sample_names <- colnames(read_counts)[4:(S+3)]
  
  tic("Step 2: creating outcome counts datatable with grouping vars")
  outcome_counts <- create_outcome_counts(read_counts)
  toc()
  
  tic("Step 3: creating a chimera counts datatable and estimating hopping rate")
  
  fit_out <- estimate_hopping_rate(outcome_counts,
                                   r_upper=25)
  
  toc()
  
  phat <- fit_out$glm_estimates$phat

  
  tic("Step 4: computing read counts distribution statistics")
  reads_dist_summary <- get_reads_dist_estimates(outcome_counts,
                                                 phat)
  toc()
  
    
  tic("Step 5: estimating pi_r matrix")
  pi_r_hat <- estimate_pi_r_hat_matrix(reads_dist_summary$conditional,
                                       S,
                                       phat)
  toc()
  
  tic("Step 6: infering the true sample of origin")
  outcome_counts <- infer_sample_of_origin( outcome_counts,
                                            pi_r_hat,
                                            S,
                                            phat)
  toc()
  
  
  tic("Step 7: estimating g and computing classification metrics")
  u <-  
    reads_dist_summary$summary_stats %>%
    pull(u)
  
  g <- get_prop_fugues(reads_dist_summary$conditional$m_bar, 
                       phat)
  
  outcome_counts <- compute_classification_metrics(outcome_counts,
                                                   g,
                                                   u)
  toc()
  
  
  tic("Step 8: determining the optimal cutoff")
  optimal_cutoff <- compute_optimal_cutoff(outcome_counts, u, g)
  toc()
  
  
  tic("Step 9: computing proportion of nonmissingness and updating summary data list")
  
  p_nonmiss_dt <- compute_prop_nonmissingness(outcome_counts)
  reads_dist_summary <- update_reads_dist_summary(reads_dist_summary,
                                                  p_nonmiss_dt,
                                                  g)
    
  toc()
  
  
  score_threshold <- 
    get_q_cutoff(optimal_cutoff,
                  outcome_counts,
                  max_fpr)
  
  
  tic(sprintf("Step 10: purging phantoms at q cutoff of  %f. Max-FPR threshold user-set %s",
              score_threshold,
              !is.null(max_fpr)))
  
  umi_counts_cell_gene <-
    purge_phantoms(
      read_counts,
      outcome_counts,
      score_threshold
    )
  toc()
  
  
  tic("Step 11: calling cells")
  
  
  names(sample_names) <- sample_names
  

  umi_count_matrices <- map(umi_counts_cell_gene, 
                            make_umi_count_matrices)
  
  called_cells <- imap(umi_count_matrices, 
                      call_cells)
  

  
  toc()
  
  
  tic("Step 12: saving purged data")
  
  purged_umi_counts <- map2(called_cells,
                            umi_count_matrices,
                            get_purged_umi_counts)
  
  imap(purged_umi_counts, 
       save_purged_umi_counts, 
       output_dir=output_dir)
  
  toc()



  
  tic("Step 13: tallying molecules by cell-barcode")
  
  called_cells_tally <- 
    imap_dfc(called_cells, get_cells_tally) %>% 
    select(c("barcode", names(called_cells)))
  
  umi_counts_cell <- map2(called_cells,  
                          umi_counts_cell_gene,   
                          get_umi_counts_cell)
  
  
  umi_counts_sample <- 
    map(umi_counts_cell,
        map_dfr,
        get_umi_counts_sample,
                           .id="split") %>%
    bind_rows( .id="sample")
  
  reads_dist_summary <- update_reads_dist_summary_2(reads_dist_summary,
                                                    umi_counts_sample)
  toc()
  
  data_list <- list(
    outcome_counts = outcome_counts  %>%
      select(outcome, n, q, qs,j, o, FPR, FNR,  everything()),
    reads_dist_summary = reads_dist_summary,
    fit_out = fit_out,
    optimal_cutoff = optimal_cutoff,
    pi_r_hat=pi_r_hat,
    umi_counts_cell=umi_counts_cell,
    called_cells_tally=called_cells_tally)

  tic("Step 14: saving results")

  saveRDS(data_list,
          results_filepath)
  toc()
  
  data_list <- c(data_list, list(read_counts=read_counts))
  
  toc()
  
  return(data_list)
}