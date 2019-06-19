
get_product_logsum <- function(x, y) {
  z_log <- log(x) + log(y)
  z <- exp(z_log)
  return(z)
}

# compute_likelihood <- function(y, S, phat){
#   
#   p_mult_vec <- rep(0, S)
#   for (i in 1:S) {
#     pvec <- rep((1 - phat) / (S - 1), S)
#     pvec[i] <- phat
#     p_mult_vec[i] <- 
#       exp(
#         log(
#           combinat::dmnom(y, prob = pvec)
#         ) +
#           log_sample_pi[i])
#   }
#   
#   p_outcome_r <- sum(unlist(p_mult_vec))
#   return(p_outcome_r)
# }

infer_sample_of_origin_outcome <- function(..., log_zi, S, phat) {
  y <- c(...)[1:S]
  sample_pi <- c(...)[(S + 1):(2 * S)]
  log_sample_pi <- log(sample_pi)
  # log_zi <- c(...)[2*S+1]
  # zi_vec <- rep(zi, S)
  # posterior_s <- zi_vec ^ y
  posterior_s <- y * log_zi + log_sample_pi

  ## normalize to the sample with maximum posterior probability
  posterior_s <- posterior_s - max(posterior_s)
  posterior_s <- exp(posterior_s)

  posterior <- posterior_s / sum(posterior_s)
  q <- max(posterior)
  s <- which.max(posterior)

 # p_outcome_r <- compute_likelihood(y, S, phat)

  return(list(
  #  p_outcome_r = p_outcome_r,
    s = s,
    q = q
  ))
}


infer_sample_of_origin <- function(outcome_counts, pi_r_hat, S, phat) {
  outcome_counts <-
    left_join(
      outcome_counts,
      pi_r_hat,
      by = c("r"),
      suffix = c(".count", ".pi")
    )

  # zi <- ((S - 1) / (1 / p - 1))
  log_zi <- log((S - 1)) - log(1 / phat - 1)


  posteriors <-
    future_pmap_dfr(
      outcome_counts %>%
        select(matches(".count|.pi$")),
      log_zi = log_zi,
      S = S,
      phat = phat,
      infer_sample_of_origin_outcome
    )

  outcome_counts <-
    bind_cols(outcome_counts, posteriors) %>%
    select(-matches(".pi$|.count$")) 
    #rename_at(vars(matches(".count$")), list(~ str_remove(., ".count$")))

  return(outcome_counts = outcome_counts)
}



reassign_reads_outcome <- function(outcome_counts, S, sample_names) {
  
  
  outcome_counts_reassigned <-
    outcome_counts %>%
    mutate(
      rr = r,
      ss = s
    ) %>%
    spread(ss, rr, fill = 0L) %>%
    rename_at(vars(as.character(1:S)), list(~ paste0(sample_names, "_hat")))  %>%
    select( -k_chimera)%>%
    select(outcome, s, q, n, everything() ) %>%
    arrange(q)

  return(outcome_counts_reassigned)
}




reassign_reads <- function(outcome_counts,
                                 pi_r_hat,
                           sample_names,
                                 phat) {
  
  
  
  
  S <- length(sample_names)
  
  outcome_counts <-
    infer_sample_of_origin(
      outcome_counts,
      pi_r_hat,
      S,
      phat
    )

  outcome_counts <- reassign_reads_outcome(outcome_counts, S, sample_names)

  return(outcome_counts)
}

get_threshold <- function(outcome_counts, summary_stats) {
  
  u <-
    summary_stats$summary_estimates %>%
    pull(u)
  
  g <-
    summary_stats$summary_estimates %>%
    pull(g)
  
  n_cugs <- summary_stats$summary_estimates$n_cugs
  
  current_thresh <-
    outcome_counts %>%
    filter(retain) %>%
    top_n(2, qr) %>%
    mutate(approach = c("torc_before", "discard_torc")) 
  
  
  next_thresh <-
    outcome_counts %>%
    filter(!retain) %>%
    top_n(-1, qr) %>%
    mutate(approach = c("torc_after")) 
  
  no_thresh <-
    outcome_counts %>%
    slice(n()) %>%
    mutate(approach = c("no_discarding")) 
  
  
  cutoff_dt <-
    bind_rows(list(current_thresh, next_thresh, no_thresh)) %>%
    select(approach, outcome, s, qr, n, tor, o,
           FP, FN, TP, TN, FPm, FNm)
  
  no_purging <- list("no_purging", NA, NA, NA, NA, NA , 1 + u, 
               round(get_product_logsum(u + g,  n_cugs)), 0, 
               round(get_product_logsum(1 - g,  n_cugs)), 0, 
               NA, NA)
  
  names(no_purging) <- colnames(cutoff_dt)
  
  cutoff_dt <-
    bind_rows(cutoff_dt, no_purging)
  
  summary_stats$cutoff_dt <- cutoff_dt
  
  return(summary_stats)
}




mark_retained_observations <- function(outcome_counts, summary_stats, torc) {
  
  
  u <-
    summary_stats$summary_estimates %>%
    pull(u)
  
  g <-
    summary_stats$summary_estimates %>%
    pull(g)
  
  n_cugs <- summary_stats$summary_estimates$n_cugs
  
  
  
  outcome_counts <-
    outcome_counts %>%
    group_by(r) %>%
    add_tally(n, name = "nn") %>%
    ungroup() %>%
    mutate(
      p_outcome = n / sum(n)
    ) %>%
    select(-nn)
  
  outcome_counts <-
    outcome_counts %>%
    mutate(
      qr = 1 - q,
      qs = -10 * log10(1e-16) + 10 * log10(qr + 1e-16),
      q=NULL) %>%
    arrange(qr) %>%
    mutate(
      o = cumsum(p_outcome),
      FPp = cumsum(get_product_logsum(
        qr,
        p_outcome
      )),
      FP=round(get_product_logsum(FPp,  n_cugs)),
      FN=round(get_product_logsum(1 - o + FPp - g,  n_cugs)),
      TP = round(get_product_logsum(o - FPp,  n_cugs)),
      TN = round(get_product_logsum(u + g - FPp,  n_cugs)),
      FPm= last(FP)- FP, # marginal
      FNm=FN-last(FN), # marginal
      tor=FNm/FPm)
      
      tor_thresh <- 
        outcome_counts %>%
        summarize(tor_thresh = tor[which.max(-pmax(tor, torc))]) %>%
        pull(tor_thresh)
      
      outcome_counts <-
        outcome_counts %>%
        mutate(
          retain= (tor >= tor_thresh) & FPm >1 ) %>%
        select(outcome, s, qr, r, n, o,
           FP, FN, TP, TN,
           FPm, FNm, tor, retain, qs,
           FPp, p_outcome, everything() ) 
  
  
  
  return(outcome_counts)
}




reassign_reads_and_mark_retained_observations <- function(outcome_counts, summary_stats, sample_names, fit_out, torc){
  
  
  
  outcome_counts <- 
    reassign_reads(
      outcome_counts,
      summary_stats$pi_r_hat,
      sample_names,
      fit_out$glm_estimates$phat
    )
  
  outcome_counts <-
    mark_retained_observations(
      outcome_counts,
      summary_stats,
      torc
    )
  
  

  
  return(outcome_counts)
  
}



purge_and_save_read_counts <- function(read_counts,
                                    dataset_name,
                                    sample_names,
                                    output_dir) {

  read_counts_purged <- 
    read_counts %>%
    select(-sample_names)  %>%
    filter(retain)%>%
    select(-retain) %>%
    rename_at(c(paste0(sample_names, "_hat")), list(~ str_remove(., "_hat$")))
  
  purged_data_filepath <-
    file.path(
      output_dir,
      sprintf(
        "%s_read_counts_purged.rds",
        dataset_name
      )
    )
  
  saveRDS(read_counts_purged, purged_data_filepath)
}
