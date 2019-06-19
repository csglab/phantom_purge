
get_prop_fugues <- function(p_cugs, phat) {
  upper_thresh <-10
  r <- 1:upper_thresh # first ten terms are sufficient
  m <- p_cugs[r]
  y <- rep((1 - phat), upper_thresh)
  g <- sum(exp(r * log(y) + log(m)))

  return(g)
}


estimate_pi_r <- function(nu, phat, S, pi_r_eps = 0.000001) {
  pi_r <- (nu * (S - 1) + (phat - 1)) / (S * phat - 1)
  pi_r <- pmax(pi_r, pi_r_eps)

  return(pi_r)
}

estimate_pi_r_hat_matrix <- function(conditional, sample_names, phat) {
  S <- length(sample_names)

  pi_r_hat <-
    conditional %>%
    mutate_at(vars(sample_names),
      estimate_pi_r,
      phat = phat,
      S = S
    ) %>%
    mutate(sum_p = rowSums(.[sample_names])) %>%
    mutate_at(vars(sample_names), ~ . / sum_p) %>%
    select(-sum_p, -max_hop_rate, -n_cugs, -p_cugs)

  return(pi_r_hat)
}


compute_summary_stats <- function(outcome_counts, phat, sample_names) {
  S <- length(sample_names)

  summary_estimates <-
    outcome_counts %>%
    mutate(
      n_cugs = 1,
      n_chimeras = (k_chimera != 1)
    ) %>%
    summarize_at(
      vars(
        r,
        k_chimera,
        sample_names,
        n_chimeras,
        n_cugs
      ),
      list(~ sum(n * .))
    ) %>%
    rename(
      n_reads = "r",
      n_molecules = "k_chimera"
    ) %>%
    mutate(
      p_chimeras = n_chimeras / n_cugs,
      u = (n_molecules - n_cugs) / n_cugs,
      n_chimeras = NULL,
      n_molecules = NULL
    ) %>%
    select(n_cugs, p_chimeras, u, n_reads, everything())


  marginal <-
    summary_estimates %>%
    select(sample_names) %>%
    gather(sample, n_reads) %>%
    mutate(prop_reads = n_reads / sum(n_reads))
  
  
  marginal <-
    outcome_counts %>%
    summarize_at(sample_names,
                 list(~ sum(n * as.integer(. >0)))
    ) %>% 
    data.table::transpose() %>%
    set_names("n_molecs")  %>%
    bind_cols(marginal,.) %>%
    mutate(p_molecs= n_molecs/sum(n_molecs),
           FRM = n_reads / n_molecs#,
           #  n_molecs = n_molecs / 1000000,
           #   n_reads = n_reads / 1000000
    )

  summary_estimates <-
    summary_estimates %>%
    select(n_cugs, p_chimeras, u, n_reads)

  conditional <-
    outcome_counts %>%
    ungroup() %>%
    mutate(
      n_reads = r,
      n_cugs = 1
    ) %>%
    select(r, n, n_reads, n_cugs, everything(), -outcome, -k_chimera) %>%
    group_by(r) %>%
    summarize_at(
      vars(n_reads, n_cugs, sample_names),
      list(~ sum(n * .))
    ) %>%
    mutate_at(
      vars(sample_names),
      list(~ . / n_reads)
    ) %>%
    ungroup() %>%
    mutate(
      p_cugs = n_cugs / sum(n_cugs),
      max_hop_rate = pmap_dbl(select(., sample_names), min) * (S - 1)
    ) %>%
    select(r, n_cugs, p_cugs, everything(), -n_reads)

  g <- get_prop_fugues(conditional$p_cugs, phat)


  summary_estimates <-
    summary_estimates %>%
    mutate(g = g,
           n_pm=round(n_cugs*(u+g)),
           frm=n_reads/n_cugs) %>%
    select(n_reads, n_cugs, n_pm, u, g, frm, p_chimeras, everything())

  pi_r_hat <-
    estimate_pi_r_hat_matrix(
      conditional,
      sample_names,
      phat
    )

  return(
    list(
      summary_estimates = summary_estimates,
      marginal = marginal,
      conditional = conditional,
      pi_r_hat = pi_r_hat
    )
  )
}