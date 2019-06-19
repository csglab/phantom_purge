create_chimera_counts <- function(outcome_counts, S) {

  # creating a data table of k_chimera counts
  chimera_counts <-
    outcome_counts %>%
    group_by(r, k_chimera) %>%
    tally(n, name = "nn") %>%
    bind_rows(tibble(
      r = rep(100000, S),
      k_chimera = rep(1:S),
      nn = 0
    )) %>%
    ungroup()

  # Converting table from long to wide format
  chimera_counts <-
    chimera_counts %>%
    complete(r,
      k_chimera,
      fill = list(nn = 0)
    ) %>%
    spread(k_chimera, nn) %>%
    filter(r != 100000) %>%
    ungroup()

  # adding tally and proportions variables
  chimera_counts <-
    chimera_counts %>%
    mutate(
      n_cugs = rowSums(.[2:(S + 1)]),
      n_chimeras = rowSums(.[3:(S + 1)]),
      p_chimeras = n_chimeras / n_cugs,
      pcum_cugs = cumsum(n_cugs) / sum(n_cugs)
    )

  return(chimera_counts)
}


fit_glm <- function(chimera_counts, S, max_r, conf_level = 0.99) {
  
  fit_dt <-
    chimera_counts %>%
    filter(r <= max_r & r > 1) %>%
    nest() %>%
    mutate(
      fit = map(data, ~ glm(
        cbind(n_cugs - n_chimeras, n_chimeras) ~ -1 + r,
        data = .x,
        family = binomial(link = log)
      )),
      tidied = map(fit, tidy),
      confint_tidied = map(fit, confint_tidy, conf.level = conf_level),
      #  glanced = map(fit, glance),
      max_r = map(data, ~ as.integer(max(.x$r)))
    )

  glm_estimates <-
    fit_dt %>%
    unnest(tidied,
      confint_tidied,
      max_r,
      .drop = TRUE
    ) %>%
    select(
      -std.error,
      -statistic,
      -p.value
    ) %>%
    mutate_if(is.double,
      .funs = list(p = ~ exp(.))
    ) %>%
    select(
      max_r,
      estimate_p,
      conf.low_p,
      conf.high_p
    ) %>%
    rename(
      phat = estimate_p,
      phat_low = conf.low_p,
      phat_high = conf.high_p
    ) %>%
    mutate(
      SIHR = 1 - phat,
      SBIHR = (1 - phat) * ((4 * S - 1) / (4 * S - 4))
    )

  return(glm_estimates)
}

update_chimera_counts <- function(chimera_counts, glm_estimates) {
  
  chimera_counts <-
    chimera_counts %>%
    mutate(
      phat_chimeras = 1 - glm_estimates$phat^r,
      phat_chimeras_low = 1 - glm_estimates$phat_low^r,
      phat_chimeras_high = 1 - glm_estimates$phat_high^r
    )
  return(chimera_counts)
}

estimate_hopping_rate <- function(outcome_counts, S, max_r = NULL) {
  
  chimera_counts <- create_chimera_counts(outcome_counts, S)

  if (is.null(max_r)) {
    max_r <-
      chimera_counts %>%
      pull(r) %>%
      max()
  }

  glm_estimates <-
    fit_glm(chimera_counts,
      S = S,
      max_r
    )

  chimera_counts <-
    update_chimera_counts(
      chimera_counts,
      glm_estimates
    )

  return(list(
    glm_estimates = glm_estimates,
    chimera_counts = chimera_counts
  ))
}