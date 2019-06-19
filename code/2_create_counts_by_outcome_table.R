
create_s_maxprop_var <- function(outcome_counts, r, sample_names, min_frac) {
  
  S <- length(sample_names)
  
  s_maxprop <-
    apply(
      cbind(outcome_counts, r),
      1,
      function(x) {
        y <- which(x[1:S] >= (x[S + 1] * min_frac))
        y <- unname(y)
        if (length(y) == 0) {
          y <- 0
        }
        return(y)
      }
    )

  return(s_maxprop)
}

create_grouping_vars <- function(outcome_counts, sample_names, min_frac) {
  S <- length(sample_names)

  outcome_counts <-
    outcome_counts %>%
    select(sample_names) %>%
    as.matrix()

  r <- as.integer(rowSums2(outcome_counts))
  k_chimera <- as.integer(S - rowCounts(outcome_counts, value = 0))


  grouping_vars <- tibble(
    r = r,
    k_chimera = k_chimera
  )

  #grouping_vars$n_hopped <- as.integer(r - rowMaxs(outcome_counts))
  

  if (!is.null(min_frac)) {
    grouping_vars$s_maxprop <- create_s_maxprop_var(
      outcome_counts,
      r,
      sample_names,
      min_frac
    )
  }
  return(grouping_vars)
}


add_vars_to_outcome_counts <- function(outcome_counts, sample_names, min_frac) {
  grouping_vars <-
    create_grouping_vars(outcome_counts, sample_names, min_frac)
  outcome_counts <-
    bind_cols(outcome_counts, grouping_vars) %>%
    arrange(r, k_chimera) %>%
    select(outcome, n, r, k_chimera, everything(), sample_names)

  return(outcome_counts)
}

create_outcome_counts <- function(read_counts, sample_names, min_frac=NULL) {
  outcome_counts <-
    read_counts %>%
    group_by(outcome) %>%
    add_tally() %>%
    slice(1) %>%
    select(-c(1:3)) %>%
    select(outcome, n, sample_names) %>%
    ungroup()

  outcome_counts <- add_vars_to_outcome_counts(
    outcome_counts,
    sample_names,
    min_frac
  )


  return(outcome_counts)
}