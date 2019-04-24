rename_var_data_list<- function(read_counts){
  
  sample_names <- names(read_counts)
  
  name_list <- function(x,y) {
    names(x[[1]])=y 
    return(x[[1]])
  }
  
  reads_vars <- 
    rep("reads", length(sample_names)) %>% 
    set_names(sample_names)%>% 
    map(as.list) %>%
    map2(., sample_names, name_list)
  
  rename_var <- function(x, y){
    x <- 
      x %>% 
      rename(!!!y)
    return(x)
  }
  
  read_counts <-
    read_counts %>%
    map2(., reads_vars, rename_var)
  
}


join_read_counts <- function(read_counts) {
  # Join transciptome mapped data by cell, umi, and gene
  read_counts <-
    read_counts %>%
    rename_var_data_list() %>%
    map(setDT) %>%
    reduce(merge,
           all = TRUE,  
           sort= FALSE,
           no.dups= TRUE,
           by = c("cell", "gene", "umi")) %>%
    replace(is.na(.), 0) 

  return(read_counts)
}


add_outcome_variable <- function(read_counts) {
  

  S <- ncol(read_counts) - 3
  sample_names <- colnames(read_counts)[4:(S + 3)]
  
  # read_counts <-
  #   read_counts %>%
  #   mutate(outcome = pmap(select(., sample_names), 
  #                                str_c,
  #                                sep=",",
  #                                collapse = "")) 
  
  setDT(read_counts)

  read_counts[,
       outcome := do.call(paste, c(.SD, sep = ",")),
       .SDcols = sample_names]

  read_counts[order(outcome)]
  
  return(read_counts)
  
}


create_grouping_vars <- function(outcome_counts) {
  
  S <- ncol(outcome_counts) -2
  
  outcome_counts <-
    outcome_counts  %>%
    select(3:(2 + S)) %>%
    as.matrix

  r <- as.integer(rowSums2(outcome_counts))
  k_chimera <- as.integer(S - rowCounts(outcome_counts, value = 0))
  #n_hopped <- as.integer(r - rowMaxs(outcome_counts))
  
  # s_orig <-
  #   apply(cbind(data,r),
  #         1,
  #         function(x) {
  #           y <- which(x[1:S] >=  (x[S + 1] * min_frac))
  #           y <- unname(y)
  #           if (length(y) == 0)
  #             y <- 0
  #           return(y)
  #         })
  
  grouping_vars <- tibble(r = r,
                          #n_hopped = n_hopped,
                          k_chimera = k_chimera)
                          
  return(grouping_vars)
}


create_joined_read_counts_dt <- function(input_dir, read_counts_filepath=NULL) {
  

  if (file.exists(read_counts_filepath)) {
    
    tic("Step 1: reading read counts from existing file")
    read_counts <- readRDS(read_counts_filepath)
    toc()
  }
  
  else{
    
    tic("Step 1: loading molecule_info files, creating read counts datatable, and saving to file")
    read_counts <- load_molecule_info_data(input_dir)
    read_counts <- join_read_counts(read_counts)
    read_counts <- add_outcome_variable(read_counts)
    saveRDS(read_counts, read_counts_filepath)
    toc()
  }

  return(read_counts)
}

add_vars_to_outcome_counts <- function(outcome_counts) {
  
  
  grouping_vars <-
    create_grouping_vars(outcome_counts) 
  outcome_counts <-
    bind_cols(outcome_counts, grouping_vars)  %>%
    arrange(r, k_chimera)%>%
    select(outcome, n, r, k_chimera, everything())
  
  return(outcome_counts)
  
}

create_outcome_counts <- function(read_counts) {
  
  outcome_counts <-
    read_counts %>%
    group_by(outcome) %>%
    add_tally() %>%
    slice(1)  %>%
    select(-c(1:3))  %>%
    select(outcome, n, everything()) %>%
    ungroup()

  outcome_counts <- add_vars_to_outcome_counts(outcome_counts)


  return(outcome_counts)
  
}


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
             fill = list(nn = 0)) %>%
    spread(k_chimera, nn) %>%
    filter(r != 100000) %>%
    ungroup()
   
  # adding tally and proportions variables
  chimera_counts <- 
    chimera_counts %>%
    mutate(
      n_obs = rowSums(.[2:(S + 1)]),
      n_chimeras = rowSums(.[3:(S + 1)]),
      p_chimeras  = n_chimeras / n_obs,
      pcum_obs = cumsum(n_obs) / sum(n_obs)
      #  p_phantoms  = n_phantoms / n_molecules,
    # n_molecules = rowSums(map2_dfr(.[2:(S + 1)], c(1:S), `*`)),
    #  n_phantoms = n_molecules - n_obs,
    #  n_reads = r * n_obs,
     # pcum_molecules = cumsum(n_molecules) / sum(n_molecules),
    #  pcum_reads = cumsum(n_reads) / sum(n_reads)
    )
  
  return(chimera_counts)
}




fit_glm <- function(chimera_counts, S,  r_upper) {
  
  
  fit_dt <-
    chimera_counts %>%
    filter (r <= r_upper & r > 1) %>%
    nest() %>%
    mutate(
      fit = map(data,  ~ glm(
        cbind(n_obs - n_chimeras, n_chimeras) ~ -1 + r,
        data = .x,
        family = binomial(link = log)
      )),
      tidied = map(fit, tidy),
      confint_tidied = map(fit, confint_tidy, conf.level = 0.99),
      #  glanced = map(fit, glance),
      max_r = map(data, ~ as.integer(max(.x$r)))
    )
  
  glm_estimates <-
    fit_dt %>%
    unnest(tidied,
           confint_tidied,
           max_r,
           .drop = TRUE) %>%
    select(-std.error,
           -statistic,
           -p.value) %>%
    mutate_if(is.double,
              .funs = list(p = ~ exp(.)))  %>%
    select(max_r,
           estimate_p,
           conf.low_p,
           conf.high_p)  %>%
    rename(phat = estimate_p,
           phat_low = conf.low_p,
           phat_high = conf.high_p) %>%
    mutate(SIHR = 1 - phat,
           SBIHR = (1 - phat) * ((4 * S - 1) / (4 * S - 4)))
  
  return(glm_estimates)
}
  
update_chimera_counts <- function(chimera_counts, glm_estimates){
  
  chimera_counts <-
    chimera_counts %>%
    mutate(
      phat_chimeras = 1 - glm_estimates$phat ^ r,
      phat_chimeras_low = 1 - glm_estimates$phat_low ^ r,
      phat_chimeras_high = 1 - glm_estimates$phat_high ^ r
      #  nhat_nonchimeras = phat_nonchimeras * n_obs,
      #  nhat_chimeras = phat_chimeras * n_obs
    )
  return(chimera_counts)
}  

estimate_hopping_rate <- function(outcome_counts, r_upper) {
  
  
  # S <- chimera_counts %>%
  #   names()  %>%
  #   str_detect("^[1-9]") %>%
  #   sum()
  S <- ncol(outcome_counts) -4
  
  chimera_counts <- create_chimera_counts(outcome_counts, S)
  
  glm_estimates <- fit_glm(chimera_counts, S=S, r_upper)  

  chimera_counts <- update_chimera_counts(chimera_counts, glm_estimates)
  
  
  return(list(glm_estimates=glm_estimates, chimera_counts=chimera_counts ))
}


estimate_pi_r <- function(nu, phat, S){
  
  pi_r <-(nu* (S - 1) + (phat - 1)) /(S *phat - 1)
  
  pi_r <- pmax(pi_r, 0.000001)
  
  return(pi_r)
  
}


get_reads_dist_estimates <- function(outcome_counts, phat) {
    
    
    S <- ncol(outcome_counts) -4
    
    
    summary_stats <-
      outcome_counts %>%
      mutate(n_obs = 1,
             n_chimeras = (k_chimera != 1)) %>%
      summarize_at(vars(c(3:(S + 5)),
                        n_chimeras, n_obs),
                   list( ~ sum(n * .))) %>%
      rename(n_reads = "r",
             n_molecules = "k_chimera") %>%
      mutate(
        p_chimeras = n_chimeras / n_obs,
        u = (n_molecules - n_obs)/ n_obs,
        n_chimeras= NULL,
        n_molecules=NULL
      ) %>%
      select(n_obs, p_chimeras, u, n_reads, everything()) 
    
    
   marginal <- 
     summary_stats  %>%
     select(5:(S + 4)) %>%
     gather(sample, n_reads) %>%
     mutate(prop_reads= n_reads/sum(n_reads))

   summary_stats <- 
     summary_stats %>%
     select(n_obs, p_chimeras, u, n_reads) 
    
    conditional <-    
      outcome_counts %>%
      ungroup()  %>%
      mutate(n_reads = r,
             n_obs = 1) %>%
      select( r, n, n_reads, n_obs,  everything(), -outcome, -k_chimera)%>%
      group_by(r) %>%
      summarize_at(vars(2:(S + 3)),
                   list( ~ sum(n * .))) %>%
      mutate_at(vars(4:(S + 3)),
                list( ~ . / n_reads)) %>%
      ungroup() %>%
      mutate(m_bar = n_obs / sum(n_obs),
             max_hop_rate= pmap_dbl( select(., 4:(S + 3)), min) * (S - 1)) %>%
      select( r, n_obs, m_bar, everything(), -n_reads) 
    
    return(
      list(
        summary_stats=summary_stats,
        marginal = marginal,
        conditional = conditional
      )
    )
  }

estimate_pi_r_hat_matrix <- function(conditional, S, phat){
  
  pi_r_hat <-
    conditional %>%
    mutate_at(vars(4:(S + 3)),  estimate_pi_r, S=S, phat=phat) %>%
    mutate(sum_p= rowSums(.[4:(S + 3)])) %>%
    mutate_at(vars(4:(S + 3)), ~ ./sum_p) %>%
    select(-sum_p, -max_hop_rate,- n_obs,-  m_bar)
  
  return(pi_r_hat)
}


get_product_logsum <- function(x, y) {
  z_log <- log(x) + log(y)
  z <- exp(z_log)
  return(z)
}

infer_sample_of_origin_outcome <- function(..., log_zi, S, phat) {
  y <- c(...)[1:S]
  sample_pi <- c(...)[(S + 1):(2 * S)]
  log_sample_pi <- log(sample_pi)
  #log_zi <- c(...)[2*S+1]
  #zi_vec <- rep(zi, S)
  #posterior_s <- zi_vec ^ y
  posterior_s <- y * log_zi + log_sample_pi
  
  ## normalize to the sample with maximum posterior probability
  posterior_s <- posterior_s - max(posterior_s)
  posterior_s <- exp(posterior_s)

  posterior <- posterior_s / sum(posterior_s)
  q <- max(posterior)
  s <- which.max(posterior)
  
  p_mult_vec <- rep(0,S)
  for(i in 1:S) {
    pvec <- rep((1 - phat) / (S - 1), S)
    pvec[i] <- phat
    p_mult_vec[i] <- exp(log(combinat::dmnom(y,
                                             prob=pvec)) +
                           log_sample_pi[i])
  }

  p_outcome_r <- sum(unlist(p_mult_vec))
  
  return(list(p_outcome_r=p_outcome_r,
              s = s,
              q = q))
}


compute_classification_metrics <- function(outcome_counts, g, u) {
  
  
  outcome_counts <- 
    outcome_counts %>% 
    group_by(r) %>% 
    add_count( wt= p_outcome_r, 
               name="p_nonmiss") %>%
    add_tally(n, name = "nn") %>%
    ungroup() %>%
    mutate(
      prop_r = nn / sum(n),
     # p_r= prop_r/p_nonmiss,
      prop_outcome_r = n / nn,
      #p_outcome=get_product_logsum(p_outcome_r, prop_r),
      prop_outcome = n / sum(n)
    ) %>%
    select(-matches(".ct|.pi$"))%>%
    select(-nn)
  
  outcome_counts <-
    outcome_counts %>%
    mutate(qr=1-q)%>%
    arrange(qr) %>%
    mutate(
      #marg_p_outcome =get_product_logsum(p_outcome,nn_prop),
      # cum_p=cumsum(marg_p_outcome),
      o = cumsum(prop_outcome),
      FP = cumsum(get_product_logsum(qr,
                                     prop_outcome)),
      FN =1- o + FP - g,
      # pcum_phantom = cumsum(get_product_logsum((k_chimera - 1) * q +
      #                                            (k_chimera) * (1 - q),
      #                                          p_outcome_m
      # )),
      #pcum_reads = cumsum(r * n) / sum(r * n),
      #pcum_phantom = last(pcum_phantom) - pcum_phantom,
      qs =  -10*log10( 1e-16) + 10 * log10(qr + 1e-16),
      FPR = FP / (u + g),
      FNR = FN / (1 - g),
      j = 1 - (FPR + FNR)
    ) 
  
  return(outcome_counts)
  
}


infer_sample_of_origin <-
  function(outcome_counts,
           pi_r_hat,
           S,
           phat) {
    

    outcome_counts <-
      left_join(
        outcome_counts,
        pi_r_hat,
        by = c("r"),
        suffix = c(".ct", ".pi")
      )
    
    
    #zi <- ((S - 1) / (1 / p - 1))
    
    log_zi <- log((S - 1)) - log(1 / phat - 1)
    
    
    posteriors <-
      future_pmap_dfr(
        outcome_counts %>%
          select(matches(".ct|.pi$")),
        log_zi = log_zi,
        S = S,
        phat=phat,
        infer_sample_of_origin_outcome
      )
    
    outcome_counts <- bind_cols(outcome_counts, posteriors)
    
    return(outcome_counts=outcome_counts)
    
  }


get_prop_fugues <- function(m_bar, phat, r_max=10) {

  r <- 1:r_max
  m <- m_bar[r]
  y <- rep((1 -phat), r_max)
  g <- sum(exp(r*log(y)+log(m)))
  
  return(g)
}

compute_prop_nonmissingness <- function(outcome_counts){
  
  p_nonmiss_dt <-
    outcome_counts %>%
    group_by(r)%>%
    slice(1)%>%
    select(r, p_nonmiss)
  
  return(p_nonmiss_dt)
  
}

update_reads_dist_summary <- function(reads_dist_summary,p_nonmiss_dt, g){
  
  reads_dist_summary$summary_stats <- 
    reads_dist_summary$summary_stats %>%
    mutate(g=g) %>%
    select(n_obs, p_chimeras, g, u, everything())
  
  reads_dist_summary$conditional <- left_join(reads_dist_summary$conditional,
                                              p_nonmiss_dt, 
                                              by="r")
  
  return(reads_dist_summary)
}

get_q_cutoff <- function(optimal_cutoff,
                          outcome_counts,
                          max_fpr = max_fpr) {
  
  if (is.null(max_fpr)) {
    
    score_threshold <-
      optimal_cutoff %>%
      filter(cutoff == "optimal") %>%
      pull(q)
  } else{
    score_threshold <-
      outcome_counts %>%
      summarize(argmax_j = q[which.max(-pmax(max_fpr - FPR, 0))]) %>%
      pull(argmax_j)
  }
  return(score_threshold)
}

get_threshold <- function(outcome_counts, argmax_j) {
  current_thresh <- 
    outcome_counts%>%
    filter(q >= argmax_j) %>%
    top_n(-2, q) %>%
    arrange(desc(j)) %>%
    mutate(cutoff = c("optimal", "above")) %>%
    select(cutoff, outcome, q, s, FPR, j, qs, o, FP)
  
  next_thresh <- 
    outcome_counts %>%
    filter(q < argmax_j) %>%
    top_n(1, q)  %>%
    mutate(cutoff = c("below")) %>%
    select(cutoff, outcome, q, s, FPR, j, qs, o, FP)
  
  no_thresh <- 
    outcome_counts %>%
    slice(n()) %>%
    mutate(cutoff = c("none"))%>%
    select(cutoff, outcome, q, s, FPR, j, qs, o, FP)
  
  threshold_dt <-
    bind_rows(list(current_thresh, next_thresh, no_thresh)) 
  
  return(threshold_dt)
}


compute_optimal_cutoff <- function(outcome_counts, u, g) {
  argmax_j  <-
    outcome_counts%>%
    summarize(argmax_j = q[which.max(j)]) %>%
    pull(argmax_j)

  threshold_dt <- get_threshold(outcome_counts, argmax_j)
  
  threshold_dt <-
    threshold_dt %>% 
      mutate( TP= o-FP, 
              FN=1-g-TP,
              TN=u+g -FP)

  return(threshold_dt)
  
}


reassign_hopped_reads <- function(read_counts,
                                  outcome_counts) {
  
  S <- ncol(read_counts) -4

  
  tic("Step 10.1: reassigning reads to sample of origin")
  outcome_counts <-
    outcome_counts %>%
    separate(
      outcome,
      into = paste0("m", 1:S),
      sep = ",",
      remove = FALSE,
      convert = TRUE
    ) %>%
    spread (s, r, fill = 0L) %>%
    rename_at(vars(as.character(1:S)), list( ~ paste0("t", .))) %>%
    select(c("outcome", "n", "q", paste0("m", 1:S), paste0("t", 1:S)))
  toc()
  
  
  tic("Step 10.2: deduplicating read counts")
  outcome_counts <-
    outcome_counts %>%
    mutate_at(c(paste0("m", 1:S), paste0("t", 1:S)),
              list( ~ as.integer(. > 0)))
  toc()
  
  return(outcome_counts)
  
}


purge_phantoms <- function(read_counts,
                           outcome_counts,
                           score_threshold) {
  
  S <- ncol(read_counts) -4
  sample_names <- colnames(read_counts)[4:(S+3)]
  
  tic("Step 10.3: reassigning hopped reads and deduplicating read counts")
  outcome_counts_dedupped <-
    reassign_hopped_reads(read_counts,
                          outcome_counts)
  
  toc()
  
  
  tic("Step 10.4: labelling phantom molecules below cutoff")
  
  outcome_counts_dedupped <-
    outcome_counts_dedupped %>%
    mutate(keep = q >= score_threshold)
  toc()
  
  
  tic("Step 10.5: adding cell-umi-gene labels")
  read_counts <-
    join_data(read_counts %>%
                select(outcome, cell, umi, gene),
              outcome_counts_dedupped) %>%
    select( -c("n", "q"))
  toc()
  
  
  tic("Step 10.6: tallying molecule counts by cell-barcode and gene ID")
  setDT(read_counts)
  umi_counts_cell_gene <-
    read_counts[, lapply(.SD, sum),
                keyby = 'cell,gene,keep',
                .SDcols = -c("umi")]
  toc()

  
  tic("Step 10.7: tranforming cell-gene molecule tally table into long format")
  
  umi_counts_cell_gene <-
    melt(
      umi_counts_cell_gene,
      id = 1:3,
      variable.name = "sample",
      measure = patterns(m = "^m", t = "^t"),
      variable.factor = FALSE
    )
  
  sample_key <-
    tibble(sample = as.character(1:S),
           sample_name = sample_names)
  
  umi_counts_cell_gene <-
    left_join(umi_counts_cell_gene,
              sample_key,
              by = "sample") 
  
  umi_counts_cell_gene <- 
    umi_counts_cell_gene %>%
    mutate(sample = sample_name,
           f=m-t,
           r_b=t*(1-keep),
           r_a=t-r_b,
           t=NULL,
           sample_name = NULL,
           keep=NULL) 
  toc()
  
  
  umi_counts_cell_gene <- split(umi_counts_cell_gene,
                                umi_counts_cell_gene$sample)
  
  
  return(umi_counts_cell_gene)
}


get_umi_counts_cell <- function( sample_called_cells, sample_umi_counts_cell_gene) {
  

  sample_umi_counts_cell_gene <- left_join(sample_called_cells %>% select(-is_cell_purged),
                                    sample_umi_counts_cell_gene,
                  by= "cell")
  
  sample_umi_counts_cell <-
    sample_umi_counts_cell_gene  %>%
    select(-c("gene"))  %>%
    group_by(is_cell_unpurged, cell) %>%
    summarize(m = sum(m),
              r_a = sum(r_a),
              f = sum(f),
              r_b = sum(r_b)) %>%
    arrange( -f) %>%
    ungroup()
  
  sample_umi_counts_cell <- 
    split(sample_umi_counts_cell %>%
            select(-is_cell_unpurged), 
          sample_umi_counts_cell$is_cell_unpurged) %>%
    set_names(c("background_cells", "called_cells"))
    
  
  return(sample_umi_counts_cell)
  
}



get_umi_counts_sample <- function(umi_counts_cell){
  
  umi_counts_sample <-
    umi_counts_cell  %>%
    ungroup() %>%
    select(c( "m", "r_a", "r_b", "f"))  %>%
    summarize(m = sum(m),
              r_a= sum(r_a),
              r_b = sum(r_b),
              f = sum(f)) 
  
}

update_reads_dist_summary_2 <- function(reads_dist_summary, umi_counts_sample){
  

  umi_counts_sample_all <- 
    umi_counts_sample %>%  
    select(-split) %>% 
    group_by(sample)%>% 
    summarize_all(list( . ~sum))%>% 
    mutate(r_a= r_a/m,
           r_b=r_b/m,
           f=f/m)%>% 
    ungroup() %>%
    mutate(prop_m=m/sum(m))
  
  umi_counts_sample_called_cells <- 
    umi_counts_sample %>%  
    filter(split=="called_cells") %>%  
    select(-split) %>% 
    group_by(sample)%>% 
    summarize_all(list( . ~sum))%>% 
    mutate(r_a= r_a/m,
           r_b=r_b/m,
           f=f/m)%>% 
    ungroup() %>%
    mutate(prop_m=m/sum(m))
  
  reads_dist_summary$marginal_called_cells <- 
    left_join(umi_counts_sample_called_cells,
              reads_dist_summary$marginal, 
              by="sample") %>%
    mutate(FRM= n_reads/m,
           m= m/1000000) %>%
    select(sample, m, r_a, r_b, f, prop_m, prop_reads, FRM)
  
  reads_dist_summary$marginal <- 
    left_join(umi_counts_sample_all,
              reads_dist_summary$marginal, 
              by="sample") %>%
    mutate(FRM= n_reads/m,
           m= m/1000000) %>%
    select(sample, m, r_a, r_b, f, prop_m, prop_reads, FRM)
  
  return(reads_dist_summary)
}


make_umi_count_matrices <- function(umi_counts_cell_gene){
  
  
  genes <- unique(umi_counts_cell_gene$gene)
  
  sample_data <- 
    umi_counts_cell_gene %>% 
    filter( m>0)
  
  
  unpurged <- DropletUtils::makeCountMatrix(sample_data$gene, 
                                                    sample_data$cell,
                                                    all.genes=genes,
                                                    value=sample_data$m)
  
  sample_data <- 
    sample_data %>% 
    filter(r_a > 0)
  
  
  purged <- DropletUtils::makeCountMatrix(sample_data$gene, 
                                                     sample_data$cell,
                                                     all.genes=genes,
                                                     value=sample_data$r_a)
  
  return(list(unpurged=unpurged,
              purged=purged))
  
  
}


call_emptydrops <- function(sample_umi_count, lower=200, fdr_thresh=0.005){
  
  sample_emptydrops <- 
    DropletUtils::emptyDrops(sample_umi_count, 
                             lower=lower)  %>%
    as.data.frame() %>% 
    rownames_to_column(var = "cell") %>% 
    select(cell, FDR)
  
  sample_emptydrops[is.na(sample_emptydrops)] <- 1
  
  sample_emptydrops <- 
    sample_emptydrops%>%
    mutate(is_cell= FDR <= fdr_thresh,
           FDR=NULL) 
  
  return(sample_emptydrops)
}

call_defaultdrops <- function(sample_umi_count){
  
  sample_defaultdrops  <- DropletUtils::defaultDrops(sample_umi_count) 
  sample_defaultdrops  <- enframe(sample_defaultdrops, 
                                  name="cell", 
                                  value="is_cell")
  
  return(sample_defaultdrops)
  
}

call_cells <- function(umi_count_matrix, sample_name){

  called_cells <- future_map(umi_count_matrix, 
                      plyr::failwith(NULL, 
                                     call_emptydrops))
    

  if(is.null(called_cells$unpurged) | is.null(called_cells$purged)){
    message(paste0("Using defaultDrops cell calling instead.
                   EmptyDrops function call failed for purged sample ",
                   sample_name))
    called_cells <- map(umi_count_matrix, call_defaultdrops)
    
  } 
  

  called_cells <- left_join(called_cells$unpurged,
                            called_cells$purged,
                            by="cell",
                            suffix=c("_unpurged", "_purged"))
  
  
  return(called_cells)
  
}


get_cells_tally <- function(called_cells, sample_name){

  #keys <- c(`FALSE`="background", `TRUE`="cell")

  cells_tally <- 
    called_cells %>%
    group_by(is_cell_unpurged,
             is_cell_purged)%>%
    tally() %>%
    ungroup()%>%
    complete(is_cell_unpurged,is_cell_purged, fill =list(n = 0)) %>%
    mutate(
           #is_cell_purged= recode(as.character(is_cell_purged), !!!keys),
          # is_cell_unpurged= recode(as.character(is_cell_unpurged), !!!keys),
           barcode= c("consensus_background",
                       "transition_cell",
                       "phantom_background",
                       "transition_background",
                       "consensus_cell" ,
                       "phantom_cell")) %>%
    select(barcode, n) %>%
    group_by(barcode)%>%
    set_names("barcode", sample_name)
    
    

return(cells_tally)

}

get_purged_umi_counts <- function(called_cells_sample, 
                                  umi_count_matrix_sample){
  
  umi_purged <- umi_count_matrix_sample$purged
  keep <- drop_na(called_cells_sample)$is_cell_purged
  umi_purged <-  umi_purged[ , keep]

  return(umi_purged)
}


join_data <- function(data, outcome_counts) {
  data <-
    left_join(data,
              outcome_counts ,
              by = c("outcome"))  %>%
    select(-outcome)
  
  return(data)
}