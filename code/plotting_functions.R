
plot_molecules_distributions <- function(data_list, dataset_name, x_lim=150) {
  
  pi_r_hat <- 
    bind_cols(data_list$pi_r_hat,
    data_list$reads_dist_summary$conditional %>%
    select(n_obs, m_bar, p_nonmiss))
  
  marginal_prop <- data_list$reads_dist_summary$marginal 
  
  
  p1 <- 
    ggplot(pi_r_hat %>%
                 select(-n_obs,
                        -m_bar,
                        - p_nonmiss) %>%
                 gather(sample, p, -c("r"))) +
    geom_line(aes(x = r,
                  y = p, 
                  colour = sample),
              size = 1,
              alpha = 1) +
    geom_point(
      data = pi_r_hat %>%
        select(r, m_bar),
      aes(x = r, y = m_bar),
      shape = 10,
      size = 1,
      alpha = .9
    ) +
    geom_point(
      data = pi_r_hat %>%
        select(r, p_nonmiss),
      aes(x = r, y = p_nonmiss),
      shape = 10,
      size = .5,
      alpha = .7
    ) +
    theme_bw() +
    geom_hline(
      data = marginal_prop,
      aes(yintercept = prop_m,
          colour = sample),
      size = 0.3,
      linetype = "longdash",
      alpha = 1
    ) +
    labs(x = "r (PCR duplicates)",
         y = "proportion") +
    xlim(0, x_lim) +
    scale_colour_viridis_d(option = "inferno",
                           direction = -1) +
    theme(axis.title.x = element_text(size = rel(1.4)),
          axis.title.y = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.1), face="bold"),
          legend.text=element_text(size=rel(1.1)))
  
  p2 <- 
    ggplot(pi_r_hat) +
    geom_line(
      aes(x = r, y = n_obs),
      size = 1,
      alpha = 1,
      colour = "coral"
    ) +
    # facet_grid(dataset~.) +
    theme_bw() +
    labs(x = "r",
         y = "count") +
    scale_y_log10() +
    theme(axis.title.x = element_text(size = rel(.7)),
          axis.title.y = element_text(size = rel(.7)),
          axis.text.y = element_blank())#+
    #geom_vline(data=sample_stats,
    #           aes(xintercept = r_thresh),
     #          linetype = "dashed") 
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position='none')
  
  
  p_rdist<-
    ggdraw() +
    draw_plot(p1, 0, 0, 1, 1) +
    draw_plot(p2, 0.75, 0.75, .23, .23)   +
    draw_label(dataset_name,
               x = 0,
               y = 1,
               vjust = 4, 
               hjust = -1.2,
               size = 12,
               fontface = 'italic')
  
  return(list(p=p_rdist, 
              legend=legend))
}


plot_fit <- function(data_list,
                     dataset_name,
                     x_lim = 150) {
  
  chimera_counts <- data_list$fit_out$chimera_counts
  glm_estimates <- data_list$fit_out$glm_estimates
  
  p1 <- ggplot(chimera_counts) +
    geom_point(aes(r,
                   p_chimeras,
                   colour = "observed"),
               size = 1) +
    geom_line(aes(r,
                  phat_chimeras,
                  colour = "predicted")) +
    geom_line(aes(r,
                  pcum_obs,
                  colour = "ECDF r")) +
    theme_bw() + 
    labs(x = "r (PCR duplicates)",
                      y = "proportion") +
    geom_vline(
      data = glm_estimates,
      mapping = aes(xintercept = max_r),
      linetype = "dashed"
    ) +
    scale_color_manual(
      name = "",
      values = c("darkgrey",
                 "red",
                 "blue",
                 "coral")
    )+
    xlim(0, x_lim)  +
    theme(axis.title.x = element_text(size = rel(1.4)),
          axis.title.y = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.1), face="bold"),
          legend.text=element_text(size=rel(1.1)))
  
  p2 <- ggplot(chimera_counts) +
    geom_point(aes(r,
                   n_chimeras,
                   colour = "observed"),
               size = 1) +
    geom_line(aes(r,
                  n_obs*phat_chimeras,
                  colour = "predicted"))  +
    theme_bw() +
    labs(x = "r",
         y = "count") +
    scale_color_manual(
      name = "",
      values = c("darkgrey","blue")
    ) +
    annotate(
      "text",
      x = x_lim - 5,
      y = max(chimera_counts$n_chimeras) - 150,
      label = sprintf("SIHR=%.4f",
                      glm_estimates$SIHR),
      size = 3,
      hjust = 1,
      vjust = 1
    ) +
    xlim(0, x_lim) +
    theme(axis.title.x =  element_text(size = rel(.7)),
          axis.title.y =  element_text(size = rel(.7))) +
    theme(legend.position = "none")
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position='none')
  
  
  p_fit <-
    ggdraw() +
    draw_plot(p1, 0, 0, 1, 1) +
    draw_plot(p2, 0.65, 0.18, .3, .3)   +
    draw_label(dataset_name,
               x = 0,
               y = 1,
               vjust = 2, 
               hjust = -1.2,
               size = 12,
               fontface = 'italic')
  
  
  return(list(p=p_fit,legend=legend) )
  
}


plot_posterior_prob <- function(data_list,
                                dataset_name) {
  
  outcome_counts <- data_list$outcome_counts
  
  optimal_cutoff <- 
    data_list$optimal_cutoff %>% 
    filter(cutoff=="optimal")
  
  p1 <- ggplot(outcome_counts) +
    geom_line(aes(x = qs,
                  y = FNR,
                  colour = "FNR"),
              alpha = 0.7)  +
    geom_line(aes(x = qs,
                  y = FPR,
                  colour = "FPR"),
              alpha = 0.7)  +
    geom_point(aes(x = qs,
                   y = o,
                   colour = "ECDF (o)"),
               size = 1,
               shape=18) + 
    geom_line(aes(x = qs,
                  y = j,
                  colour="J"),
              size = 1,
              alpha = 0.7)  +
    geom_vline(data= optimal_cutoff, 
               aes(xintercept = qs),
               linetype = "dashed",
               color = "grey") +
    annotate(
      "text",
      x = optimal_cutoff$qs + 5,
      y = .1,
      label = sprintf("qs*=%.1f",
                      optimal_cutoff$qs),
      size = 4
    )  +
    labs(x = "qs",
         y = "proportion") +
    theme_bw() +
    scale_color_manual(name = "",
                       values = c("red",
                                  "green",
                                  "darkblue",
                                  "purple"
                                  
                       ))   +
    theme(axis.title.x = element_text(size = rel(1.4)),
          axis.title.y = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.3)),
          legend.title = element_text(size = rel(1.1), face="bold"),
          legend.text=element_text(size=rel(1.1)))
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position='none')
  
  p_fit <-
    ggdraw() +
    draw_plot(p1, 0, 0, 1, 1) +
    draw_label(dataset_name,
               x = 0,
               y = 1,
               vjust = 2, 
               hjust = -1.3,
               size = 12,
               fontface = 'italic')
  
  return(list(p=p_fit,legend=legend) )
  
}



make_plots <- function(data_list, dataset_name, x_lim = 160) {
  p_rdist <- plot_read_distributions(data_list$reads_dist_summary)
  
  
  p_fit <- plot_fit(data_list$fit_out$chimera_counts,
                    data_list$fit_out$glm_estimates,
                    x_lim = x_lim)
  
  
  p_class <- plot_posterior_prob(data_list$outcome_counts,
                                 # data_list$fit_out$glm_estimates,
                                 data_list$optimal_cutoff,
                                 x_lim = x_lim,
                                 y_lim = 1)
  
  p_phantoms <- plot_phantoms(data_list$umi_counts_cell,
                              data_list$umi_counts_gene,
                              data_list$umi_counts_sample)
  
  
  plot_list <- list(
    p_rdist = p_rdist,
    p_fit = p_fit,
    p_class = p_class,
    p_phantoms = p_phantoms
  )
  saveRDS(plot_list,
          file.path(output_dir,
                    sprintf("%s_ggplot_objects.rds",
                            dataset_name)))
  
  
  return(plot_list)
  
}