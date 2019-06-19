

plot_molecules_distributions <- function(data_list, dataset_name, x_lim = 150) {
  pi_r_hat <-
    bind_cols(
      data_list$summary_stats$pi_r_hat,
      data_list$summary_stats$conditional %>%
        select(n_cugs, p_cugs)
    )

  marginal_prop <- data_list$summary_stats$marginal


  p1 <-
    ggplot(pi_r_hat %>%
      select(
        -n_cugs,
        -p_cugs
      ) %>%
      gather(sample, p, -c("r"))) +
    geom_line(aes(
      x = r,
      y = p,
      colour = sample
    ),
    size = 1,
    alpha = 1
    ) +
    geom_point(
      data = pi_r_hat %>%
        select(r, p_cugs),
      aes(x = r, y = p_cugs),
      shape = 10,
      size = 1,
      alpha = .9
    ) +
    theme_bw() +
    geom_hline(
      data = marginal_prop,
      aes(
        yintercept = p_molecs,
        colour = sample
      ),
      size = 0.3,
      linetype = "longdash",
      alpha = 1
    ) +
    labs(
      x = "r (PCR duplicates)",
      y = "proportion"
    ) +
    xlim(0, x_lim) +
    scale_colour_viridis_d(
      option = "inferno",
      direction = -1
    ) +
    theme(
      axis.title.x = element_text(size = rel(1.4)),
      axis.title.y = element_text(size = rel(1.4)),
      axis.text = element_text(size = rel(1.3)),
      legend.title = element_text(size = rel(1.1), face = "bold"),
      legend.text = element_text(size = rel(1.1))
    )

  p2 <-
    ggplot(pi_r_hat) +
    geom_line(
      aes(x = r, y = n_cugs),
      size = 1,
      alpha = 1,
      colour = "coral"
    ) +
    # facet_grid(dataset~.) +
    theme_bw() +
    labs(
      x = "r",
      y = "count"
    ) +
    scale_y_log10() +
    theme(
      axis.title.x = element_text(size = rel(.7)),
      axis.title.y = element_text(size = rel(.7)),
      axis.text.y = element_blank()
    ) #+
  # geom_vline(data=sample_stats,
  #           aes(xintercept = r_thresh),
  #          linetype = "dashed")

  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")


  p_rdist <-
    ggdraw() +
    draw_plot(p1, 0, 0, 1, 1) +
    draw_plot(p2, 0.75, 0.75, .23, .23) +
    draw_label(dataset_name,
      x = 0,
      y = 1,
      vjust = 4,
      hjust = -1.2,
      size = 12,
      fontface = "italic"
    )

  return(list(
    p = p_rdist,
    legend = legend
  ))
}


plot_fit <- function(data_list,
                     dataset_name,
                     x_lim = 150) {
  chimera_counts <- data_list$fit_out$chimera_counts
  glm_estimates <- data_list$fit_out$glm_estimates

  p1 <- ggplot(chimera_counts) +
    geom_point(aes(r,
      p_chimeras,
      colour = "observed"
    ),
    size = 1
    ) +
    geom_line(aes(r,
      phat_chimeras,
      colour = "predicted"
    )) +
    geom_line(aes(r,
      pcum_cugs,
      colour = "cumsum cugs"
    )) +
    theme_bw() +
    labs(
      x = "r (PCR duplicates)",
      y = "proportion"
    ) +
    geom_vline(
      data = glm_estimates,
      mapping = aes(xintercept = max_r),
      linetype = "dashed"
    ) +
    scale_color_manual(
      name = "",
      values = c(
        "darkgrey",
        "red",
        "blue",
        "coral"
      )
    ) +
    xlim(0, x_lim) +
    theme(
      axis.title.x = element_text(size = rel(1.4)),
      axis.title.y = element_text(size = rel(1.4)),
      axis.text = element_text(size = rel(1.3)),
      legend.title = element_text(size = rel(1.1), face = "bold"),
      legend.text = element_text(size = rel(1.1))
    )

  p2 <- ggplot(chimera_counts) +
    geom_point(aes(r,
      n_chimeras,
      colour = "observed"
    ),
    size = 1
    ) +
    geom_line(aes(r,
      n_cugs * phat_chimeras,
      colour = "predicted"
    )) +
    theme_bw() +
    labs(
      x = "r",
      y = "count"
    ) +
    scale_color_manual(
      name = "",
      values = c("darkgrey", "blue")
    ) +
    annotate(
      "text",
      x = x_lim - 5,
      y = max(chimera_counts$n_chimeras) - 150,
      label = sprintf(
        "SIHR=%.4f",
        glm_estimates$SIHR
      ),
      size = 3,
      hjust = 1,
      vjust = 1
    ) +
    xlim(0, x_lim) +
    theme(
      axis.title.x = element_text(size = rel(.7)),
      axis.title.y = element_text(size = rel(.7))
    ) +
    theme(legend.position = "none")

  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")


  p_fit <-
    ggdraw() +
    draw_plot(p1, 0, 0, 1, 1) +
    draw_plot(p2, 0.65, 0.18, .3, .3) +
    draw_label(dataset_name,
      x = 0,
      y = 1,
      vjust = 2,
      hjust = -1.2,
      size = 12,
      fontface = "italic"
    )


  return(list(p = p_fit, legend = legend))
}





plot_fp_reduction <- function(data_list,
                     dataset_name) {
  
 outcome_counts <- data_list$outcome_counts
  
  
 p1 <-  ggplot(outcome_counts) + 
    
    geom_line(
      aes(x = FP,
          y = FN),
      colour="grey")+
    geom_point(
      aes(x = FP,
          y = FN),
      size=1)+
    labs(x="False Positives",
         y="False Negatives") +
    geom_point(data=data_list$summary_stats$cutoff_dt %>% 
                 filter(approach %in% c( "discard_torc", "no_discarding", "no_purging")),
               aes(x = FP,
                   y = FN,
                   shape=approach),
               size=2,
               colour="darkblue")+
    scale_x_sqrt() +
    scale_y_sqrt() +
    theme_bw()  +
   theme(
     legend.title = element_text(face = "bold")
   )
    
  
  p1 <- ggdraw() +
    draw_plot(p1, 0, 0, 1, 1) +
    draw_label(dataset_name,
               x = 0,
               y = 1,
               vjust = 2,
               hjust = -1.3,
               size = 12,
               fontface = "italic"
    )
  
  
  return( p1)
}


plot_fp_tradoff <- function(data_list,
                              dataset_name) {
  
  outcome_counts <- data_list$outcome_counts
  
  
  p1 <-  
    ggplot(data_list$outcome_counts) + 
    geom_point(
      aes(x = FPm,
          y = FNm),
      size=.5)+
    labs(x="Marginal Decrease in False Positives (reduce phantom molecs) ",
         y="Marginal Increase in False Negatives (discard real molecs)") +
    geom_point(data=data_list$summary_stats$cutoff_dt %>% 
                 filter(approach %in% c( "discard_torc")),
               aes(x = FPm,
                   y = FNm),
               size=2,
               colour="purple") +
    geom_line(
      aes(x = FPm,
          y = FPm,
          colour="1")
    ) +
    geom_line(
      aes(x = FPm,
          y = 2*FPm,
          colour="2"))+
    geom_line(
      aes(x = FPm,
          y = 3*FPm,
          colour="3"))+
    geom_line(
      aes(x = FPm,
          y = 4*FPm,
          colour="4"))+
    geom_line(
      aes(x = FPm,
          y = 5*FPm,
          colour="5"))+
    geom_line(
      aes(x = FPm,
          y = 9*FPm,
          colour="9"))+
    scale_y_log10() +
    theme_bw()  +
      theme(
        legend.title = element_text(face = "bold")) + 
      scale_color_discrete(name = "TORC") 
    
  
  
  p1 <- ggdraw() +
    draw_plot(p1, 0, 0, 1, 1) +
    draw_label(dataset_name,
               x = 0,
               y = 1,
               vjust = 2,
               hjust = -1.8,
               size = 12,
               fontface = "italic"
    )
  
  
  return(p1)
}

make_plots <- function(data_list, dataset_name, x_lim = 160) {
  p_rdist <- plot_read_distributions(data_list$reads_dist_summary)


  p_fit <- plot_fit(data_list$fit_out$chimera_counts,
    data_list$fit_out$glm_estimates,
    x_lim = x_lim
  )


  p_class <- plot_posterior_prob(data_list$outcome_counts,
    # data_list$fit_out$glm_estimates,
    data_list$optimal_cutoff,
    x_lim = x_lim,
    y_lim = 1
  )

  p_phantoms <- plot_phantoms(
    data_list$umi_counts_cell,
    data_list$umi_counts_gene,
    data_list$umi_counts_sample
  )


  plot_list <- list(
    p_rdist = p_rdist,
    p_fit = p_fit,
    p_class = p_class,
    p_phantoms = p_phantoms
  )
  saveRDS(
    plot_list,
    file.path(
      output_dir,
      sprintf(
        "%s_ggplot_objects.rds",
        dataset_name
      )
    )
  )


  return(plot_list)
}

plot_posterior_prob <- function(data_list,
                                dataset_name) {
  
  outcome_counts <-data_list$outcome_counts
  
  optimal_cutoff <-
    data_list$summary_stats$cutoff_dt %>%
    filter(cutoff == "torc")   %>%
    mutate(
      qs = -10 * log10(1e-16) + 10 * log10(qr + 1e-16)
    )
  
  p1 <- ggplot(outcome_counts) +
    geom_line(aes(
      x = qs,
      y = FPp,
      colour = "FPp"
    ),
    alpha = 0.7
    )  +
    geom_point(aes(
      x = qs,
      y = o,
      colour = "ECDF (o)"
    ),
    size = 1,
    shape = 18
    )  +
    geom_vline(
      data = optimal_cutoff,
      aes(xintercept = qs),
      linetype = "dashed",
      color = "grey"
    ) +
    annotate(
      "text",
      x = optimal_cutoff$qs + 5,
      y = .1,
      label = sprintf(
        "qs*=%.1f",
        optimal_cutoff$qs
      ),
      size = 4
    ) +
    labs(
      x = "qs",
      y = "proportion"
    ) +
    theme_bw() +
    scale_color_manual(
      name = "",
      values = c(
        "red",
        "green",
        "darkblue",
        "purple"
      )
    ) +
    theme(
      axis.title.x = element_text(size = rel(1.4)),
      axis.title.y = element_text(size = rel(1.4)),
      axis.text = element_text(size = rel(1.3)),
      legend.title = element_text(size = rel(1.1), face = "bold"),
      legend.text = element_text(size = rel(1.1))
    )
  
  legend <- get_legend(p1)
  p1 <- p1 + theme(legend.position = "none")
  
  p_fit <-
    ggdraw() +
    draw_plot(p1, 0, 0, 1, 1) +
    draw_label(dataset_name,
               x = 0,
               y = 1,
               vjust = 2,
               hjust = -1.3,
               size = 12,
               fontface = "italic"
    )
  
  return(list(p = p_fit, legend = legend))
}

# save_reassigned_read_counts <- function(read_counts,
#                                         data_outcome,
#                                         S,
#                                         reassigned_counts_filepath,
#                                         sample_names) {
#   data_outcome <-
#     data_outcome  %>%
#     select(c("outcome", paste0("t", 1:S))) %>%
#     setNames(c("outcome", sample_names))
#
#   reassigned_read_counts <-
#     join_data(read_counts %>%
#                 select(outcome, cell, umi, gene),
#               data_outcome)
#
#   saveRDS(reassigned_read_counts, reassigned_counts_filepath)
#
# }
# marg_p_outcome =get_product_logsum(p_outcome,nn_prop),
# cum_p=cumsum(marg_p_outcome),
# pcum_phantom = cumsum(get_product_logsum((k_chimera - 1) * q +
#                                            (k_chimera) * (1 - q),
#                                          p_outcome_m
# )),
# pcum_reads = cumsum(r * n) / sum(r * n),
# pcum_phantom = last(pcum_phantom) - pcum_phantom

# read_counts <-
#   read_counts %>%
#   mutate(outcome = pmap(select(., sample_names),
#                                str_c,
#                                sep=",",
#                                collapse = ""))
###########################################################


#
# compute_prop_nonmissingness <- function(outcome_counts) {
#   p_nonmiss_dt <-
#     outcome_counts %>%
#     group_by(r) %>%
#     slice(1) %>%
#     select(r, p_nonmiss)
#
#   return(p_nonmiss_dt)
# }
#
# update_reads_dist_summary <- function(reads_dist_summary, p_nonmiss_dt, g) {
#   reads_dist_summary$summary_stats <-
#     reads_dist_summary$summary_stats %>%
#     mutate(g = g) %>%
#     select(n_cugs, p_chimeras, g, u, everything())
#
#   reads_dist_summary$conditional <- left_join(reads_dist_summary$conditional,
#     p_nonmiss_dt,
#     by = "r"
#   )
#
#   return(reads_dist_summary)
# }

################################################