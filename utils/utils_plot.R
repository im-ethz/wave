# Coefficients of explanatory analysis
get_varlabel <- function(name) {
  if (name == "HEART_RATE") {
    expression(HR[])
  } else if (name == "HRV") {
    expression(HRV[])
  } else if (name == "RESPIRATION") {
    expression(RF[])
  }
}

get_univar_plot <- function(mtype) {
  mtype_plot <- m_fits_plotdata %>%
    filter(measure_type == mtype) %>%
    mutate(measure_type_label = case_when(
      measure_type == "HEART_RATE" ~ "Heart rate",
      measure_type == "HRV" ~ "Heart rate variability",
      measure_type == "RESPIRATION" ~ "Respiration frequency"
    )) %>%
    ggplot(aes(y = .variable, x = .value)) +
    stat_halfeye(aes(
      fill = stat(ifelse(cdf < 0.025 | cdf > 0.975, "a", x < 0)),
      alpha = stat(f),
    ),
    .width = c(.8, .95),
    point_interval = mean_qi, orientation = "horizontal",
    point_size = 1.4, interval_size_range = c(0.5, 1.1)
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(name = "Interval", values = c(NA, "#E64B35FF", "#4DBBD5FF"), na.translate = FALSE) +
    scale_y_discrete(name = get_varlabel(mtype), labels = function(x) {
      return(str_remove(x, "HEART_RATE_|RESPIRATION_|hrv_"))
    }) +
    scale_x_continuous(name = "Standardized coefficient") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6),
      text = element_text(family = "Times New Roman", size = 8)
    ) +
    coord_cartesian(xlim = c(-1.4, 1.4)) #+ facet_wrap(vars(measure_type_label))

  if (mtype == "HEART_RATE") {
    mtype_plot <- mtype_plot +
      annotate(geom = "text", label = expression(HR[]), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, fontface = "bold", size = 5)
  }
  if (mtype == "HRV") {
    mtype_plot <- mtype_plot +
      annotate(geom = "text", label = expression(HRV[]), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, fontface = "bold", size = 5)
  }
  if (mtype == "RESPIRATION") {
    mtype_plot <- mtype_plot +
      annotate(geom = "text", label = expression(RF[]), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, fontface = "bold", size = 5)
  }

  return(mtype_plot)
}

# PCA loadings in risk score
get_loadings_plot <- function(df) {
  loadings_plot <- df %>%
    ggplot(aes(y = feature.name, x = value, fill = direction_influence, alpha = abs(value))) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(measure_type ~ PC, scales = "free_y", space = "free_y") +
    scale_y_discrete(position = "left") +
    scale_x_continuous(breaks = c(-0.5, 0, 0.5)) +
    scale_fill_manual(values = c("#4DBBD5FF", "#E64B35FF")) +
    scale_alpha_continuous(range = c(0.6, 1), trans = "pseudo_log") +
    theme_bw() +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none",
      strip.background.x = element_blank(),
      strip.text.x = element_blank(),
      strip.background.y = element_rect("white"),
      text = element_text(family = "Times New Roman")
    ) +
    xlab("Loading") +
    coord_cartesian(xlim = c(-0.65, 0.65))

  return(loadings_plot)
}

# Coefficients of risk score model
get_multivar_plot <- function(df) {
  plotdata <- df %>%
    gather_draws(`b_.*`, regex = T) %>%
    mutate(.variable = str_remove(.variable, "b_")) %>%
    mutate(.variable = str_replace(.variable, "sexmale", "Sex\n(male)")) %>%
    mutate(.variable = str_replace(.variable, "age", "Age")) %>%
    mutate(.variable = str_replace(.variable, "rel_day", "Current\nlength of stay")) %>%
    filter(!str_contains(.variable, "Intercept")) %>%
    filter(str_contains(.variable, "PC")) %>%
    mutate(.variable = factor(.variable, ordered = T, levels = c("Age", "Sex\n(male)", "Current\nlength of stay", paste0("PC", 1:100))))

  return(plotdata %>%
    ggplot(aes(x = .variable, y = .value)) +
    stat_eye(aes(fill = stat(ifelse(cdf < 0.025 | cdf > 0.975, "a", y < 0)), alpha = stat(f)),
      orientation = "vertical", .width = c(0.8, 0.95), point_interval = mean_qi
    ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~.variable, nrow = 1, scales = "free_x") +
    scale_fill_manual(name = "Interval", values = c(NA, "#E64B35FF", "#4DBBD5FF"), na.translate = FALSE) +
    coord_cartesian(ylim = c(-1.5, 2)) +
    theme_bw() +
    ylab("\nStandardized coefficient\n ") +
    xlab("PC") +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "white"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      text = element_text(family = "Times New Roman")
    ))
}

# Probabilities of risk score
inv_cloglog <- function(x) 1 - exp(-exp(x))
discharge_prob <- function(pred, tau1, tau2) inv_cloglog(tau1 - pred)
stay_prob <- function(pred, tau1, tau2) inv_cloglog(tau2 - pred) - inv_cloglog(tau1 - pred)
icu_prob <- function(pred, tau1, tau2) 1 - inv_cloglog(tau2 - pred)

get_probability_plot <- function(df) {
  main_graph <- df %>%
    ggplot(aes(x = pred, y = value)) +
    stat_lineribbon(aes(alpha = stat(level)), fill = discharge_color, color = discharge_color, .width = c(.95, 0), data = df %>% filter(name == "discharge_prob"), point_interval = mean_qi) +
    stat_lineribbon(aes(alpha = stat(level)), fill = icu_color, color = icu_color, .width = c(.95, 0), data = df %>% filter(name == "icu_prob"), point_interval = mean_qi) +
    scale_alpha_discrete(name = "alpha", range = c(0.4, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    xlab("Risk score") +
    ylab("Probability of outcome") +
    theme_bw() +
    theme(legend.position = "none", text = element_text(family = "Times New Roman"))

  graph_legend <- get_legend(
    data.frame(a = 1:2, b = 1:2, c = c("Hospital discharge", "ICU admission")) %>%
      ggplot(aes(x = a, y = b, fill = c)) +
      geom_area() +
      scale_fill_manual(name = "Outcome", values = c(discharge_color, icu_color)) +
      theme(legend.position = "bottom") +
      guides(
        fill = guide_legend(nrow = 1, label.hjust = 1)
      )
  )

  total_plot <- plot_grid(plot_grid(white_plot, graph_legend, white_plot, ncol = 3, rel_widths = c(0.225, 0.5, 0.225)), main_graph, ncol = 1, rel_heights = c(0.8, 7))

  return(total_plot)
}

# Performance evaluation
get_auc_plot_multi <- function(aucdata1, aucdata2, color1, fill1, color2, fill2) {
  ggplot() +
    geom_ribbon(aes(
      x = aucdata1[["auc_smooth"]]$x,
      ymin = aucdata1[["auc_smooth"]]$nne - 1.96 * sqrt(aucdata1[["auc_smooth"]]$var),
      ymax = pmin(1, aucdata1[["auc_smooth"]]$nne + 1.96 * sqrt(aucdata1[["auc_smooth"]]$var))
    ), fill = fill1) +
    geom_point(aes(x = aucdata1[["auc"]]$time, y = aucdata1[["auc"]]$mean.rank), color = color1, pch = 16, size = 2) +
    geom_line(aes(x = aucdata1[["auc_smooth"]]$x, aucdata1[["auc_smooth"]]$nne), color = color1, size = 1.2) +
    geom_ribbon(aes(
      x = aucdata2[["auc_smooth"]]$x,
      ymin = aucdata2[["auc_smooth"]]$nne - 1.96 * sqrt(aucdata2[["auc_smooth"]]$var),
      ymax = pmin(1, aucdata2[["auc_smooth"]]$nne + 1.96 * sqrt(aucdata2[["auc_smooth"]]$var))
    ), fill = fill2) +
    geom_point(aes(x = aucdata2[["auc"]]$time, y = aucdata2[["auc"]]$mean.rank), color = color2, pch = 21, size = 2) +
    geom_line(aes(x = aucdata2[["auc_smooth"]]$x, aucdata2[["auc_smooth"]]$nne), color = color2, size = 1.2, linetype = "dashed") +
    geom_hline(yintercept = 0.5, linetype = "dashed") +
    theme_bw() +
    xlab("Length of hospital stay (days)") +
    ylab("AUROC") +
    scale_x_continuous(breaks = 1:10, expand = expansion(add = 0.05)) +
    scale_y_continuous(breaks = seq(0, 1, 0.1), expand = expansion(add = 0.05))
}

get_legend_multi <- function(auc1_label, auc2_label, color1, fill1, color2, fill2) {
  get_legend(tibble(a = 1:2, `Risk score` = fct_inorder(c(auc1_label, auc2_label), ordered = T)) %>%
    ggplot(aes(x = a, y = a, color = `Risk score`, pch = `Risk score`), size = 3) +
    geom_line(aes(linetype = `Risk score`)) +
    geom_ribbon(aes(ymin = a - 0.5, ymax = a + 0.5, fill = `Risk score`), color = NA) +
    geom_point() +
    scale_color_manual(values = c(color1, color2)) +
    scale_fill_manual(values = c(fill1, fill2)) +
    scale_shape_manual(values = c(16, 21)) +
    theme_bw() +
    theme(legend.position = "top"))
}

plot_auc_multi <- function(aucdata1, auc1_label, aucdata2, auc2_label, color1 = "#3C5488FF", fill1 = "#3C548833", color2 = "#9c9c9cFF", fill2 = "#9c9c9c33") {
  performance_plot <- get_auc_plot_multi(aucdata1, aucdata2, color1, fill1, color2, fill2) + coord_cartesian(ylim = c(0, 1), xlim = c(1, max(aucdata1$auc$time)))
  legend <- get_legend_multi(auc1_label, auc2_label, color1, fill1, color2, fill2)

  plot_grid(plot_grid(white_plot, legend, white_plot, ncol = 3, rel_widths = c(0.16, 0.6, 0.16)), performance_plot, ncol = 1, rel_heights = c(1, 12))
}
