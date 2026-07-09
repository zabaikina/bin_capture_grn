base_dir <- here::here("DATA", "3S_periodogram")
p_values <- c("p_0.5", "p_0.25", "p_0.1")

make_density_plot_for_p <- function(pname, bins = 80, xlim = NULL) {
  act.for_density <- fread(file.path(base_dir, "act_for_density.csv"))
  map.for_density <- fread(file.path(base_dir, pname, "map_for_density.csv"))
  obs.for_density <- fread(file.path(base_dir, pname, "obs_for_density.csv"))
  
  levs <- c("actual","observed","mapped")
  act.for_density[, type := factor(type, levels = levs)]
  map.for_density[, type := factor(type, levels = levs)]
  obs.for_density[, type := factor(type, levels = levs)]
  
  ggplot() + 
    geom_histogram(
      data = obs.for_density,
      aes(x = concentration, y = after_stat(density), fill = type),
      bins = bins, position = "identity", alpha = 1, colour = NA
    ) +
    stat_density(
      data = map.for_density,
      aes(x = concentration, y = after_stat(density), colour = type),
      geom = "line", position = "identity", linewidth = 1
    ) +
    # stat_density(
    #   data = act.for_density,
    #   aes(x = concentration, y = after_stat(density), colour = type),
    #   geom = "line", position = "identity", linewidth = 1
    # ) +
    labs(
      title = paste0("p = ", sub("^p_", "", pname)),
      x = "Protein concentration, x",
      y = "Protein distribution"
    ) +
    coord_cartesian(xlim = c(0,70), ylim = c(0, 0.151), expand = FALSE) +
    scale_color_manual(
      name = NULL,
      aesthetics = c("colour", "fill"),
      values = c(observed = "#79CEDB", mapped = "#873DC2", actual = "#f78b59"),
      labels = c(
        observed = "Observed via kernel",
        mapped   = "Model with mapped parameters",
        actual   = "True model"
      )
    ) +
    gg_base0
}

g1 <- make_density_plot_for_p(pname = p_values[1]) +
  labs(y = NULL)
g2 <- make_density_plot_for_p(pname = p_values[2]) +
  labs(y = NULL) +
  theme(axis.title.y = element_blank())
g3 <- make_density_plot_for_p(pname = p_values[3]) +
  labs(y = NULL) +
  theme(axis.title.y = element_blank())

combined <- ggarrange(
  g1,g2,g3,p1, p2, p3,
  ncol = 3, nrow = 2,
  legend = "none",
  heights = c(0.85, 1)
)

annotate_figure(
  combined,
  left = text_grob("Variance-normalized spectrum Probability distribution", 
                   size = 10, rot = 90)
)
