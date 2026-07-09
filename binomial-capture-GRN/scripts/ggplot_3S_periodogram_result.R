source(here::here("R", "0preambule.R"))
source(here::here("R", "simulation_periodogram_functions.R"))

base_dir <- here::here("DATA", "3S_periodogram")
p_values <- c("p_0.5", "p_0.25", "p_0.1")

read_periodogram_csv <- function(path, type) {
  df <- read.csv(path, check.names = FALSE)
  per_df <- data.frame(
    freq  = df[["periodogram_df.freq"]],
    power = df[["periodogram_df.power"]],
    type  = factor(type),
    check.names = FALSE
  )
  list(
    per_df   = per_df,
    f_star   = df[["peak_freq"]][1],
    T_star   = df[["peak_period"]][1]
  )
}

make_plot_for_p <- function(pname) {

  act   <- read_periodogram_csv(file.path(base_dir, "act_periodogram_df.csv"), "actual")
  obs   <- read_periodogram_csv(file.path(base_dir, pname, "obs_periodogram_df.csv"), "observed")
  map   <- read_periodogram_csv(file.path(base_dir, pname, "map_periodogram_df.csv"), "mapped")
  
  levs <- c("actual","observed","mapped")
  
  common.period_df <- rbind(act$per_df, obs$per_df, map$per_df)
  common.period_df$type <- factor(common.period_df$type, levels = levs)

  common.peaks_df <- data.frame(
    type   = factor(levs, levels = levs),
    f_star = c(act$f_star, obs$f_star, map$f_star),
    T_star = c(act$T_star, obs$T_star, map$T_star),
    y_max  = c(
      max(act$per_df$power, na.rm = TRUE),
      max(obs$per_df$power, na.rm = TRUE),
      max(map$per_df$power, na.rm = TRUE)
    )
  )
  common.peaks_df$label <- paste0(
    "Peak: f ≈ ", signif(common.peaks_df$f_star, 3),
    "\n Period: T ≈ ",   signif(common.peaks_df$T_star, 3)
  )
  
  ggplot(common.period_df, aes(x = freq, y = power, color = type)) +
    geom_vline(data = common.peaks_df, aes(xintercept = f_star),
               color = "black", linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
    # geom_text(data = common.peaks_df, aes(x = f_star, y = y_max, label = label),
    #           size = 3.75, show.legend = FALSE) +
    geom_line(data = subset(common.period_df, type %in% c("actual","mapped")), linewidth = 1) +
    geom_point(data = subset(common.period_df, type == "observed"), size = 1.5) +
    geom_line(data = subset(common.period_df, type == "observed"), lty = "dashed") +
    labs(title = NULL,
         x = "Frequency (cycles/time unit), f",
         y = "Variance-normalized power spectrum density") +
    coord_cartesian(xlim = c(-0.025, 0.75),
                    ylim = c(-0.5, max(common.period_df$power, na.rm = TRUE)+0.5),
                    expand = F) +
    scale_color_manual(
      name = NULL, aesthetics = c("colour","fill"),
      values = c(observed = "#79CEDB", mapped = "#873DC2", actual = "#f78b59"),
      labels = c(observed = "Observed via kernel",
                 mapped   = "Model with mapped parameters",
                 actual   = "True model")
    ) +
    gg_base0
}

p1 <- make_plot_for_p(pname = p_values[1]) +
  labs(y = NULL)
p2 <- make_plot_for_p(pname = p_values[2]) +
  labs(y = NULL) +
  theme(axis.title.y = element_blank())
p3 <- make_plot_for_p(pname = p_values[3]) +
  labs(y = NULL) +
  theme(axis.title.y = element_blank())

combined <- ggarrange(
  p1, p2, p3,
  ncol = 3, nrow = 1,
  common.legend = TRUE, 
  legend = "bottom"      
)

annotate_figure(
  combined,
  left = text_grob("Variance-normalized spectrum", rot = 90)
)
#8.67

