source(here::here("R", "ggplot_ggbase0.R"))

# BASE FOR PLOTS ---------------------------------------------------------------
ymin10powlog <- -4

gg_base <- ggplot() + gg_base0 + 
  labs(x = "Protein concentration, x", y = "", 
                title = "Observed vs. mapped protein concentration distributions",
                fill = "", color = ""
       )+
  scale_color_manual(
    aesthetics = c("colour", "fill"),
    values = c(observed = "#79CEDB", mapped = "#873DC2"),
    labels = c(observed = "Observed via kernel \u03C6", 
               mapped = "Model with mapped parameters")
  ) +
  guides(fill   = guide_legend(nrow = 1, byrow = TRUE),
         colour = guide_legend(nrow = 1, byrow = TRUE))

#-------------------------------------------------------------------------------

# DRAW ALL DATA IN GG_DF_1S ----------------------------------------------------
gg_base + 
  geom_histogram(
    data = subset(gg_df_1S, kind %in% c("observed")),
    aes(x = val, y = after_stat(density), fill = kind),
    bins = 80, position = "identity", alpha = 1, colour = NA
  ) +
  geom_line(              #logspline
    data = ls_df_1S,
    aes(x = x, y = y, colour = kind,
        group = interaction(p,t) ),
    linewidth = 0.8,
  ) +
  # stat_bin(           #compare to histogram
  #   data = subset(gg_df_1S, kind == "observed"),   
  #   aes(y = after_stat(density), fill = kind, colour = kind),
  #   bins = 80, geom = "area", position = "identity",
  #   color = "#000000", fill = NA,
  #   pad = TRUE, alpha = 0.15, linewidth = 0.7
  # ) +
  # stat_density(         #compare logspline to the built-in density
  #   data = subset(gg_df_1S, kind == "observed"),
  #   aes(y = after_stat(density), colour = kind),
  #   fill = NA, alpha = 0.75, n = 200, color = "green",
  #   kernel = "cosine",
  #   linewidth = 0.7, adjust = 1
  #) +
  facet_grid(rows = vars(p), cols = vars(t), 
             scales = "fixed",
             switch = "y"
  ) +
  geom_vline(
    data = means_1S,
    aes(xintercept = mean, colour = kind),
    linewidth = 0.5, linetype = "dashed", alpha = 1
  ) +
  geom_text(
    data = means_1S,
    aes(x = mean, y = Inf, label = sprintf("μ = %.1f", mean),
        colour = kind, vjust = vjust),
    nudge_x = 50,
    size = 2.8, show.legend = FALSE
  ) +
  # scale_y_continuous(
  #   trans   = pseudo_log_trans(base = 10, sigma = 10^ymin10powlog),
  #   breaks  = 10^(ymin10powlog:0),
  #   labels  = function(b) parse(text = paste0("10^", round(log10(b)))),
  #   position = "right"
  # ) +
  scale_y_continuous(  #for normal (linear) plot
    breaks = seq(0, 0.07, by = 0.01),
    position = "right"
  ) +
  coord_cartesian(ylim = c(0, 0.046), 
                  xlim = c(10, 290),
                  expand = F) 

#-------------------------------------------------------------------------------



# DRAW ALL DATA IN GG_DF_2S ----------------------------------------------------

gg_base + 
  geom_histogram(
    data = subset(gg_df_2S, kind %in% c("observed")),
    aes(x = val, y = after_stat(density), fill = kind),
    binwidth = 1, position = "identity", alpha = 1, colour = NA
  ) +
  geom_line(
    data = ls_df_2S,
    aes(x = x, y = y, colour = kind,
        group = interaction(kind,p,t) ),
    linewidth = 0.8
  )+
  # stat_bin(        #compare to histogram
  #   data = subset(gg_df_2S, kind == "observed"),
  #   aes(y = after_stat(density), fill = kind, colour = kind),
  #   bins = 90, geom = "area", position = "identity",
  #   color = "#000000",
  #   pad = TRUE, alpha = 0.15, linewidth = 0.7
  # ) +
  # geom_density(    #compare to the built-in density
  #   data = subset(gg_df_2S, kind == "observed"),
  #   aes(y = after_stat(density), colour = kind), color = "red",
  #   fill = NA, alpha = 0.75,
  #   linewidth = 0.7, adjust = 1
  # ) +
  facet_grid(rows = vars(p), cols = vars(t), 
             scales = "free",
             switch = "y"
  ) +
  scale_y_continuous(
    trans   = pseudo_log_trans(base = 10, sigma = 10^ymin10powlog),
    breaks  = 10^(ymin10powlog:0),
    labels  = function(b) parse(text = paste0("10^", round(log10(b)))),
    position = "right"
  ) +
  # scale_y_continuous(  #for normal (linear) plot
  #   breaks = seq(0, 0.07, by = 0.01),
  #   position = "right"
  # ) +
  coord_cartesian(xlim = c(-2,300), 
                  ylim = c(10^(-5), 1), #for logplot
                  expand = F) 
#-------------------------------------------------------------------------------


# DRAW GG_DF_2S ONLY SELECTED VALUES OF P ------------------------------------------------

ggplot() + gg_base + 
  geom_histogram(
    data = subset(gg_df_2S_chop, kind %in% c("observed")),
    aes(x = val, y = after_stat(density), fill = kind),
    binwidth = 1, position = "identity", alpha = 1, colour = NA
  ) +
  geom_line(
    data = ls_df_2S_chop,
    aes(x = x, y = y, colour = kind,
        group = interaction(kind,p,t) ),
    linewidth = 0.8
  )+
  facet_grid(rows = vars(p), cols = vars(t), 
             scales = "fixed",
             switch = "y"
  ) +
  # scale_y_continuous(
  #   trans   = pseudo_log_trans(base = 10, sigma = 10^ymin10powlog),
  #   breaks  = 10^(ymin10powlog:0),
  #   labels  = function(b) parse(text = paste0("10^", round(log10(b)))),
  #   position = "right"
  # ) +
  scale_y_continuous(  #for normal (linear) plot
    n.breaks = 4,
    position = "right"
  ) +
  coord_cartesian(xlim = c(-2,300), 
                  ylim = c(0, 0.06), #for logplot
                  #ylim = c(0, 0.13),    #for normal (linear) plot
                  expand = F) 
#7.41 x 5.14
#-------------------------------------------------------------------------------




