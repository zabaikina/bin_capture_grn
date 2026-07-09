gg_base0 <- theme_minimal() +
  theme(
    title = element_text(colour = "#000000", size = 10),
    panel.spacing.x      = unit(4, "mm"),  # space between columns
    panel.spacing.y      = unit(4,  "mm"),  # space between rows
    strip.background.y   = element_rect(fill = "#F2F2F2", colour = NA), # rows (p)
    strip.text.y         = element_text(colour = "#000000", size = 10),
    strip.background.x   = element_rect(fill = "#F2F2F2", colour = NA), # cols (t)
    strip.text.x         = element_text(colour = "#000000", size = 10),
    strip.placement      = "outside",
    panel.border         = element_rect(colour = "#9e9e9e", fill = NA, linewidth = 1),
    axis.text            = element_text(colour = "#000000", size = 8),
    legend.position      = "bottom",
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.text          = element_text(colour = "#000000", size = 10),
  )