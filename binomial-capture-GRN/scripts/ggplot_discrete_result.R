source(here::here("R", "0preambule.R"))

files_work <- list.files(path = here::here("DATA", "discrete_works"), pattern = "[.]csv$", 
                    full.names = TRUE, include.dirs = F)

df_work <- rbindlist(discr.read_csv_to_df(files_work), use.names = TRUE, fill = TRUE)

files_fail <- list.files(path = here::here("DATA", "discrete_breaks"), pattern = "[.]csv$", 
                         full.names = TRUE, include.dirs = F)

df_fail <- rbindlist(discr.read_csv_to_df(files_fail), use.names = TRUE, fill = TRUE)

pal <- c(
  "1"    = "#F78B59",
  "0.5"  = "#FFE033",
  "0.25" = "#20C7B6",
  "0.1"  = "#1872D9"
)

gg1 <- gg_base_local +
  geom_line(data = df_work[kind == "actual"],
            aes(x = x, y = density, color = p_fac, 
                linetype = kind),
            linewidth=1.2) +
  geom_line(data = df_work[kind == "mapped"],
            aes(x = x, y = density, 
                colour = p_fac, group = p_fac, 
                linetype = kind),
            linewidth = 0.9) +
  geom_line(data = df_work[kind == "observed"],
            aes(x = x, y = density, 
                colour = p_fac, group = p_fac,
                linetype = kind),
            linewidth = 0.4) +
  geom_point(data = df_work[kind == "observed"],
             aes(x = x, y = density, colour = p_fac, 
                 shape = kind),
             fill = "white", stroke = 0.9, size = 1) +
  scale_colour_manual(
    name   = "p = ",
    breaks = names(pal),
    values = pal,
    guide  = guide_legend(order = 1, override.aes = list(
      linetype = "solid", shape = NA, linewidth = c(1.7, 1.2, 1.2, 1.2)))
  )+
  coord_cartesian(xlim = c(-1,80), ylim = c(-0.0025, 0.151), expand = F); gg1

gg2 <- gg_base_local +
  geom_line(data = df_fail[kind == "actual"],
            aes(x = x, y = density, 
                color = p_fac, linetype = kind),
            linewidth=1.2) +
  geom_line(data = df_fail[kind == "mapped"],
            aes(x = x, y = density, 
                colour = p_fac, group = p_fac, 
                linetype = kind),
            linewidth = 0.9) +
  geom_line(data = df_fail[kind == "observed"],
            aes(x = x, y = density, colour = p_fac, group = p_fac,
                linetype = kind),
            linewidth = 0.4) +
  geom_point(data = df_fail[kind == "observed"],
             aes(x = x, y = density, colour = p_fac, 
                 shape = kind),
             fill = "white", stroke = 0.9, size = 1) +
  scale_colour_manual(
    name   = "p = ",
    breaks = names(pal),
    values = pal,
    guide  = guide_legend(order = 1, override.aes = list(
      linetype = "solid", shape = NA, linewidth = 1.1))
  )+
  coord_cartesian(xlim = c(0,20), ylim = c(-0.0025,0.11), expand = FALSE); gg2

ggarrange(
  gg1, gg2,
  ncol = 2, nrow = 1,
  align = "hv",
  common.legend = TRUE,  # one legend for both
  legend = "right"       # place it on the right
)



# FUNCTIONS --------------------------------------------------------------------
discr.parse_file_name <- function(f) {
  file_name <- basename(f)
  file_name <- substr(file_name, start = 1, stop = nchar(file_name)-4)
  kind <- if (str_detect(file_name, "(?i)observed|obs")) "observed"
  else if (str_detect(file_name, "(?i)mapped|map")) "mapped"
  else if (str_detect(file_name, "(?i)true|actual")) "actual"
  else "error"
  #(?i) = case 
  p <- str_match(file_name, "p[=_]?([0-9.]+)")[,2]
  # p is letter, then optionally capture "=" or "_" -- [=_]?, then capture ([0-9.]+) 
  p <- if (!is.na(p)) as.numeric(p) else NA_real_
  list(kind = kind, p = p)
}

discr.read_csv_to_df <- function(files) {
  rlist <- vector("list", length(files))
  
  for (i in seq_along(files)) {
    m <- discr.parse_file_name(files[i])
    dt <- fread(files[i], col.names = c("x", "density"))
    dt$kind <- factor(m$kind)
    dt$p_val <- rep(m$p)
    dt$p_fac <- factor(m$p)
    rlist[[i]] <- dt
  }
  rlist
}

gg_base_local <- ggplot()+
  gg_base0+
  scale_linetype_manual(name="Model",
                        breaks = c("actual","mapped","observed"),
                        values = c(actual = "longdash", mapped = "solid", observed = "dashed")) +
  scale_shape_manual(name = "Model",
                     breaks = c("actual","mapped","observed"),
                     values = c(actual = NA, mapped = NA, observed = 21), drop = FALSE) +
  guides(
    linetype = guide_legend(order = 2, 
                            override.aes = list(
                              colour    = c("#F78B59","grey50","grey50"),
                              linewidth = c(1.7, 1.2, 0.4),
                              shape     = c(NA, NA, 21),
                              fill      = c(NA, NA, "white"),
                              size      = c(NA, NA, 1.4),
                              stroke    = c(NA, NA, 0.9)
                              ) ),
    shape = guide_legend(order = 2)
  ) +
  theme(legend.position = "right")+
  labs(x = "Protein concentration, x", y = "Protein distribution, P(x,t)")
