

files <- list.files(path = here::here("DATA", "discrete_breaks"), 
                    pattern = "[.]csv$", 
                    full.names = TRUE,
                    include.dirs = F); files

parse_name <- function(f) {
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

rlist <- vector("list", length(files))

for (i in seq_along(files)) {
  m <- parse_name(files[i])
  dt <- fread(files[i], col.names = c("x", "density"))
  dt$kind <- factor(m$kind)
  dt$p_val <- rep(m$p)
  dt$p_fac <- factor(m$p)
  rlist[[i]] <- dt
}

df_all <- rbindlist(rlist, use.names = TRUE, fill = TRUE)

ggplot() +
  geom_line(data = df_all[kind=="actual"],
            aes(x = x, y = density),
            colour = "#F78B59", linewidth = 1.2, inherit.aes = FALSE) +
  geom_line(data = df_all[kind=="mapped"],
            aes(x = x, y = density, colour = p_val, group = p_fac),
            linewidth = 0.9) +
  geom_point(data = df_all[kind=="observed"],
             aes(x = x, y = density, colour = p_val),
             shape = 21, fill = "white", stroke = 0.9, size = 1.1) +
  geom_line(data = df_all[kind=="observed"],
            aes(x = x, y = density, colour = p_val, group = p_fac),
            linewidth = 0.5, lty = "dashed") +
  scale_colour_viridis_c(
    option = "viridis",
    begin = 0.3, end = 0.8,
    limits = c(0.1, 0.5),                 # gradient only for p<1
    breaks  = c(0.5, 0.25, 0.1),
    name = "p = "
  ) +
  labs(x = "x", y = "density") +
  coord_cartesian(xlim = c(0,20), ylim = c(-0.0025, 0.12), expand = F)+
  gg_base0






