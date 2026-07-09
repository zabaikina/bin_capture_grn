source(here::here("R", "0preambule.R"))

# GENERAL PARAMS ---------------------------------------------------------------

runtime <- 25
rec_times <- c(0.25, 1.5)
y0_mean <- 60; y0_sd <- 5
n_sims <- 10^5
ssize <<- 1

k    <- c(1)      # decay rates
pars <- list(
  rho_u   = c(12), 
  rho_b   = c(52), 
  sigma_u = c(1000), 
  sigma_b = c(1), 
  beta = c(4), 
  n = matrix( c(3), nrow = 1, byrow = T)
)

y_mean_meanfield(pars, k)

p_vals_1S   <- c(1, 0.75, 0.5)              # e.g., c(0.2, 0.5, 0.8)
p_labels_1S <- paste0("p=", format(p_vals_1S)) # labels for facets/legend

# COMPARISON FOR P_i -----------------------------------------------------------


#actual conc Y
tic()
act_all_1S <- simulate_many_mcl(N = ssize, n_sims = n_sims,
  k = k, alpha = alpha, burst.size = burst.size,
  Tmax = runtime, y0_mean = y0_mean, y0_sd = y0_sd, pars = pars, seed = NULL,
  rec_times = rec_times); toc()

ggplot(act_all_1S, aes(x = y_actual, fill = species)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ t, nrow = 1) +
  geom_vline(
    data = aggregate(y_actual ~ t+species, data = act_all_1S, FUN = mean, na.rm = TRUE), 
    aes(xintercept = y_actual),
    linewidth = 0.6, linetype = "dashed") +
  geom_text(data = aggregate(y_actual ~ t+species, data = act_all_1S, FUN = mean, na.rm = TRUE),
            aes(x = y_actual, y = Inf, #top
                label = sprintf("μ = %.1f", y_actual)), #position of the label
            vjust = 5, hjust = 0, nudge_x = 10,
            size = 3, show.legend = FALSE)+
  labs(x = "concentration (y_1)", y = "count",
       title = "Distribution of y_1 across times") +
  theme_minimal()

#observed conc X through norm kernel
obs_list_1S <- vector("list", length = length(p_vals_1S))
for (i in seq_along(p_vals_1S)) {
  tmp <- observe_x_normal(act_all_1S, p = c(y_1 = p_vals_1S[i]))
  tmp$p <- p_labels_1S[i]
  obs_list_1S[[i]] <- tmp
}
obs_all_1S <- do.call(rbind, obs_list_1S)
obs_all_1S$p <- factor(obs_all_1S$p, levels = p_labels_1S)

#mapped X (original model w rescaled params)
map_list_1S <- vector("list", length = length(p_vals_1S))
for (i in seq_along(p_vals_1S)) {
  pars_map <- pars
  pars_map$beta <- pars$beta * p_vals_1S[i]
  pars_map$sigma_u <- pars$sigma_u * 
    as.numeric(ps_for_sigma_u_mapping(N = ssize, pars, p_vals_1S[i]))
  
  y0_sd_map <- as.numeric(sd_for_normIC_mapping(y0_mean, y0_sd, p_vals_1S[i]))
  
  tmp <- simulate_many_mcl(N = ssize, n_sims = n_sims,
    k = k, alpha = alpha, burst.size = burst.size,
    Tmax = runtime, y0_mean = y0_mean*p_vals_1S[i], y0_sd = y0_sd_map,
    pars = pars_map, seed = NULL, rec_times = rec_times)
  names(tmp) <- c("t","species", "x_mapped")
  tmp$p <- p_labels_1S[i]
  
  map_list_1S[[i]] <- tmp
}
map_all_1S <- do.call(rbind, map_list_1S)
map_all_1S$p <- factor(map_all_1S$p, levels = p_labels_1S)


gg_df_1S <- rbind(
  data.frame(t = map_all_1S$t, species = map_all_1S$species, 
             kind = "mapped",   val = map_all_1S$x_mapped,   p = map_all_1S$p),
  data.frame(t = obs_all_1S$t, species = obs_all_1S$species, 
             kind = "observed", val = obs_all_1S$x_observed, p = obs_all_1S$p)
)
gg_df_1S$kind <- factor(gg_df_1S$kind, levels = c("observed","mapped"))

#write.csv(gg_df_1S, here::here("DATA", "gg_df_1S_S1_Ru12_Rb52_Su1000_B4_n3.csv"), row.names = F)
#gg_df_1S <- read.csv(here::here("DATA", "gg_df_1S_S1_Ru12_Rb52_Su1000_B4_n3.csv"))

means_1S <- ggplot_means_long(gg_df_1S, include = c("observed","mapped"), 
                         vjust = 1.5, p = p_vals_1S)
ls_df_1S <- make_logspline_df(gg_df_1S, kind_name = "mapped",
                           lbound = 0, ubound = max(gg_df_1S$val),
                           n = 80, nkn = 1)

# # Simulate and draw one trajectory ---------------------------------------------
# # ymf <- y_mean_meanfield(pars, k)
# # y0_traj <- 100
# # 
# # res <- simulate_once_traj(k, alpha, burst.size, 
# #                           Tmax = 10, 
# #                           y0 = y0_traj,
# #                           pars, dt = 0.02)
# # 
# # ggplot(res, aes(t, y_actual, color = species)) +
# #   geom_line(linewidth = 0.6) +
# #   labs(x = "time", y = "concentration", color = "species") +
# #   coord_cartesian(ylim = c(0, max(res$y_actual)))+
# #   theme_minimal()
# # 
