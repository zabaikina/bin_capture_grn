source(here::here("R", "0preambule.R"))
source(here::here("R", "simulation_periodogram_functions.R"))

# MODEL PARAMETERS SET 1 (FIG) -------------------------------------------------
k    <- c(1, 1, 1)      # decay rates
y0   <- c(95, 75, 100)             #IC
pars <- list(
  rho_u   = rep(130, times = 3),
  rho_b   = rep(2.9, times = 3),
  sigma_u = rep(7*1e4, times = 3),
  sigma_b = rep(1, times = 3),
  beta    = rep(2, times = 3) ,
  n = matrix( c( 0, 0, 3,
                 3, 0, 0,
                 0, 3, 0), nrow = 3, byrow = T)
)
ssize <- 3

# MODEL PARAMETERS SET 2 (TEST) ------------------------------------------------
# k    <- c(1, 1, 1)      # decay rates
# y0   <- c(95, 75, 100)             #IC
# pars <- list(
#   rho_u   = rep(150, times = 3),
#   rho_b   = rep(2.9, times = 3),
#   sigma_u = rep(7*1e4, times = 3),
#   sigma_b = rep(1, times = 3),
#   beta    = rep(10, times = 3) ,
#   n = matrix( c( 0, 0, 3,
#                  3, 0, 0,
#                  0, 3, 0), nrow = 3, byrow = T)
# )
# ssize <- 3

#SIMULATION PARAMETERS ---------------------------------------------------------
Tmax <- 30
nsims <- 10000
dt    <- 0.01
sp <- c("y_1")

# TRAJECTORIES OF ALL SPECIES ---------------------------------------------------

act.sample_traj <- simulate_once_traj(k, alpha, burst.size, Tmax = Tmax, 
                                  y0 = y0, pars, dt = dt)


ggplot(act.sample_traj, aes(t, y_actual, color = species)) +
  geom_line(linewidth = 0.6) +
  labs(x = "time", y = "concentration", color = "species") +
  coord_cartesian(
    #xlim = c(0.8*Tmax, Tmax),
    ylim = c(0, max(act.sample_traj$y_actual))) +
  theme_minimal()

# MEAN and SD -------------------------------------------------------------------
act.stats_df <- simulate_many_grid_stats(
  nsims = nsims,
  k = k, alpha = alpha, burst.size = burst.size,
  Tmax = Tmax, y0 = y0, pars = pars, dt = dt,
  mc.cores = max(1L, parallel::detectCores() - 1L)
)

#draw mean and sd for given species
act.df_sp <- act.stats_df[act.stats_df$species == sp, ]
ggplot(act.df_sp, aes(x = t, y = mean)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.2) +
  geom_line(linewidth = 0.7) +
  labs(title = paste0(sp, ": mean ± sd (", nsims, " sims)"),
       x = "Time, t", y = "Concentration, x") +
  theme_minimal(base_size = 12)

#(useless) chunk that was here is in simulation_result_3S_draft.R

# TRUE MODEL -------------------------------------------------------------------

#these are simple trajectories on grid
#return dataframe: sim (id of sim), t (time), species, y_actual

act.traj_df <- simulate_many_grid_traj(
  nsims = nsims,
  k = k, alpha = alpha, burst.size = burst.size,
  Tmax = Tmax, y0 = y0, pars = pars, dt = dt,
  mc.cores = max(1L, parallel::detectCores() - 1L)
) 


#draw some trajectories for fun
act.some_traj <- act.traj_df[(act.traj_df$sim %in% seq(1,10,by=1)) & (act.traj_df$species == sp), ]
ggplot(act.some_traj, aes(x = t, y = y_actual, group = sim)) +
  geom_line(alpha = 0.15, linewidth = 0.3) +
  labs(title = paste0(sp, ": trajectories (", nsims, " sims)"),
       x = "time", y = "concentration") +
  theme_minimal(base_size = 12)


act.traj_df_edit <- act.traj_df
act.traj_df_edit$species <- "y_1"
act.traj_df_edit$sim <- rep(seq_len(3*nsims), each = 3001)

act.traj_df_edit <- fread(here::here("DATA", "3S_periodogram/act_traj_df.csv"))

#df_sp to list of time series for acf and fft
act.series_from_traj_list <- extract_series_list_fast(act.traj_df_edit, species = sp)

#autocorrelation coefficient 
act.acf_df <- mean_acf_for_species(series_list = act.series_from_traj_list, 
                                   species = sp, dt = dt, max_lag_time = Tmax)

ggplot(act.acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  labs(title = paste0(sp, " — mean autocorrelation across sims"),
       x = "lag (time units)", y = "ACF") +
  theme_minimal(base_size = 12)

# periodogram (mean across simulations)
act.periodogram_df <- mean_periodogram_for_species_slow(series_list = act.series_from_traj_list, 
                                                   species = sp, dt = dt,
                                                   include_0freq = T, demean = T, 
                                                   normalize_area = T)


act.per_df  <- act.periodogram_df$periodogram_df
act.per_df$type <- factor("actual")
act.f_star  <- act.periodogram_df$peak_freq
act.T_star  <- act.periodogram_df$peak_period

ggplot(act.per_df, aes(x = freq, y = power)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = act.f_star, linetype = "dashed") +
  annotate("text", x = act.f_star, y = max(act.per_df$power),
           label = paste0("Peak: f ≈ ", signif(act.f_star, 3),
                          "  (T ≈ ", signif(act.T_star, 3), ")"),
           hjust = -0.1, vjust = 2, size = 3.75) +
  labs(title = paste0(sp, " — mean raw periodogram"),
       x = "Frequency (cycles / time unit)", y = "Power (energy)") +
  coord_cartesian(xlim = c(0,1.5))+
  theme_minimal(base_size = 12)

# act.per_df$omega <- 2*pi*act.per_df$freq   # rad / time-unit
# 
# #decorations
# br_all <- c(0, pi/2, pi, 3*pi/2, 2*pi, 5*pi/2, 3*pi)
# ok <- br_all <= max(act.per_df$omega) + 1e-9
# br <- br_all[ok]
# lab <- expression(0, pi/2, pi, 3*pi/2, 2*pi, 5*pi/2, 3*pi)[ok]
# 
# ggplot(act.per_df, aes(x = omega, y = power)) +
#   geom_line(linewidth = 0.7) +
#   scale_x_continuous(breaks = br, labels = lab) +
#   coord_cartesian(xlim = c(0,2*pi))+
#   geom_vline(xintercept = f_star*2*pi, linetype = "dashed") +
#   annotate("text", x = f_star, y = max(act.per_df$power),
#            label = paste0("Peak: w ≈ ", signif(act.f_star*2*pi, 3),
#                           "  (T ≈ ", signif(act.T_star, 3), ")"),
#            hjust = -1, vjust = 2, size = 3.75) +
#   labs(x = expression("Frequency "~omega~"(rad / time)"), y = "Power")+
#   theme_minimal(base_size = 12)

# OBSERVED MODEL ---------------------------------------------------------------
p <- c(y_1 = 0.1, y_2 = 0.1, y_3 = 0.1)

obs.traj_df <- observe_x_normal(act.traj_df_tmp, p)
obs.keep_col <- c("sim", "t", "species", "x_observed")
obs.traj_df <- obs.traj_df[,obs.keep_col]

obs.df_sp <- obs.traj_df[(obs.traj_df$sim %in% seq(1,10,by=1)) & (obs.traj_df$species == sp), ]
ggplot(obs.df_sp, aes(x = t, y = x_observed, group = sim)) +
  geom_line(alpha = 0.05, linewidth = 0.3) +
  labs(title = paste0(sp, ": trajectories (", nsims, " sims)"),
       x = "time", y = "concentration") +
  theme_minimal(base_size = 12)

obs.traj_df_edit <- obs.traj_df
obs.traj_df_edit$species <- "y_1"
obs.traj_df_edit$sim <- rep(seq_len(3*nsims), each = 3001)

obs.traj_df_edit <- fread(here::here("DATA", "3S_periodogram/p_0.1/obs_traj_df.csv"))

obs.series_from_traj_list <- extract_series_list_fast(obs.traj_df_edit, 
                                                      species = sp, observed = T)

#autocorrelation coefficient 
obs.acf_df <- mean_acf_for_species(series_list = obs.series_from_traj_list, 
                                    species = sp, dt = dt, max_lag_time = Tmax/2)

ggplot(obs.acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  labs(title = paste0(sp, " — mean autocorrelation across sims"),
       x = "Lag (time units)", y = "ACF") +
  theme_minimal(base_size = 12)

obs.periodogram_df <- mean_periodogram_for_species_slow(
                            series_list = obs.series_from_traj_list, 
                            species = sp, dt = dt,
                            include_0freq = T, demean = T, 
                            normalize_area = T)


obs.per_df  <- obs.periodogram_df$periodogram_df
obs.per_df$type <- factor("observed")
obs.f_star  <- obs.periodogram_df$peak_freq
obs.T_star  <- obs.periodogram_df$peak_period

ggplot(obs.per_df, aes(x = freq, y = power)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = obs.f_star, linetype = "dashed") +
  annotate("text", x = obs.f_star, y = max(obs.per_df$power),
           label = paste0("Peak: f ≈ ", signif(obs.f_star, 3),
                          "  (T ≈ ", signif(obs.T_star, 3), ")"),
           hjust = -0.1, vjust = 2, size = 3.75) +
  labs(title = paste0(sp, " — mean raw periodogram across sims"),
       x = "Frequency (cycles / time unit)", y = "Power") +
  coord_cartesian(xlim = c(0,1.5))+
  theme_minimal(base_size = 12)

#MAPPED (RESCALED) MODEL -------------------------------------------------------
#p is given above

map.pars <- pars
map.pars$beta <- pars$beta * p
map.pars$sigma_u <- pars$sigma_u * as.numeric(ps_for_sigma_u_mapping(N = ssize, pars, p))

map.traj_df <- simulate_many_grid_traj(
  nsims = nsims,
  k = k, alpha = alpha, burst.size = burst.size,
  Tmax = Tmax, y0 = y0, pars = map.pars, dt = dt,
  mc.cores = max(1L, parallel::detectCores() - 1L)
) 

map.some_traj <- map.traj_df[(map.traj_df$sim %in% seq(1,10,by=1)) & (map.traj_df$species == sp), ]
ggplot(map.some_traj, aes(x = t, y = y_actual, group = sim)) +
  geom_line(alpha = 0.15, linewidth = 0.3) +
  labs(title = paste0(sp, ": trajectories (", nsims, " sims)"),
       x = "time", y = "concentration") +
  theme_minimal(base_size = 12)

map.traj_df_edit <- map.traj_df
map.traj_df_edit$species <- "y_1"
map.traj_df_edit$sim <- rep(seq_len(3*nsims), each = 3001)

map.traj_df_edit <- fread(here::here("DATA", "3S_periodogram/p_0.1/map_traj_df.csv"))

map.series_from_traj_list <- extract_series_list_fast(map.traj_df_edit, species = sp)

#autocorrelation coefficient 
map.acf_df <- mean_acf_for_species(series_list = map.series_from_traj_list, 
                                   species = sp, dt = dt, max_lag_time = Tmax)

ggplot(map.acf_df, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  labs(title = paste0(sp, " — mean autocorrelation across sims"),
       x = "lag (time units)", y = "ACF") +
  theme_minimal(base_size = 12)

# periodogram (mean across simulations)
map.periodogram_df <- mean_periodogram_for_species_slow(
                              series_list = map.series_from_traj_list, 
                              species = sp, dt = dt,
                              include_0freq = T, demean = T,
                              normalize_area = T)

map.per_df  <- map.periodogram_df$periodogram_df
map.per_df$type <- factor("mapped")
map.f_star  <- map.periodogram_df$peak_freq
map.T_star  <- map.periodogram_df$peak_period

ggplot(map.per_df, aes(x = freq, y = power)) +
  geom_line(linewidth = 0.7) +
  geom_vline(xintercept = map.f_star, linetype = "dashed") +
  annotate("text", x = map.f_star, y = max(map.per_df$power),
           label = paste0("Peak: f ≈ ", signif(map.f_star, 3),
                          "  (T ≈ ", signif(map.T_star, 3), ")"),
           hjust = -0.1, vjust = 2, size = 3.75) +
  labs(title = paste0(sp, " — mean raw periodogram"),
       x = "Frequency (cycles / time unit)", y = "Power (energy)") +
  coord_cartesian(xlim = c(0,1.5))+
  theme_minimal(base_size = 12)

# SAVE DATA --------------------------------------------------------------------
# 
# write.csv(act.traj_df_edit,
#           here::here("DATA", "3S_periodogram/act_traj_df2.csv"), row.names = F)
# 
write.csv(act.periodogram_df,
          here::here("DATA", "3S_periodogram/act_periodogram_df.csv"), row.names = F) #+
# 
# write.csv(obs.traj_df_edit,
#           here::here("DATA", "3S_periodogram/p_0.1/obs_traj_df.csv"), row.names = F)
# 
write.csv(obs.periodogram_df,
          here::here("DATA", "3S_periodogram/p_0.1/obs_periodogram_df.csv"), row.names = F) #+++

# write.csv(map.traj_df_edit,
#           here::here("DATA", "3S_periodogram/p_0.1/map_traj_df.csv"), row.names = F)
# 
write.csv(map.periodogram_df,
          here::here("DATA", "3S_periodogram/p_0.1/map_periodogram_df.csv"), row.names = F) #+++

# FINAL PLOTS ------------------------------------------------------------------

levs <- c("actual", "observed", "mapped")

common.period_df <- rbind(act.per_df, obs.per_df, map.per_df)

common.peaks_df <- data.frame(
  type   = factor(c("actual", "observed", "mapped"), levels = levs),
  f_star = c(act.f_star,  obs.f_star,  map.f_star),
  T_star = c(act.T_star,  obs.T_star,  map.T_star),
  y_max  = c(max(act.per_df$power, na.rm = TRUE),
             max(obs.per_df$power, na.rm = TRUE),
             max(map.per_df$power, na.rm = TRUE))
)

common.peaks_df$label <- paste0("Peak: f \u2248 ", signif(common.peaks_df$f_star, 3),
                         "\nPeriod: T \u2248 ", signif(common.peaks_df$T_star, 3))

gg1<- ggplot(common.period_df, aes(x = freq, y = power, color = type)) +
  geom_vline(data = subset(common.peaks_df), 
             aes(xintercept = f_star), color = "black",
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  geom_text(data = subset(common.peaks_df),
            aes(x = f_star, y = y_max, label = label),
            #hjust = -0.1, vjust = -3.6, color = "black",
            size = 3.75, show.legend = FALSE) +
  geom_line(data = subset(common.period_df, type %in% c("actual", "mapped") ),
            linewidth = 1) +
  geom_point(data = subset(common.period_df, type == "observed"),
             size = 2) +
  geom_line(data = subset(common.period_df, type == "observed"),
             lty = "dashed") +
  labs(title = "",
       x = "Frequency (cycles / time unit), f", y = "Variance-normalized power spectrum density") +
  coord_cartesian(xlim = c(0, 0.75), 
                  ylim = c(0, max(common.period_df$power)), expand = T) +
  scale_color_manual(
    name   = NULL,
    aesthetics = c("colour", "fill"),
    values = c(observed = "#79CEDB", mapped = "#873DC2", actual = "#f78b59"),
    labels = c(observed = "Observed via kernel",
               mapped   = "Model with mapped parameters",
               actual   = "True model"),
  ) + 
  gg_base0; gg1

# DISTRIBUTIONS ----------------------------------------------------------------

act.traj_df_edit$type <- factor("actual")
names(act.traj_df_edit)[names(act.traj_df_edit) == 'y_actual'] <- 'concentration'

obs.traj_df_edit$type <- factor("observed")
names(obs.traj_df_edit)[names(obs.traj_df_edit) == 'x_observed'] <- 'concentration'

map.traj_df_edit$type <- factor("mapped")
names(map.traj_df_edit)[names(map.traj_df_edit) == 'y_actual'] <- 'concentration'

act.for_density <- act.traj_df_edit[act.traj_df_edit$t == 30, c("species", "concentration", "type")]
obs.for_density <- obs.traj_df_edit[obs.traj_df_edit$t == 30, c("species", "concentration", "type")]
map.for_density <- map.traj_df_edit[map.traj_df_edit$t == 30, c("species", "concentration", "type")]


gg2 <- ggplot() + 
  geom_histogram(
    data = obs.for_density,
    aes(x = concentration, y = after_stat(density), fill = type),
    bins = 80, position = "identity", alpha = 1, colour = NA
  ) +
  stat_density(           #compare to histogram
    data = map.for_density,
    aes(x = concentration, y = after_stat(density), colour = type),
    geom = "line", position = "identity", linewidth = 1
  ) +
  stat_density(           #compare to histogram
    data = act.for_density,
    aes(x = concentration, y = after_stat(density), colour = type),
    geom = "line", position = "identity", linewidth = 1
  ) + 
  labs(title = "",
       x = "Protein concentration, x", y = "Protein distribution") +
  coord_cartesian(expand = F) +
  scale_color_manual(
    name   = NULL,
    aesthetics = c("colour", "fill"),
    values = c(observed = "#79CEDB", mapped = "#873DC2", actual = "#f78b59"),
    labels = c(observed = "Observed via kernel",
               mapped   = "Model with mapped parameters",
               actual   = "True model"),
  ) + 
  gg_base0; gg2

ggarrange(gg1, gg2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)  
