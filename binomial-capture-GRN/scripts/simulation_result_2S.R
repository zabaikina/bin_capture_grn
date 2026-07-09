source(here::here("R", "0preambule.R"))

# 2D PARAMS --------------------------------------------------------------------
B <- 3
L <- 40 
H <- 200
K1 <- (1-0.1)*H*L
Rb <- (L*H-K1)/(B*(L+H))
Ru <-( K1*(L^2+L*H+H^2) + (H*L)^2 ) / (B*K1*(H+L))

runtime <- 40
rec_times <- c(0.05*runtime, 0.1125*runtime); rec_times
y0_mean <- rep((L/2+H/2), 2); y0_sd <- c(10,10)
n_sims <- 10^5
ssize <<- 2

k    <- c(1, 1)      # decay rates
pars <- list(
  rho_u   = c(Ru, Ru),
  rho_b   = c(Rb, Rb),
  sigma_u = c(K1, K1),
  sigma_b = c(1, 1),
  beta = c(B, B),
  n = matrix( c(0 , 2,
                2, 0), nrow = 2, byrow = T)
)


p_vals   <- matrix(c( 1, 1, 
                      0.75, 0.75,
                      0.5, 0.5,
                      0.25, 0.25,
                      0.1, 0.1, 
                      0.05, 0.05, 
                      0.01, 0.01), ncol = 2, byrow = T)     
pN <- nrow(p_vals) #number of pairs of p
p_labels <- paste0("p=", format(p_vals[,1])) 

# COMPARISON FOR P_i -----------------------------------------------------------


tic()
act_all_2S <- simulate_many_mcl(
  N = ssize,
  n_sims = n_sims,
  k = k, alpha = alpha, burst.size = burst.size,
  Tmax = runtime, 
  y0_mean = y0_mean, y0_sd = y0_sd,
  pars = pars, seed = NULL,
  rec_times = rec_times
); toc()

ggplot(act_all_2S, aes(x = y_actual)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 position = "identity", alpha = 0.75, color = NA, binwidth = 3)+
  stat_bin(
    aes(y = after_stat(density), color = species),
    bins = 50, geom = "area", position = "identity",
    pad = TRUE, alpha = 0.15, linewidth = 0.7
  ) +
  facet_wrap(~ t, nrow = 1) +
  labs(x = "concentration", y = "density",
       fill = "", color = "") +
  theme_minimal()

#observed conc X through norm kernel
obs_list_2S <- vector("list", length = pN)
for (i in 1:pN) {
  tmp <- observe_x_normal(act_all_2S, p = c(y_1 = p_vals[i,1], y_2 = p_vals[i,2]))
  tmp$p <- p_labels[i]
  obs_list_2S[[i]] <- tmp
}
obs_all_2S <- do.call(rbind, obs_list_2S)
obs_all_2S$p <- factor(obs_all_2S$p, levels = p_labels)

#mapped X (original model w rescaled params)
map_list_2S <- vector("list", length = pN)
for (i in 1:pN) {
  pars_map <- pars
  pars_map$beta <- pars$beta * p_vals[i,]
  pars_map$sigma_u <- pars$sigma_u * 
    as.numeric(ps_for_sigma_u_mapping(N = ssize, pars, p_vals[i,]))
  
  tmp <- simulate_many_mcl(N = ssize, n_sims = n_sims,
                           k = k, alpha = alpha, burst.size = burst.size,
                           Tmax = runtime, 
                           y0_mean = y0_mean*p_vals[i,], y0_sd = y0_sd*p_vals[i,],
                           pars = pars_map, seed = NULL, rec_times = rec_times)
  names(tmp) <- c("t","species", "x_mapped")
  tmp$p <- p_labels[i]
  
  map_list_2S[[i]] <- tmp
}
map_all_2S <- do.call(rbind, map_list_2S)
map_all_2S$p <- factor(map_all_2S$p, levels = p_labels)

gg_df_2S <- rbind(
  data.frame(t = map_all_2S$t, species = map_all_2S$species, 
             kind = "mapped", val = map_all_2S$x_mapped, p = map_all_2S$p),
  data.frame(t = obs_all_2S$t, species = obs_all_2S$species, 
             kind = "observed", val = obs_all_2S$x_observed, p = obs_all_2S$p)
)
gg_df_2S$kind <- factor(gg_df_2S$kind, levels = c("observed","mapped"))

means_2S <- ggplot_means_df(obs_all_2S, map_all_2S, include = c("observed","mapped"), 
                         vjust = 5, p = p_vals)

#write.csv(gg_df_2S, here::here("DATA", "gg_df_2S_L40_H200_B3_n2_case3.csv"), row.names = F)
#gg_df_2S <- read.csv(here::here("DATA", "gg_df_2S_L40_H200_B3_n2_case3.csv"))


gg_df_2S <- subset(gg_df_2S, (val >= 0))
ls_df_2S <- make_logspline_df(gg_df_2S, kind_name = "mapped",
                           lbound = 0, ubound = max(gg_df_2S$val),
                           n = 400, nkn = 2)


gg_df_2S_chop <- subset(gg_df_2S, p %in% c("p=1.00", "p=0.50", "p=0.25"))
ls_df_2S_chop <- make_logspline_df(gg_df_2S_chop, kind_name = "mapped",
                                   lbound = 0, ubound = max(gg_df_2S_chop$val),
                                   n = 400, nkn = 2)


# # Simulate and draw one trajectory ---------------------------------------------
# 
# sample_traj <- simulate_once_traj(k, alpha, burst.size,
#                                   Tmax = 100, y0 = y0_mean, pars, dt = 0.02)
# 
# ggplot(sample_traj, aes(t, y_actual, color = species)) +
#   geom_line(linewidth = 0.6) +
#   labs(x = "time", y = "concentration", color = "species") +
#   coord_cartesian(ylim = c(0, max(sample_traj$y_actual)))+
#   theme_minimal()