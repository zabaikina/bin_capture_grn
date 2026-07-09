# Simulate many trajectories on a common grid and return LONG trajectories
# Output columns: sim, t, species, y_actual
simulate_many_grid_traj <- function(
    nsims, k, alpha, burst.size, Tmax, y0, pars = list(),
    dt = 0.02, mc.cores = max(1L, parallel::detectCores() - 1L)
) {
  N <- length(y0)
  obs_times <- unique(c(seq(0, Tmax, by = dt), Tmax))
  nT <- length(obs_times)
  
  decay_to <- function(y, k, delta) {
    y * exp(-k * delta)
  }
  
  sim_one_to_grid <- function() {
    Y <- matrix(NA_real_, nrow = nT, ncol = N)
    Y[1, ] <- y0
    
    t <- 0.0
    y <- y0
    idx <- 2L
    
    while (t < Tmax) {
      a <- alpha(y, t, pars)
      
      # draw next reaction times; nonpositive/invalid rates -> no event (Inf)
      tau_vec <- rep(Inf, N)
      for (j in seq_len(N)) {
        if (is.finite(a[j]) && a[j] > 0) {
          tau_vec[j] <- rexp(1, rate = a[j])
        }
      }
      
      i <- which.min(tau_vec)
      t_next <- t + tau_vec[i]
      t_seg  <- if (is.finite(t_next)) min(t_next, Tmax) else Tmax
      
      # fill grid points between current t and t_seg via pure decay
      while (idx <= nT && obs_times[idx] <= t_seg) {
        delta <- obs_times[idx] - t
        Y[idx, ] <- decay_to(y, k, delta)
        idx <- idx + 1L
      }
      
      # stop if no finite event or reached Tmax
      if (!is.finite(tau_vec[i]) || t_next >= Tmax) break
      
      # jump: advance to t_next, then add burst to component i
      y <- decay_to(y, k, t_next - t)
      t <- t_next
      y[i] <- y[i] + burst.size(i, y, t, pars)
    }
    
    # if time remains to Tmax, continue pure decay on the grid
    if (idx <= nT) {
      while (idx <= nT) {
        delta <- obs_times[idx] - t
        Y[idx, ] <- decay_to(y, k, delta)
        idx <- idx + 1L
      }
    }
    
    Y
  }

  Y_list <- parallel::mclapply(
    X = seq_len(nsims),
    FUN = function(s) sim_one_to_grid(),
    mc.cores = mc.cores
  )
  
  species_names <- paste0("y_", seq_len(N))
  df_list <- vector("list", length = nsims)
  
  for (s in seq_len(nsims)) {
    M <- Y_list[[s]]  # nT x N
    df_list[[s]] <- data.frame(
      sim      = rep.int(s, nT * N),
      t        = rep(obs_times, times = N),
      species  = rep(species_names, each = nT),
      y_actual = as.vector(M),
      check.names = FALSE
    )
  }
  
  traj_df <- do.call(rbind, df_list)
  traj_df
}


simulate_many_grid_stats <- function(
    nsims, k, alpha, burst.size, Tmax, y0, pars = list(),
    dt = 0.02, mc.cores = max(1L, detectCores() - 1L)
) {
  N <- length(y0)
  obs_times <- unique(c(seq(0, Tmax, by = dt), Tmax))
  nT <- length(obs_times)
  
  decay_to <- function(y, k, delta) {
    y * exp(-k * delta)
  }
  
  sim_one_to_grid <- function() {
    Y <- matrix(NA_real_, nrow = nT, ncol = N)
    Y[1, ] <- y0
    
    t <- 0.0
    y <- y0
    idx <- 2L  # next grid row to fill
    
    while (t < Tmax) {
      a <- alpha(y, t, pars)
      
      tau_vec <- rep(Inf, N)
      for (j in seq_len(N)) {
        if (is.finite(a[j]) && a[j] > 0) tau_vec[j] <- rexp(1, rate = a[j])
      }
      i <- which.min(tau_vec)
      t_next <- t + tau_vec[i]
      t_seg  <- min(t_next, Tmax)
      
      # fill grid points
      while (idx <= nT && obs_times[idx] <= t_seg) {
        delta <- obs_times[idx] - t
        Y[idx, ] <- decay_to(y, k, delta)
        idx <- idx + 1L
      }
      
      # stop if no finite event or we've reached Tmax
      if (!is.finite(tau_vec[i]) || t_next >= Tmax) break
      
      
      y <- decay_to(y, k, t_next - t)
      t <- t_next
      y[i] <- y[i] + burst.size(i, y, t, pars)
    }
    
    # if anything remains to Tmax (no more jumps), keep decaying on the grid
    if (idx <= nT) {
      while (idx <= nT) {
        delta <- obs_times[idx] - t
        Y[idx, ] <- decay_to(y, k, delta)
        idx <- idx + 1L
      }
    }
    
    Y
  }
  
  
  Y_list <- mclapply(seq_len(nsims), function(s) sim_one_to_grid(), mc.cores = mc.cores)
  
  # --- aggregate mean and sd without pipes ---
  # sum and sum of squares over sims
  sum_mat   <- Reduce(`+`, Y_list)
  sq_list   <- lapply(Y_list, function(M) M * M)
  sumsq_mat <- Reduce(`+`, sq_list)
  
  mean_mat <- sum_mat / nsims
  var_mat  <- pmax(sumsq_mat / nsims - mean_mat * mean_mat, 0)
  sd_mat   <- sqrt(var_mat)
  
  species_names <- paste0("y_", seq_len(N))
  out_list <- vector("list", length = N)
  for (j in seq_len(N)) {
    out_list[[j]] <- data.frame(
      t       = obs_times,
      species = species_names[j],
      mean    = mean_mat[, j],
      sd      = sd_mat[, j],
      check.names = FALSE
    )
  }
  stats_df <- do.call(rbind, out_list)
  stats_df
}
