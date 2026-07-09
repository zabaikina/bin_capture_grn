source(here::here("R", "simulation_functions.R"))

append_segment <- function(y0, t0, t1, dt, k, N){
  
  if (t1 <= t0) 
    return(list(times = numeric(0), states = matrix(numeric(0), ncol = N)))
  
  ts <- seq(t0, t1, by = dt) #grid #if t0-t1<dt then ts = t0
  ts <- c(ts, t1)
  
  S <- matrix(NA_real_, nrow = length(ts), ncol = N)
  for (r in seq_along(ts)) 
  {
    S[r, ] <- y0 * exp(-k * (ts[r] - t0))
  }
  
  list(times = ts, states = S)
}


simulate_once_traj <- function(
    k, alpha, burst.size, Tmax, 
    y0, pars = list(), dt = 0.01, seed = NULL)
  {
  
  if (!is.null(seed)) set.seed(seed)
  
  y <- y0; t <- 0 #loop variables
  
  N <- length(y)

  times <- c(0)
  states <- matrix(y, nrow = 1, ncol = N)
  
  while (t < Tmax) {
    
    a <- alpha(y, t, pars)
    tau_vec <- -log(runif(N)) / a
    i <- which.min(tau_vec)
    
    t_next <- t + tau_vec[i]
    
    #graphic part 
    t_segment <- min(t_next, Tmax)
    seg <- append_segment(y, t, t_segment, dt, k, N)
    if (length(seg$times)) {
      times  <- c(times,  seg$times)
      states <- rbind(states, seg$states)
    }
    
    if (t_next >= Tmax) break
    
    y <- decay.rate(y, k, t_next-t)
    t <- t_next
    y[i] <- y[i] + burst.size(i, y, t, pars)
    
  }

  colnames(states) <- paste0("y_", seq_len(N))
  traj_wide <- data.frame(t = times, states, check.names = FALSE)
  
  # tidy format for ggplot
  species_names <- paste0("y_", seq_len(N))
  traj_long <- data.frame(
    t       = rep(traj_wide$t, times = N),
    species = rep(species_names, each = nrow(traj_wide)),
    y_actual   = as.vector(as.matrix(traj_wide[ , -1, drop = FALSE])),
    check.names = FALSE
  )
  
  traj_long
  
}

