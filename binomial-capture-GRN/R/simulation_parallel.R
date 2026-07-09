

simulate_many_mcl <- function(N, n_sims, k, alpha, burst.size, Tmax, y0_mean, y0_sd,
                              pars = list(), seed = NULL, rec_times,
                              mc.cores = max(1L, parallel::detectCores() - 1L)) {
  # create per-run seeds no matter what
  if (!is.null(seed)) set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_sims)
  
  pieces <- parallel::mclapply(
    X = seq_len(n_sims),
    FUN = function(i)
      simulate_once_mcl(N, k, alpha, burst.size, Tmax, y0_mean, y0_sd, pars,
                        seed = seeds[i], rec_times = rec_times),
    mc.cores = mc.cores,
    mc.set.seed = FALSE  # we control seeding ourselves
  )
  
  runs <- data.table::rbindlist(pieces, use.names = TRUE)
  runs[, t := factor(t, levels = sort(unique(t)))]
  runs[]
}


simulate_once_mcl <- function(
    N,
    k,                 # degradation rates (scalar or length-N), first-order
    alpha,             # function(y, t, pars) -> length-N rates
    burst.size,        # function(i, y, t, pars) -> jump size for species i
    Tmax,              # time of observation 
    y0_mean, y0_sd,    # initial conditions ~ N(mean, sd), abs() applied
    pars = list(),
    seed = NULL, 
    rec_times = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # init
  y <- abs(stats::rnorm(N, mean = y0_mean, sd = y0_sd))
  t <- 0.0
  
  # record times (append Tmax)
  rec_times <- c(rec_times, Tmax)
  R <- length(rec_times)
  
  # prealloc
  rec_y <- matrix(NA_real_, nrow = R, ncol = N)
  s <- 1L
  species <- paste0("y_", seq_len(N))
  
  # vectorized first-order decay snapshots from current (y, t)
  decay_block <- function(y_now, k_vec, dt_vec) {
    ef <- exp(-outer(dt_vec, k_vec))          # [len(dt) x N]
    ef * rep(y_now, each = length(dt_vec))    # broadcast y across rows
  }
  
  while (t < Tmax) {
    a <- alpha(y, t, pars)
    tau_vec <- stats::rexp(N, rate = a)
    i <- which.min(tau_vec)
    tau <- tau_vec[i]
    
    # fill all snapshots up to t + tau from current state (no mutation)
    upto <- findInterval(t + tau, rec_times)  # number of rec_times <= t+tau
    if (upto >= s) {
      dt <- rec_times[s:upto] - t
      rec_y[s:upto, ] <- decay_block(y, k, dt)
      s <- upto + 1L
    }
    
    # if next event happens after Tmax, finish snapshots and exit
    if (t + tau >= Tmax) {
      if (s <= R) {
        dt <- rec_times[s:R] - t
        rec_y[s:R, ] <- decay_block(y, k, dt)
      }
      break
    }
    
    # advance to event (in-place decay) and apply jump
    y[] <- y * exp(-k * tau)     # k can be scalar or length-N
    t   <- t + tau
    y[i] <- y[i] + burst.size(i, y, t, pars)
  }
  
  data.frame(
    t         = rep(rec_times, each = N),
    species   = factor(rep(species, times = R), levels = species),
    y_actual  = as.vector(t(rec_y)),
    check.names = FALSE
  )
}
