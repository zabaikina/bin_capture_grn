

# AUXILIARY FUNCTIONS ----------------------------------------------------------

ps_for_sigma_u_mapping <- function(N, pars, p_vec) {

  base <- matrix(p_vec, nrow = N, ncol = N, byrow = TRUE)

  p_to_n <- base^pars$n

  p_for_sig <- apply(p_to_n, 1L, prod)
  p_for_sig
}

sd_for_normIC_mapping <- function(y0_mean, y0_sd, p_vec){
  #new_sd <- sqrt( (p_vec^2)*(y0_sd^2) + (1-p_vec)*p_vec*y0_mean )
  new_sd <- y0_sd*p_vec
  new_sd
}


y_mean_meanfield <- function(pars, k) {
  
  cubic_coef <- c(1, 
                 - pars$beta * pars$rho_b/k, 
                 pars$sigma_u/pars$sigma_b, 
                 -pars$beta * pars$sigma_u/pars$sigma_b * pars$rho_u/k)
  roots <- cubic(cubic_coef)
  tolerate <- 1e-15
  re_roots <- Re(roots[abs(Im(roots)) < tolerate])
  valid_root <- re_roots[re_roots > 0] 
      #NOTE: if we have 3 real roots then we may need better criteria
  re_roots
}

# means helper ---------------------------------------------------------------
ggplot_means_long <- function(
    df_long,
    groups     = c("t", "species", "p"),
    value_col  = "val",                  # numeric column to average
    value_name = "mean", 
    include    = c("observed", "mapped", "actual"),
    vjust      = 1,
    vjust_step = 1.5,
    p          = NULL
){
  
  if (!is.null(p) && length(p) == 1) df_long$p <- factor(p)
  
  df_long <- subset(df_long, kind %in% include)
  
  df_long$kind <- factor(df_long$kind, levels = c("actual","observed","mapped"))
  
  groups <- intersect(groups, names(df_long))
  
  by_vars <- c(groups, "kind")
  fml <- reformulate(by_vars, response = value_col)
  out <- aggregate(fml, data = df_long, FUN = mean, na.rm = TRUE)
  
  names(out)[names(out) == value_col] <- value_name
  
  out$vjust <- vjust + (as.integer(out$kind) - 1) * vjust_step
  out
}



#-------------------------------------------------------------------------------

# BASIC FUNCTIONS --------------------------------------------------------------

decay.rate <- function(y, k, dt){ 
  y * exp(-k *dt)
} #returns vector

burst.size <- function(i, y, t, pars){
  rexp(1, rate = 1 / pars$beta[i])
} #returns burst size for given species i

alpha <- function(y, t, pars){ #this is in general form just in case 
  
  N <- length(y)
  
  val <- function(x) if (is.function(x)) x(t) else x
  rho_u   <- val(pars$rho_u)    
  rho_b   <- val(pars$rho_b)  
  sigma_u <- val(pars$sigma_u)    
  sigma_b <- val(pars$sigma_b)
  K       <- sigma_u/sigma_b    
  nmat    <- pars$n            # matrix of n_{ij}
  
 
  prod_term <- numeric(N)      # this is product of y_j^{n_{ij}} in denom
  for (i in 1:N) prod_term[i] <- prod(y ^ nmat[i, ])
  
  (rho_u * K + rho_b * prod_term) / (K + prod_term)
}

#-------------------------------------------------------------------------------

# SIMULATE AND TRACK Y AT MULTIPLE TIME POINTS ---------------------------------

simulate_once <- function(
    N,
    k,                 # degradation rates
    alpha,             # functions of freq at t given vec(y)
    burst.size,        # function of burst size
    Tmax,              # time of observation 
    y0_mean, y0_sd,                # IC 
    pars = list(),
    seed = NULL, 
    rec_times = NULL)
  {
  
  #if (!is.null(seed) || !is.na(seed)) set.seed(seed) #reproduce the same results if needed
  
  #initiate loop variables
  y0 <- abs(rnorm(N, mean = y0_mean, sd = y0_sd))
  y <- y0
  t <- 0
  
  #intermediate records
  
  rec_times <- c(rec_times, Tmax)
  
  rec_y <- matrix(NA_real_, 
                  nrow = length(rec_times), 
                  ncol = N,
                  dimnames = list(NULL, paste0("y_", seq_len(N)))) 
  s <- 1L
  
  #run
  while (t < Tmax) {

    a <- alpha(y, t, pars)      # reaction rates CURRENT state

    tau_vec <- -log(runif(N))/a   #vec of waiting times
    i   <- which.min(tau_vec)     #next reaction number
    tau <- tau_vec[i]             #next reaction wt
    
    while (s <= length(rec_times) && (rec_times[s] <= (t+tau))) {
      rec_y[s, ] <- decay.rate(y, k, rec_times[s] - t)   # snapshot from current (y,t) → no state mutation
      s <- s + 1L
    }
    
    if (t + tau >= Tmax) break
    
    #update loop variables
    y <- decay.rate(y, k, tau)
    t <- t + tau
    y[i] <- y[i] + burst.size(i, y, t, pars)
    
  }
  
  species <- paste0("y_", seq_len(N))
  
  data.frame(
    t       = rep(rec_times, each = N),
    species = factor(rep(species, times = length(rec_times)), levels = species),
    y_actual   = as.vector(t(rec_y)),   # row-wise flatten (times x species)
    check.names = FALSE )
  
}

simulate_many <- function(
    N, #system size
    n_sims,  #sim number 
    k, 
    alpha, 
    burst.size,
    Tmax, #runtime
    y0_mean, y0_sd, #mean for normal distr at t=0
    pars = list(), 
    seed = NULL,
    rec_times){
  
  if (!is.null(seed)) 
  {
    set.seed(seed) #fix random generator 
    seeds <- sample.int(1e9, n_sims) #reproduce the same random sequences for each single simulation
  } else {
    seeds <- rep(NULL, times = n_sims)
  }
  
  runs <- data.frame(
    t       = numeric(),
    species = factor(),
    y_actual   = numeric())
  
  for (i in seq_len(n_sims)) {
    tmp <- simulate_once(N, k, alpha, burst.size, 
                         Tmax, y0_mean, y0_sd, pars, seed = seeds[i], rec_times)
    runs <- rbind(runs, tmp)
  }
  
  runs <- runs[which(runs$y_actual > 0),] #remove nonpositive y
  
  runs$t <- factor(runs$t, levels = sort(unique(runs$t)))
  #times as factors for grouping
  
  runs 
}

observe_x_normal_1D <- function(df, p, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  y_actual <- df$y_actual
  
  if (p == 0) {
    x <- rep(0, times = length(y_actual))
  } else if (p == 1) {
    x <- y_actual
  } else {
    mu_y <- p * y_actual
    sd_y <- sqrt(p * (1 - p) * y_actual)
    x <- rnorm(length(y_actual), mean = mu_y, sd = sd_y)
  }
  
  df$x_observed <- x
  df
}

observe_x_normal <- function(df, p, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  y_row <- df$y_actual
  p_row <- p[as.character(df$species)]
  
  mu  <- p_row * y_row
  sdv <- sqrt(p_row * (1 - p_row) * y_row)
  
  x <- rnorm(length(y_row), mean = mu, sd = sdv)
  
  i0 <- (p_row == 0)
  if (any(i0)) x[i0] <- 0
  i1 <- (p_row == 1)
  if (any(i1)) x[i1] <- y_row[i1]
  
  df$x_observed <- x
  df
}

#-------------------------------------------------------------------------------

# SIMULATE AND TRACK Y AT Tmax ONLY --------------------------------------------

simulate_end <- function(
    k,                 # degradation rates
    alpha,             # functions of freq at t given vec(y)
    burst.size,        # function of burst size
    Tmax,              # time of observation 
    y0,                # IC 
    pars = list(),
    seed = NULL)
{
  
  if (!is.null(seed)) set.seed(seed) #reproduce the same results if needed
  
  #initiate loop variables
  y <- as.numeric(y0)
  t <- 0
  N <- length(y)

  #run
  while (t < Tmax) {
    
    a <- alpha(y, t, pars)      # reaction rates CURRENT state
    
    tau_vec <- -log(runif(N))/a   #vec of waiting times
    i   <- which.min(tau_vec)     #next reaction number
    tau <- tau_vec[i]             #next reaction wt
    
    if (t + tau >= Tmax) {             # fin
      y <- decay.rate(y, k, Tmax - t)
      break
    }
    
    #update loop variables
    y <- decay.rate(y, k, tau)
    t <- t + tau
    y[i] <- y[i] + burst.size(i, y, t, pars)
    
  }
  y
}

simulate_many_ends <- function(
    n_sims, 
    k, 
    alpha, 
    burst.size,
    Tmax, 
    y0, 
    pars = list(), 
    seed = NULL
){
  
  if (!is.null(seed)) 
  {
    set.seed(seed) #fix random generator 
    seeds <- sample.int(1e9, n_sims) #reproduce the same random sequences for each single simulation
  } else {
    seeds <- rep(NULL, times = n_sims)
  }
  
  N <- length(y0)
  
  runs <- matrix(NA_real_, nrow = n_sims, ncol = N)
  
  for (i in seq_len(n_sims)) {
    runs[i, ] <- simulate_end(k, alpha, burst.size, Tmax, y0, pars, seed = seeds[i])
  }
  
  df <- as.data.frame(runs, 
                      row.names = NULL)
  names(df) <- paste0("y_", seq_len(N))

  df
}






#-------------------------------------------------------------------------------
