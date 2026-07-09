extract_series_list <- function(traj_df, species) {
  sim_ids <- sort(unique(traj_df$sim))
  L <- vector("list", length = length(sim_ids))
  for (i in seq_along(sim_ids)) {
    s <- sim_ids[i]
    df_s <- traj_df[traj_df$sim == s & traj_df$species %in% species, c("t","y_actual")]
            #sub dataframe for given simulation
    df_s <- df_s[order(df_s$t), ]
            #insurance that the data is ordered by time 
    L[[i]] <- df_s$y_actual #write time series to list
  }
  L
}

extract_series_list_fast <- function(traj_df, species, observed = F) {
  sims <- sort(unique(traj_df$sim))
  o    <- order(traj_df$sim, traj_df$t)          # one global sort by sim,t
  if (observed) {
    y    <- traj_df$x_observed[o]
  } else {
    y    <- traj_df$y_actual[o]
  }
  s    <- traj_df$sim[o]
  sp   <- traj_df$species[o]
  
  # filter species, then split by sim; keep empty sims and drop names
  unname(split(y[sp == species],
               factor(s[sp == species], levels = sims),
               drop = FALSE))
}


mean_acf_for_species <- function(series_list, species, dt, max_lag_time,
                                 mc.cores = max(1L, parallel::detectCores() - 1L)) {
  
  max_lag_N <- as.integer(max(1L, floor(max_lag_time / dt)))
  
  acf_list <- parallel::mclapply(
    X = series_list,
    FUN = function(x) {
      x <- x - mean(x)
      stats::acf(x, lag.max = max_lag_N, plot = F)$acf[, 1, 1]
    },
    mc.cores = mc.cores
  )
  
  acf_mat <- do.call(cbind, acf_list)  # (max_lag+1) x nsims
  lag_vec <- 0:max_lag_N
  
  acf_mean <- rowMeans(acf_mat)
  data.frame(lag = lag_vec * dt, acf = acf_mean, check.names = FALSE)
}

mean_periodogram_for_species <- function(series_list, species, dt,
                                         include_0freq = TRUE,
                                         demean = TRUE,
                                         normalize_area = TRUE,
                                         mc.cores = max(1L, parallel::detectCores() - 1L)) {
  
  n    <- length(series_list[[1]])
  kmax <- floor(n / 2)
  df   <- 1 / (n * dt)  # frequency step on your grid
  
  per_list <- parallel::mclapply(
    X = series_list,
    FUN = function(x) {
      if (demean) x <- x - mean(x)
      X <- fft(x)
      start_idx <- if (include_0freq) 1L else 2L
      idx <- start_idx:(kmax + 1L)
      P <- (dt / n) * (Mod(X[idx])^2)   
      
      lo <- start_idx + 1L
      hi <- (kmax + 1L) - if (n %% 2L == 0L) 1L else 0L
      if (hi >= lo) P[lo:hi] <- 2 * P[lo:hi]
      
      if (normalize_area) {
        B <<- 1
        A <- sum(P, na.rm = TRUE) * df   # area under this spectrum
        if (A > 0) P <- P / A            # unit-area (shape-only)
      }
      list(power = P)
    },
    mc.cores = mc.cores
  )
  
  f0   <- if (include_0freq) 0L else 1L
  freq <- (f0:kmax) / (n * dt)
  
  Pmat  <- sapply(per_list, `[[`, "power")
  Pmean <- rowMeans(Pmat)
  
  # peak search (still skip the first bin to avoid near-DC)
  skip <- 1L
  pk_idx <- which.max(Pmean[(skip + 1L):length(Pmean)]) + skip
  
  list(
    periodogram_df = data.frame(freq = freq, power = Pmean, check.names = FALSE),
    peak_freq   = freq[pk_idx],
    peak_period = 1 / freq[pk_idx]
  )
}




mean_periodogram_for_species_slow <- function(series_list, species, dt,
                                              include_0freq = TRUE,
                                              demean = TRUE,
                                              normalize_area = TRUE) {
  N    <- length(series_list[[1]])
  kmax <- floor(N / 2)  #N/2 Nyquist
  df   <- 1 / (N * dt)  # frequency bin 1/T
  
  start_k <- if (include_0freq) 0L else 1L  # k index 
  idx     <- (start_k:kmax) + 1L            # the same but for R indexes
  
  f0   <- if (include_0freq) 0L else 1L
  freq <- (f0:kmax) / (N * dt)    #frequency grid
  
  
  per_list <- lapply(series_list, function(x) {
   
    if (demean) x <- x - mean(x)
    
    X <- fft(x)  
    
    I <- (dt/N) * (Mod(X[idx])^2)   #  dt^2/T * |X_k|^2 that comes from continuous FT 
    
    # one sided periodogram
    L <- length(I)
    if (L > 0L) {
      if (include_0freq) {
        # P indices: 1=DC, last=Nyquist if n even
        if (N %% 2L == 0L) {         # even n: DC .. Nyquist
          if (L >= 3L) I[2:(L-1)] <- 2 * I[2:(L-1)]
        } else {                     # odd n: DC .. highest positive (no Nyquist)
          if (L >= 2L) I[2:L] <- 2 * I[2:L]
        }
      } else {
        # P indices: 1=first positive k=1, last=Nyquist if n even
        if (N %% 2L == 0L) {         # even n: k=1 .. Nyquist
          if (L >= 2L) I[1:(L-1)] <- 2 * I[1:(L-1)]
        } else {                     # odd n: k=1 .. highest positive
          if (L >= 1L) I[1:L] <- 2 * I[1:L]
        }
      }
    }
    
    if (normalize_area) {
      #pop variance 
      pop_var <- mean(x^2)
      I <- I/pop_var
    }
    
    list(power = I)
  })
  
  Pmat  <- do.call(cbind, lapply(per_list, `[[`, "power"))
  Pmean <- rowMeans(Pmat)
  
  # peak (skip first bin heuristically to avoid near-DC)
  skip <- 1L
  if (length(Pmean) > skip + 1L) {
    pk_idx <- which.max(Pmean[(skip + 1L):length(Pmean)]) + skip
  } else {
    pk_idx <- which.max(Pmean)
  }
  
  list(
    periodogram_df = data.frame(freq = freq, power = Pmean, check.names = FALSE),
    peak_freq   = freq[pk_idx],
    peak_period = 1 / freq[pk_idx]
  )
}






