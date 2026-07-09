make_logspline_df <- function(df, kind_name,
                              xcol = "val", split_cols = c("p","t"),
                              lbound = 0, ubound = 400,
                              n = 300, expand_max = 1.05, nkn = 2) {
  
  d <- df[df$kind == kind_name, , drop = FALSE]
  if (nrow(d) == 0) return(data.frame())
  
  groups <- unique(d[, split_cols, drop = FALSE])
  out <- vector("list", nrow(groups))
  eps <- 1e-8  # NEW: small nudge
  
  for (i in seq_len(nrow(groups))) {
    p_i <- groups[i, split_cols[1]][[1]]
    t_i <- groups[i, split_cols[2]][[1]]
    
    sel <- d[[split_cols[1]]] == p_i & d[[split_cols[2]]] == t_i
    x   <- d[sel, xcol]
    
    if (!length(x)) { out[[i]] <- NULL; next }
    
    x <- x[is.finite(x)]
    x <- x[x >= lbound]
    if (!length(x)) { out[[i]] <- NULL; next }
    
    x[x <= lbound] <- lbound + eps               # NEW: nudge off lower bound
    x_max  <- max(x)
    
    ub_fit <- min(ubound, x_max * 1.005)         # ← a hair more headroom (0.5%)  [CHANGED]
#    nkn    <- 2                                   # ← fewer knots (stronger smoothing) [CHANGED]
    
    fit <- logspline(x, lbound = lbound, ubound = ub_fit, nknots = nkn)  # NEW
    
    x_high <- min(ubound, x_max * expand_max)
    x_grid <- seq(lbound + eps, x_high, length.out = n)   # NEW: avoid exact lbound
    
    dens <- dlogspline(x_grid, fit)
    
    out[[i]] <- data.frame(
      x = x_grid, y = dens,
      kind = kind_name, p = p_i, t = t_i
    )
  }
  do.call(rbind, out)
}
