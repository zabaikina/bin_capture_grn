#this is a helper that extracts data from "*_traj_df.csv*

base_dir <- here::here("DATA", "3S_periodogram")
p_values <- c("p_0.5", "p_0.25", "p_0.1")

n_threads <- max(1L, parallel::detectCores() - 1L)
T_target  <- 30

# tiny helper: read -> keep t==30 -> standardize -> save
simple_process <- function(in_path, out_path, conc_col, type_label) {
  DT <- fread(in_path, select = c("species", "t", conc_col), nThread = n_threads)
  setnames(DT, conc_col, "concentration")
  DT <- DT[t == T_target, .(species, concentration)]
  DT[, type := factor(type_label, levels = c("actual","observed","mapped"))]
  fwrite(DT, out_path)
}

# actual (in base dir)
simple_process(
  file.path(base_dir, "act_traj_df.csv"),
  file.path(base_dir, "act_for_density.csv"),
  conc_col   = "y_actual",
  type_label = "actual"
)

# mapped + observed (in each p_* subfolder)
for (p in p_values) {
  dir_p <- file.path(base_dir, p)
  
  simple_process(
    file.path(dir_p, "map_traj_df.csv"),
    file.path(dir_p, "map_for_density.csv"),
    conc_col   = "y_actual",
    type_label = "mapped"
  )
  
  simple_process(
    file.path(dir_p, "obs_traj_df.csv"),
    file.path(dir_p, "obs_for_density.csv"),
    conc_col   = "x_observed",
    type_label = "observed"
  )
}

gc()
