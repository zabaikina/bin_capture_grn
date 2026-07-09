# Install all packages required to run the code in this repository.
pkgs <- c(
  "here", "tictoc", "stringr",
  "ggplot2", "cowplot", "ggpubr", "wesanderson", "scales", "forcats",
  "parallel", "future", "future.apply", "data.table",
  "pracma", "RConics", "logspline", "deSolve"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
