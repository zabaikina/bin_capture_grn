# LIBRARIES --------------------------------------------------------------------
library(here)
library(tictoc)
library(stringr)

library(ggplot2)
library(cowplot)
library(ggpubr)
library(wesanderson)

library(parallel)
library(future)
library(future.apply)
library(data.table)

library(pracma)
library(RConics)
library(logspline)
library(forcats)
library(scales)
library(deSolve)

# FILES ------------------------------------------------------------------------------
source(here::here("R", "simulation_functions.R"))
source(here::here("R", "simulation_parallel.R"))
source(here::here("R", "simulation_trajectories_correct.R"))
source(here::here("R", "simulation_trajectories_on_grid.R"))
source(here::here("R", "ggplot_line_logspline.R"))
source(here::here("R", "ggplot_ggbase0.R"))




