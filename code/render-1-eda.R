#!/usr/bin/env Rscript
# render-1-eda.R

reseed <- 42
set.seed(seed = reseed)
library(here)
library(workflowr)
analysis_dir <- here("analysis")
workflowr::wflow_publish(
  seed = reseed,
  verbose = TRUE
)
