#!/usr/bin/env Rscript
# render-1-eda.R

reseed <- 42
set.seed(seed = reseed)
library(here)
library(workflowr)
analysis_dir <- here("analysis")
workflowr::wflow_publish(
  files = here(analysis_dir, "*.Rmd"),
  all = TRUE,
  seed = reseed,
  verbose = TRUE
)
