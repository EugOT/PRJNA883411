GetDoubletRate <- function(ncells) {
  ref <-
    tibble(
      multiplate_rate = c(
        .004, .008, .016, .023, .031, .039, .046, .054,  .061,  .069,  .076),
      cells_loaded = c(
        800,  1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
      cells_recovered = c(
        500,  1000, 2000, 3000, 4000, 5000, 6000, 7000,  8000,  9000,  10000))
  cls <- ifelse(ncells > 700, round(ncells / 1000) * 1000, 500)
  rate <-
    ref %>%
    dplyr::filter(cells_recovered <= cls) %>%
    dplyr::top_n(n = 1, wt = multiplate_rate) %>%
    dplyr::mutate(multiplate_rate = multiplate_rate / cells_recovered * ncells) %>%
    pluck("multiplate_rate")

  return(rate)
}

AddExpectedDoublets <- function(dat) {
  doublets <-
    dat$NTotalCells %>%
    map(GetDoubletRate) %>%
    simplify()
  dd <- tibble(NExpectedDoubletRate = doublets)
  df <- bind_cols(
    dat, dd
  )
  write_tsv(
    x = df,
    file = "samples.tsv")
  return(df)
}

library(vroom)
library(tidyverse)
library(magrittr)


read_tsv("samples.tsv") %>% AddExpectedDoublets()
