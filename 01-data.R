### Downloading and processing the data for this tutorial.
### This can be downloaded from the following URL: https://www.stata-press.com/data/fpsaus.html
### We also include the dataset in Stata format in the data-raw/ folder for convenience.

# Packages
library(haven)
library(tidyverse)

# Load data
rott2 <- read_dta(file = "data-raw/rott2.dta")

# Remove all Stata attributes/labels
rott2 <- rott2 |>
  zap_formats() |>
  zap_label() |>
  zap_labels() |>
  zap_missing()

# Generate follow-up time time as the time to death or time to relapse (depending on which event occurred first)
# Note: There are 43 subjects who have died without relapse, but their time of death is greater than the censoring time for relapse.
# How to handle this censoring is not straightforward but for simiplicity (as this dataset is only used to demonstrate the methods) we will assume that these individuals remained relapse fee after their censoring time for relapse until their available death time.
rott2 <- rott2 |>
  mutate(exit = os) |>
  mutate(exit = ifelse(rfi == 1, rf, exit)) |>
  mutate(failure = ifelse(osi == 1 | rfi == 1, 1, 0))

# Administrative censoring at 10 years (120 months), to replicate the following stset call:
# . stset exit, failure(failure==1) scale(12) exit(time 120)
rott2 <- rott2 |>
  mutate(failure = ifelse(exit > 120, 0, failure)) |>
  mutate(exit = pmin(exit, 120))

# Then, recode time to years
rott2 <- rott2 |>
  mutate(exit = exit / 12)

# Export in R format
saveRDS(object = rott2, file = "data/rott2.RDS")
