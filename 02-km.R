### Exploring survival differences between treatment groups with the Kaplan-Meier estimator

# Packages
library(tidyverse)
library(survival)
library(ggplot2)
library(pammtools)
library(ragg)

# Load data
rott2 <- readRDS(file = "data/rott2.RDS")

# Obtain Kaplan-Meier survival probabilities by treatment arm (hormon)
# It is important to note that the term survival probabilities does not necessarily refer to being alive or not.
# Instead it refers to being event free.
# In our example, since time to relapse or death (whatever occurred first) is under study, survival probabilities correspond to the probability of being alive without a relapse
km <- survfit(Surv(exit, failure == 1) ~ hormon, data = rott2)

# Extract KM data
km_data <- with(summary(km), tibble(time = time, surv = surv, lower = lower, upper = upper, hormon = strata))

# Add time zero
km_data <- bind_rows(
  km_data,
  distinct(km_data, hormon) |> mutate(time = 0, surv = 1, lower = 1, upper = 1)
)

# Recode treatment factor for plotting
km_data <- km_data |>
  mutate(hormon = factor(
    x = hormon,
    levels = c("hormon=1", "hormon=0"),
    labels = c("Hormonal therapy", "No hormonal therapy")
  ))

# Plot using {ggplot2}
km_plot <- ggplot(km_data, aes(x = time, y = surv, group = hormon)) +
  geom_stepribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_step(aes(linetype = hormon)) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = "Survival probability", x = "Years from surgery", color = "", fill = "", linetype = "") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0, 0),
    legend.justification = c(0, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
km_plot

# Export plot
ggsave(filename = "output/02-km-plot.pdf", plot = km_plot, height = 4, width = 6)
ggsave(filename = "output/02-km-plot.png", plot = km_plot, device = agg_png, height = 4, width = 6, dpi = 300)
