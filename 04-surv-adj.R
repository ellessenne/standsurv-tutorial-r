### Adjusted survival curves
### Adjusted survival curves for those who had hormonal therapy and those who did not have hormonal therapy can also be obtained in addition to the hazard ratio obtained from the model output
### A naive way to do this is by setting all adjusting variables to their mean value, so our estimates correspond to the survival of an "average" individual if this "average" individual received hormonal therapy and an "average" individual if they didn't receive hormonal therapy

# Packages
library(tidyverse)
library(rstpm2)
library(ggplot2)
library(ragg)

# Load data
rott2 <- readRDS(file = "data/rott2.RDS")

# Manually create indicator variables for factors
rott2 <- rott2 |>
  mutate(size2 = as.numeric(size == 2), size3 = as.numeric(size == 3))

# Let's now account for imbalances between those who received hormonal therapy and those who didn't by adjusting for several factors in the regression model
# Fit FPM adjusting for size, differantiation grade, number of nodes, progesterone level and age allowing for time-dependent effects for some variables (relaxing the proportional hazards assumption)
# This is not necessarily the best model for this data but it's used for the demonstartion of the methods
# Time-dependent effects are allowed via the `tvc` argument
fpm <- stpm2(Surv(exit, failure == 1) ~ hormon + size2 + size3 + grade + enodes + pr_1 + age, df = 4, data = rott2)
summary(fpm)

# After fitting a FPM, adjusted survival curves by hormonal therapy status (when all other values set to the mean value in the population) is done as follows

# First, we create a new dataset with average values for the model covariates
nd <- summarise(rott2, across(.cols = c("age", "grade", "size2", "size3", "enodes", "pr_1"), .fns = mean))

# Add 100 time points between 0 and 10 years
nd <- crossing(nd, exit = seq(0, 10, length.out = 100))

# Estimate adjusted survival probability for the event of relapse/death for those who did not receive hormonal therapy
nd <- crossing(nd, hormon = c(0, 1))
nd <- cbind(nd, predict(fpm, newdata = nd, type = "surv", se.fit = TRUE))

# Adjusted relapse-free survival probability estimates for breast cancer patients at 10 years
filter(nd, exit == 10)

# Recode treatment factor for plotting
nd <- nd |>
  mutate(hormon = factor(
    x = hormon,
    levels = c(1, 0),
    labels = c("Hormonal therapy", "No hormonal therapy")
  ))

# Plot the survival curves from the unadjusted model, using {ggplot2}
adjs_plot <- ggplot(nd, aes(x = exit, y = Estimate, group = hormon)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_line(aes(linetype = hormon)) +
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
adjs_plot

# Export plot
ggsave(filename = "output/04-adjs-plot.pdf", plot = adjs_plot, height = 4, width = 6)
ggsave(filename = "output/04-adjs-plot.png", plot = adjs_plot, device = agg_png, height = 4, width = 6, dpi = 300)

# Compare flexible parametric model and Cox model (for simplicity we assume no time-dependent effects for this comparison)
# Fit Cox model
cox <- coxph(Surv(exit, failure == 1) ~ hormon + size2 + size3 + grade + enodes + pr_1 + age, data = rott2)

# Predictions for the Cox model
coxnd <- summarise(rott2, across(.cols = c("age", "grade", "size2", "size3", "enodes", "pr_1"), .fns = mean))
coxnd <- crossing(coxnd, hormon = c(0, 1))
coxpreds <- survfit(cox, newdata = coxnd, se.fit = TRUE)

# Create new dataset with predictions from the Cox model
coxpredsdt <- lapply(X = 1:nrow(coxnd), FUN = function(i) {
  crossing(
    coxnd[i, ],
    tibble(
      exit = coxpreds$time,
      Estimate = coxpreds$surv[, i],
      lower = coxpreds$lower[, i],
      upper = coxpreds$upper[, i]
    )
  )
})
# Stack predictions and recode hormon for plotting
coxpredsdt <- bind_rows(coxpredsdt) |>
  mutate(hormon = factor(
    x = hormon,
    levels = c(1, 0),
    labels = c("Hormonal therapy", "No hormonal therapy")
  ))

# Plot comparison
adjs_comp_plot <- ggplot() +
  geom_ribbon(data = nd, aes(x = exit, ymin = lower, ymax = upper, fill = "FPM", group = hormon), alpha = 0.1) +
  geom_ribbon(data = coxpredsdt, aes(x = exit, ymin = lower, ymax = upper, fill = "Cox", group = hormon), alpha = 0.1) +
  geom_line(data = nd, aes(x = exit, y = Estimate, color = "FPM", linetype = hormon)) +
  geom_line(data = coxpredsdt, aes(x = exit, y = Estimate, color = "Cox", linetype = hormon)) +
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
adjs_comp_plot

# Export plot
ggsave(filename = "output/04-adjs-comp-plot.pdf", plot = adjs_comp_plot, height = 4, width = 6)
ggsave(filename = "output/04-adjs-comp-plot.png", plot = adjs_comp_plot, device = agg_png, height = 4, width = 6, dpi = 300)
