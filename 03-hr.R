### Estimating hazard ratios using regression models
### We will now fit a regression model including only hormonal therapy
### The most commonly applied statistical regression model when studying time-to-event outcomes is Cox's proportional hazards model but here we will fit a flexible parametric model (FPM)
### FPMs have advantages in terms of modelling time-dependent effects (i.e. relaxing the assumption of proportional hazards) and also in terms of predictions
### FPMs explicitly estimate the baseline log cumulative hazard by using restricted cubic splines for the logarithm of time t rather than assuming linearity with time
### The complexity of the available data dictates the number of knots (or the number of degrees of freedom (df) that is equal to the number of knots minus 1) needed to accurately capture the shape for the underlying hazard.

# Packages
library(tidyverse)
library(rstpm2)
library(ggplot2)

# Load data
rott2 <- readRDS(file = "data/rott2.RDS")

# We fit the FPM model using the {rstpm2} package
fpm <- stpm2(Surv(exit, failure == 1) ~ hormon, df = 5, data = rott2)
# Argument `df` denotes the degrees of freedom used for the splines

# Model table
summary(fpm)
# In the output, nsx1-nsx5 refer to the splines used to model the baseline hazard

# Obtain survival probabilites from the unadjusted model at 1 year after diagnosis
nd <- crossing(exit = 1, hormon = c(0, 1))
predict(fpm, newdata = nd, type = "surv")

# Obtain survival probabilites from the unadjusted model at 10 years after diagnosis
nd <- crossing(exit = 10, hormon = c(0, 1))
predict(fpm, newdata = nd, type = "surv")

# We also create timevar variable with 100 observations taking values from 0 to 10 years to use for the survival estimates across the total follow-up
nd <- crossing(exit = seq(0, 10, length.out = 100), hormon = c(0, 1))
nd <- cbind(nd, predict(fpm, newdata = nd, type = "surv", se.fit = TRUE))

# Recode treatment factor for plotting
nd <- nd |>
  mutate(hormon = factor(
    x = hormon,
    levels = c(1, 0),
    labels = c("Hormonal therapy", "No hormonal therapy")
  ))

# Plot the survival curves from the unadjusted model, using {ggplot2}
unadjs_plot <- ggplot(nd, aes(x = exit, y = Estimate, group = hormon)) +
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
unadjs_plot

# Export plot
ggsave(filename = "output/03-unadjs-plot.pdf", plot = unadjs_plot, height = 4, width = 6)
