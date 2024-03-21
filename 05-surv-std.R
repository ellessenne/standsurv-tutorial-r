### Standardised survival curves
### A more useful way of obtaining adjusted survival probabilities is to apply regression standardisation
### To do so, we estimate for each individual in the study population their survival probability given the individuals covariate pattern and then average the individual specific estimates to obtain the standardised survival probability
### In Stata this can be done with command -standsurv-; here, we implement standardised predictions by hand

# Packages
library(tidyverse)
library(rstpm2)
library(ggplot2)
library(cowplot)

# Load data
rott2 <- readRDS(file = "data/rott2.RDS")

# Manually create indicator variables for factors
rott2 <- rott2 |>
  mutate(size2 = as.numeric(size == 2), size3 = as.numeric(size == 3))

# Fit the same flexible parametric model from the previous script
fpm <- stpm2(Surv(exit, failure == 1) ~ hormon + size2 + size3 + grade + enodes + pr_1 + age, df = 4, data = rott2)

# To obtain standardised predictions, we need to predict for hormon=0 and hormon=1, for all time points
# We start by predicting at 10 years, which we later generalise.
pred_h1 <- mean(predict(fpm, newdata = mutate(rott2, hormon = 1, exit = 10), type = "surv"))
pred_h0 <- mean(predict(fpm, newdata = mutate(rott2, hormon = 0, exit = 10), type = "surv"))
cbind(pred_h1 = pred_h1, pred_h0 = pred_h0, diff = pred_h1 - pred_h0)
# Looks consistent with the paper!
# Now, we define a function that can do that for us, for any number of timepoints
# This is specific to this application, so will have to be adjusted/redefined for other prediction settings
standsurv <- function(fit, newdata, timevar, hormon) {
  nd <- newdata
  nd$hormon <- hormon
  nd <- lapply(X = timevar, FUN = function(x) {
    this <- nd
    this$exit <- x
    this
  })
  nd <- do.call(rbind.data.frame, nd)
  nd$s <- predict(fit, newdata = nd, type = "surv")
  nd <- group_by(nd, exit) |>
    summarise(s = mean(s)) |>
    pull(s)
  return(nd)
}
# Note that this function is not really general nor optimised, so it is somewhat slow, but it does the job for now and is easy to adapt to other settings.
# Then, we use the function above to create predictions for time 0 to 10:
timevar <- seq(0, 10, length.out = 100)
pred_h0 <- standsurv(fit = fpm, newdata = rott2, timevar = timevar, hormon = 0)
pred_h1 <- standsurv(fit = fpm, newdata = rott2, timevar = timevar, hormon = 1)
# We also define a function for predicting the difference
standsurv_diff <- function(value, ref, ...) {
  sv <- standsurv(..., hormon = value)
  sr <- standsurv(..., hormon = ref)
  return(sv - sr)
}
pred_diff <- standsurv_diff(value = 1, ref = 0, fit = fpm, newdata = rott2, timevar = timevar)
# Plot the standardised survival probabilities under treatment arms, using {ggplot2}
predictions <- bind_rows(
  tibble(exit = timevar, s = pred_h0, hormon = 0),
  tibble(exit = timevar, s = pred_h1, hormon = 1)
) |>
  mutate(hormon = factor(
    x = hormon,
    levels = c(1, 0),
    labels = c("Hormonal therapy", "No hormonal therapy")
  ))
# Then, plot:
std_plot <- ggplot(predictions, aes(x = exit, y = s, linetype = hormon)) +
  geom_line() +
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
std_plot
# Same, for the survival difference:
stdd_plot <- tibble(exit = timevar, difference = pred_diff) |>
  ggplot(aes(x = exit, y = difference)) +
  geom_line(linetype = "dotdash") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = "Difference in survival probability", x = "Years from surgery", color = "", fill = "", linetype = "") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0, 0),
    legend.justification = c(0, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
stdd_plot

# We combine the two plots in one via the {cowplot} package
std_stdd_plot <- plot_grid(std_plot, stdd_plot, align = "hv", axis = "tblr", nrow = 1, labels = letters)
std_stdd_plot

# Export plots
ggsave(filename = "output/05-std-plot.pdf", plot = std_plot, height = 4, width = 6)
ggsave(filename = "output/05-stdd-plot.pdf", plot = stdd_plot, height = 4, width = 6)
ggsave(filename = "output/05-std-stdd-plot.pdf", plot = std_stdd_plot, height = 4, width = 7)

# Note that the predictions above do not include confidence intervals, as we did not calculate them.
# For this, we can use the predictnl() function from the {rstpm2} package.
# This uses the numerical delta method:
pred_h0_se <- predictnl(fpm, newdata = rott2, fun = standsurv, timevar = timevar, hormon = 0)
pred_h0_se <- cbind(pred_h0_se, confint(pred_h0_se))
pred_h1_se <- predictnl(fpm, newdata = rott2, fun = standsurv, timevar = timevar, hormon = 1)
pred_h1_se <- cbind(pred_h1_se, confint(pred_h1_se))
pred_diff_se <- predictnl(fpm, newdata = rott2, fun = standsurv_diff, value = 1, ref = 0, timevar = timevar)
pred_diff_se <- cbind(pred_diff_se, confint(pred_diff_se))

# Now, recreate the plots but with confidence intervals:
predictions_se <- bind_rows(
  tibble(exit = timevar, s = pred_h0_se$fit, lower = pred_h0_se$`2.5 %`, upper = pred_h0_se$`97.5 %`, hormon = 0),
  tibble(exit = timevar, s = pred_h1_se$fit, lower = pred_h1_se$`2.5 %`, upper = pred_h1_se$`97.5 %`, hormon = 1)
) |>
  mutate(hormon = factor(
    x = hormon,
    levels = c(1, 0),
    labels = c("Hormonal therapy", "No hormonal therapy")
  ))
# Then, plot:
stdse_plot <- ggplot(predictions_se, aes(x = exit, y = s, linetype = hormon)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_line() +
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
stdse_plot
# Same, for the survival difference:
stddse_plot <- tibble(exit = timevar, difference = pred_diff_se$fit, lower = pred_diff_se$`2.5 %`, upper = pred_diff_se$`97.5 %`) |>
  ggplot(aes(x = exit, y = difference)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_line(linetype = "dotdash") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = "Difference in survival probability", x = "Years from surgery", color = "", fill = "", linetype = "") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0, 0),
    legend.justification = c(0, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
stddse_plot
# Same for the combined plot
stdse_stddse_plot <- plot_grid(stdse_plot, stddse_plot, align = "hv", axis = "tblr", nrow = 1, labels = letters)
stdse_stddse_plot

# Export these new plots
ggsave(filename = "output/05-stdse-plot.pdf", plot = stdse_plot, height = 4, width = 6)
ggsave(filename = "output/05-stddse-plot.pdf", plot = stddse_plot, height = 4, width = 6)
ggsave(filename = "output/05-stdse-stddse-plot.pdf", plot = stdse_stddse_plot, height = 4, width = 7)

# We can also calculate the ratio instead of the difference:
standsurv_ratio <- function(value, ref, ...) {
  sv <- standsurv(..., hormon = value)
  sr <- standsurv(..., hormon = ref)
  return(sv / sr)
}
pred_ratio_se <- predictnl(fpm, newdata = rott2, fun = standsurv_ratio, value = 1, ref = 0, timevar = timevar)
pred_ratio_se <- cbind(pred_ratio_se, confint(pred_ratio_se))

# ... and plot it
stdrse_plot <- tibble(exit = timevar, ratio = pred_ratio_se$fit, lower = pred_ratio_se$`2.5 %`, upper = pred_ratio_se$`97.5 %`) |>
  ggplot(aes(x = exit, y = ratio)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_line(linetype = "dotdash") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = "Ratio of survival probability", x = "Years from surgery", color = "", fill = "", linetype = "") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0, 0),
    legend.justification = c(0, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
stdrse_plot
ggsave(filename = "output/05-stdrse-plot.pdf", plot = stdrse_plot, height = 4, width = 6)

### Finally, we show how to standardise within a subset of the study population
### Above standardisation was performed using the empirical covariate distribution in the whole population i.e. we estimated the average relapse-free survival probability within the whole population, if everyone was treated compared to if no one was treated
### However, in certain cases it will be more relevant to apply the empirical covariate distribution of a subset of the total study population such as the distribution among treated
### For this, we standardise to the subset of study subjects with hormon=1:
predsub_h0_se <- predictnl(fpm, newdata = filter(rott2, hormon == 1), fun = standsurv, timevar = timevar, hormon = 0)
predsub_h0_se <- cbind(predsub_h0_se, confint(predsub_h0_se))
predsub_h1_se <- predictnl(fpm, newdata = filter(rott2, hormon == 1), fun = standsurv, timevar = timevar, hormon = 1)
predsub_h1_se <- cbind(predsub_h1_se, confint(predsub_h1_se))
predsub_diff_se <- predictnl(fpm, newdata = filter(rott2, hormon == 1), fun = standsurv_diff, value = 1, ref = 0, timevar = timevar)
predsub_diff_se <- cbind(predsub_diff_se, confint(predsub_diff_se))

# Plots of these:
predictionssub_se <- bind_rows(
  tibble(exit = timevar, s = predsub_h0_se$fit, lower = predsub_h0_se$`2.5 %`, upper = predsub_h0_se$`97.5 %`, hormon = 0),
  tibble(exit = timevar, s = predsub_h1_se$fit, lower = predsub_h1_se$`2.5 %`, upper = predsub_h1_se$`97.5 %`, hormon = 1)
) |>
  mutate(hormon = factor(
    x = hormon,
    levels = c(1, 0),
    labels = c("Hormonal therapy", "No hormonal therapy")
  ))
# Then, plot:
stdsesub_plot <- ggplot(predictionssub_se, aes(x = exit, y = s, linetype = hormon)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_line() +
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
stdsesub_plot
# Same, for the survival difference:
stddsesub_plot <- tibble(exit = timevar, difference = predsub_diff_se$fit, lower = predsub_diff_se$`2.5 %`, upper = predsub_diff_se$`97.5 %`) |>
  ggplot(aes(x = exit, y = difference)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) +
  geom_line(linetype = "dotdash") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  labs(y = "Difference in survival probability", x = "Years from surgery", color = "", fill = "", linetype = "") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0, 0),
    legend.justification = c(0, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
stddsesub_plot
# Same for the combined plot
stdsesub_stddsesub_plot <- plot_grid(stdsesub_plot, stddsesub_plot, align = "hv", axis = "tblr", nrow = 1, labels = letters)
stdsesub_stddsesub_plot

# Export these plots
ggsave(filename = "output/05-stdsesub-plot.pdf", plot = stdsesub_plot, height = 4, width = 6)
ggsave(filename = "output/05-stddsesub-plot.pdf", plot = stddsesub_plot, height = 4, width = 6)
ggsave(filename = "output/05-stdsesub-stddsesub-plot.pdf", plot = stdsesub_stddsesub_plot, height = 4, width = 7)

# Compare predictions standardised over the entire population vs predictions standardised over the subset with hormon=1:
predictions_comp <- bind_rows(
  predictions_se |> mutate(std = "Standardised to entire population"),
  predictionssub_se |> mutate(std = "Standardised to treated population")
)
# Then, plot:
stdcomp_plot <- ggplot(predictions_comp, aes(x = exit, y = s, linetype = hormon)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = std), alpha = 0.1) +
  geom_line(aes(color = std)) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
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
stdcomp_plot
# Difference
stddcomp_plot <- bind_rows(
  tibble(exit = timevar, difference = pred_diff_se$fit, lower = pred_diff_se$`2.5 %`, upper = pred_diff_se$`97.5 %`, std = "Standardised to entire population"),
  tibble(exit = timevar, difference = predsub_diff_se$fit, lower = predsub_diff_se$`2.5 %`, upper = predsub_diff_se$`97.5 %`, std = "Standardised to treated population")
) |>
  ggplot(aes(x = exit, y = difference)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = std), alpha = 0.1) +
  geom_line(aes(color = std), linetype = "dotdash") +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "Difference in survival probability", x = "Years from surgery", color = "", fill = "", linetype = "") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(1, 0),
    legend.justification = c(1, 0),
    legend.key = element_blank(),
    legend.background = element_blank()
  )
stddcomp_plot
# Same for the combined plot
stdcomp_stddcomp_plot <- plot_grid(stdcomp_plot, stddcomp_plot + theme(legend.position = "none"), align = "hv", axis = "tblr", nrow = 1, labels = letters)
stdcomp_stddcomp_plot

# Export these plots
ggsave(filename = "output/05-stdcomp-plot.pdf", plot = stdcomp_plot, height = 4, width = 6)
ggsave(filename = "output/05-stddcomp-plot.pdf", plot = stddcomp_plot, height = 4, width = 6)
ggsave(filename = "output/05-stdcomp-stddcomp-plot.pdf", plot = stdcomp_stddcomp_plot, height = 5, width = 7)
