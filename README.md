# Standardised Survival Probabilities, in R

This repository contains an R implementation of the Stata code from the manuscript titled [Standardised survival probabilities: a useful and informative tool for reporting regression models for survival data](https://www.nature.com/articles/s41416-022-01949-6), by Syriopoulou et al., and published in 2022 in the British Journal of Cancer.

The content of this repository is organised as follows.

* Code to clean the data is included in the `01-data.R` file;
* Code to re-create the Kaplan-Meier plots is included in the `02-km.R` file;
* Code to re-create the unadjusted predictions and hazard ratios is included in the `03-hr.R` file;
* Code to re-create the adjusted survival probability predictions is included in the `04-surv-adj.R` file;
* Finally, code to re-create the standardised survival probability predictions (and contrast thereof) is included in the `05-surv-std.R` file.

Plots replicating those from the manuscript are included in the `output/` folder, and the original Stata code is included in the `source/` folder as well, for completeness.

# License

This software is released under the [MIT license](https://opensource.org/license/mit).
