/*
// Stata code for the illustrative example of the paper "Standardised survival probabilities: a useful and informative tool for reporting regression models for survival data"
// The analysis was run using Stata 17

// This analysis uses the rotterdam data set (rott2.dta) that includes 2982 primary breast cancers patients whose data whose records were included in the Rotterdam tumor bank.
// Data can be downloaded using the following commands: 

net from http://www.stata-press.com/data/fpsaus/
net get fpsaus-dta

// Finally, in the regression models we will fit below, we adopted the same approach as in previous analyses of this data that have found that transformation of the nodes variable (exp(-0.12*nodes)) and the progesterone variable (log(pr + 1)) model the non-linear effects of these variables fairly well

// For the analysis, we use some user-written Stata commands. These can be installed within Stata from the Boston College Statistical Software Components (SSC) archive as follows:
// To fit the flexible parametric survival models
ssc install stpm2, replace
// To generate the restricted cubic spline functions (this command is required to be able to use stpm2)
ssc install rcsgen, replace 
// To obtain standardised survival probabilities using regression standardisation 
ssc install standsurv, replace
*/


// Load data
use rott2.dta, clear

// Generate follow-up time time as the time to death or time to relapse (depending on which event occurred first)
// Note: There are 43 subjects who have died without relapse, but their time of death is greater than the censoring time for relapse. 
// How to handle this censoring is not straightforward but for simiplicity (as this dataset is only used to demonstrate the methods) we will assume that these individuals remained relapse fee after their censoring time for relapse until their available death time.
gen exit=os
replace exit=rf if rfi==1
// Generate failure variable 
gen failure=0
replace failure=1 if osi==1 | rfi==1


// Declare survival data with the outcome of interest being time to relapse or death (whichever occurred first) and exit time is either relapse time, death time or 120 months (10 years) since surgery (whatever occurred first)
stset exit, failure(failure==1) scale(12) exit(time 120)



//------- Exploring survival differences between treatment groups with the Kaplan-Meier estimator -------//

// Obtain Kaplan-Meier survival probabilities by treatment arm (hormon)
// It is important to note that the term survival probabilities does not necessarily refer to being alive or not. Instead it refers to being event free. In our example, since time to relapse or death (whatever occurred first) is under study, survival probabilities correspond to the probability of being alive without a relapse
sts list, by(hormon) compare 
// Plot Kaplan-Meier survival estimates
sts graph, by(hormon) name(KM, replace) plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	legend(order( 2 "hormonal therapy" 1 "no hormonal therapy") ring(0) pos(7) cols(1) region(lcolor(white))) ylabel(, angle(0)) ///
	title("") ytitle("Survival probability") ///
	xtitle("Years from surgery") ///
	plot1opts(lpattern(dash) lcolor(black)) ///
	plot2opts(lpattern(solid) lcolor(black))
	

	
//------- Estimating hazard ratios using regression models -------//
	
// We will now fit a regression model including only hormonal therapy
// The most commonly applied statistical regression model when studying time-to-event outcomes is Cox's proportional hazards model but here we will fit a flexible parametric model (FPM) 
// FPMs have advantages in terms of modelling time-dependent effects (i.e. relaxing the assumption of proportional hazards) and also in terms of predictions
// FPMs explicitly estimate the baseline log cumulative hazard by using restricted cubic splines for the logarithm of time t rather than assuming linearity with time 
// The complexity of the available data dictates the number of knots (or the number of degrees of freedom (df) that is equal to the number of knots minus 1) needed to accurately capture the shape for the underlying hazard. 
// Option df() denotes the degrees of freedom used for the splines
stpm2 hormon, scale(hazard) df(5) eform
// In the output, rcs1 – rcs5 refers to the splines used to model the baseline hazard
// obtain survival probabilites from the unadjusted model at 1 year after diagnosis
gen t1=1 if _n==1
predict s_unadj_h0_t1, surv at(hormon 0) timevar(t1)  
predict s_unadj_h1_t1, surv at(hormon 1) timevar(t1)  
li  s_unadj_h0_t1 s_unadj_h1_t1 if t1==1
// obtain survival probabilites from the unadjusted model at 10 years after diagnosis
gen t10=10 if _n==1
predict s_unadj_h0_t10, surv at(hormon 0) timevar(t10)  
predict s_unadj_h1_t10, surv at(hormon 1) timevar(t10)  
li  s_unadj_h0_t10 s_unadj_h1_t10 if t10==10
// We also create timevar variable with 100 observations taking values from 0 to 10 years to use for the survival estimates across the total follow-up
range timevar 0 10 100
// obtain survival probabilites from the unadjusted model
predict s_unadj_h0, surv at(hormon 0) timevar(timevar) ci 
predict s_unadj_h1, surv at(hormon 1) timevar(timevar) ci 
// Plot the survival curves from the unadjusted model
twoway (line s_unadj_h0 timevar, sort lcolor(black) lpattern(dash)) ///
	(line s_unadj_h1 timevar, sort lcolor(black) lpattern(solid)) ///
	(rarea s_unadj_h0_lci s_unadj_h0_uci timevar, sort color(black%10)) ///
	(rarea s_unadj_h1_lci s_unadj_h1_uci timevar, sort color(black%10)) ///
	, xtitle("Years from surgery") ///
	ytitle("Survival probbaility") ///
	ylabel(0(.25)1,angle(h)) ///
	name(adj1_naive, replace) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	legend(order( 2 "hormonal therapy" 1 "no hormonal therapy") ring(0) pos(7) cols(1) region(lcolor(white))) ylabel(0(0.25)1, format(%6.2f) angle(0)) 
	



// Let's now account for imbalances between those who received hormonal therapy and those who didn't by adjusting for several factors in the regression model
// First we need to create dummy variables for the size variable - stpm2 command does not support factor variable and requires categorical variables with more than 2 levels to be included in the model as dummy variables
tab size, gen(size)

// Fit FPM adjusting for size, differantiation grade, number of nodes, progesterone level and age allowing for time-dependent effects for some variables (relaxing the proportional hazards assumption)
// This is not necessarily the best model for this data but it's used for the demonstartion of the methods
// Time-dependent effects are allowed for the variables included in option tvc()
// Option dftvc() should also be given for the degrees of freedom required to the time-dependet effects 
stpm2 hormon size2 size3 grade enodes pr_1 age, scale(hazard) df(4) eform tvc(enodes grade) dftvc(3)

//------- Adjusted survival curves -------//

// Adjusted survival curves for those who had hormonal therapy and those who did not have hormonal therapy can also be obtained in addition to the hazard ratio obtained from the model output
// A naive way to do this is by setting all adjusting variables to their mean value, so our estimates correspond to the survival of an "average" individual if this "average" individual received hormonal therapy and an "average" individual if they didn't receive hormonal therapy

// After fitting a FPM, adjusted survival curves by hormonal therapy status (when all other values set to the mean value in the population) is done as follows
// store average covariate value
foreach var in age grade size2 size3 enodes pr_1 {
         summ `var', meanonly
         local atopt `atopt' `var' `r(mean)'
		 di "`r(mean)'"
 }
// Estimate adjusted survival probability for the event of relapse/death for those who did not receive hormonal therapy
predict s_naive_h0, surv at(hormon 0 `atopt') timevar(timevar) ci 
// Estimate adjusted survival probability for the event of relapse/death for those who received hormonal therapy
predict s_naive_h1, surv at(hormon 1 `atopt') timevar(timevar) ci 

// Adjusted relapse-free survival probability estimates for breast cancer patients at 10 years
li s_naive_h0* s_naive_h1* if timevar==10, abb(30)

// Plot the survival curves
twoway (line s_naive_h0 timevar, sort lcolor(black) lpattern(dash)) ///
	(line s_naive_h1 timevar, sort lcolor(black) lpattern(solid)) ///
	(rarea s_naive_h0_lci s_naive_h0_uci timevar, sort color(black%10)) ///
	(rarea s_naive_h1_lci s_naive_h1_uci timevar, sort color(black%10)) ///
	, xtitle("Years from surgery") ///
	ytitle("Survival probbaility") ///
	ylabel(0(.25)1,angle(h)) ///
	name(adj1_naive, replace) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	legend(order( 2 "hormonal therapy" 1 "no hormonal therapy") ring(0) pos(7) cols(1) region(lcolor(white))) ylabel(0(0.25)1, format(%6.2f) angle(0)) 
	

	
/*
// Compare flexible parametric model and Cox model (for simplicity we assume no time-dependent effects for this comparison)
// Fit a FPM
stpm2 hormon age grade size2 size3 enodes pr_1, scale(hazard) df(5) eform
// Obtain survival estimates for those who had hormonal therapy using the average covariate value (i.e. for the average individual)
predict s_fpm, surv at(hormon 1 `atopt') timevar(timevar)
// Fit Cox model
stcox hormon age i.grade i.size enodes pr_1
// Compare the survival probaility for those who had hormonal therapy obtained from the FPM to the one obtained from the Cox model 
stcurve, survival at(hormon=1) addplot(line s_fpm timevar, sort)
*/
	
//------- Standardised survival curves -------//

// A more useful way of obtaining adjusted survival probabilities is to apply regression standardisation 
// To do so, we estimate for each individual in the study population their survival probability given the individuals covariate pattern and then average the individual specific estimates to obtain the standardised survival probability
// In Stata this can be done with command standsurv
standsurv, at1(hormon 0) at2(hormon 1) timevar(timevar) ci atvars(s_h0 s_h1) contrast(difference) contrastvars(dif)
// Standardised relapse-free survival probabilities for breast cancer patients at 10 years
li s_h0* s_h1* dif* if timevar==10, abb(30)

// Plot the standardised survival probabilities under treatment arms 
twoway (line s_h0 timevar, sort lcolor(black) lpattern(dash)) ///
	(line s_h1  timevar, sort lcolor(black) lpattern(solid)) ///
	(rarea s_h0_lci s_h0_uci timevar, sort color(black%10)) ///
	(rarea s_h1_lci s_h1_uci timevar, sort color(black%10)) ///
	, xtitle("Years from surgery") ///
	ytitle("Survival probability") ///
	ylabel(0(.25)1,angle(h)) ///
	 name(adj1, replace) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	legend(order( 2 "hormonal therapy" 1 "no hormonal therapy") ring(0) pos(7) cols(1) region(lcolor(white))) ylabel(0(0.25)1, format(%6.2f) angle(0)) ///
	title("a", span color(black) pos(11))
// Plot the difference in standardised survival probabilities under hormonal theray and no hormonal therapy		
twoway (line dif timevar, sort  lcolor(black) lpattern(shortdash)) ///
	(rarea dif_lci dif_uci timevar, sort color(black%10)) ///
	, xtitle("Years from surgery") ///
	ytitle("Difference in survival probabilities") ///
	ylabel(0(.02)0.14, format(%6.2f)angle(h)) ///
	legend(off) name(difference, replace) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	title("b", span color(black) pos(11))


graph combine adj1 difference, plotregion(fcolor(white)) graphregion(fcolor(white))



// We can also calculate the ratio instead of the difference
standsurv, at1(hormon 0) at2(hormon 1) timevar(timevar) ci atvars(surv_h0 surv_h1) contrast(ratio) contrastvars(ratio)
// Standardised relapse-free survival probabilities for breast cancer patients at 10 years
li ratio* if timevar==10, abb(30)
// Plot the ratio of standardised survival probabilities under hormonal theray and no hormonal therapy		
twoway (line ratio timevar, sort  lcolor(black) lpattern(shortdash)) ///
	(rarea ratio_lci ratio_uci timevar, sort color(black%10)) ///
	, xtitle("Years from surgery") ///
	ytitle("Ratio of survival probabilities") ///
	ylabel(, format(%6.2f)angle(h)) ///
	legend(off) name(ratio, replace) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) 
	


//------- Standardising within a subset of the study population -------//
// Above standardisation was performed using the empirical covariate distribution in the whole population i.e. we estimated the average relapse-free survival probability within the whole population, if everyone was treated compared to if no one was treated
// However, in certain cases it will be more relevant to apply the empirical covariate distribution of a subset of the total study population such as the distribution among treated
// This can easily be done in standsurv with the option if hormon==1
standsurv if hormon==1, at1(hormon 0) at2(hormon 1) timevar(timevar) ci atvars(s_h0_sub s_h1_sub) contrast(difference) contrastvars(dif_sub)
// Standardised relapse-free survival probabilities among breast cancer patients who received hormonal therapy at 10 years
li s_h0_sub* s_h1_sub* dif_sub* if timevar==10, abb(30)



// Plot the survival probabilities under treatment arms when standarised among those who had hormonal therapy
twoway (line s_h0_sub timevar, sort lcolor(black) lpattern(dash)) ///
	(line s_h1_sub  timevar, sort lcolor(black) lpattern(solid)) ///
	(rarea s_h0_sub_lci s_h0_sub_uci timevar, sort color(black%10)) ///
	(rarea s_h1_sub_lci s_h1_sub_uci timevar, sort color(black%10)) ///
	, xtitle("Years from surgery") ///
	ytitle("Survival probability") ///
	ylabel(0(.25)1,angle(h)) ///
	 name(adj1_sub, replace) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	legend(order( 2 "hormonal therapy" 1 "no hormonal therapy") ring(0) pos(7) cols(1) region(lcolor(white))) ylabel(0(0.25)1, format(%6.2f) angle(0)) ///
	title("a", span color(black) pos(11))
	
// Plot the difference in survival probabilities under hormonal therapy and no hormonal therapy when standarised among those who had hormonal therapy
// This corresponds to the difference in relapse-free survival probabilities of the patients who had hormonal therapy that is due to hormonal therapy
twoway (line dif_sub timevar, sort  lcolor(black) lpattern(shortdash)) ///
	(rarea dif_sub_lci dif_sub_uci timevar, sort color(black%10)) ///
	, xtitle("Years from surgery") ///
	ytitle("Difference in survival probabilities") ///
	ylabel(0(.02)0.14, format(%6.2f)angle(h)) ///
	legend(off) name(difference_sub, replace) ///
	plotregion(fcolor(white)) graphregion(fcolor(white)) ///
	title("b", span color(black) pos(11))


graph combine adj1_sub difference_sub, plotregion(fcolor(white)) graphregion(fcolor(white))


