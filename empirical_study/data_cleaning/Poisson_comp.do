clear all
set more off 
 
*set up global variables containing file path to data directory*
global projdir ~\Serious_infection\
global logdir $projdir/Logfiles
global dodir $projdir/Dofiles
global datadir $projdir/Datasets
 
capture log close 
 
cd "$datadir"
pwd

log using "$logdir/Dofile_PoissonComp.log", text replace

*~~~~~~~~~~~~
*CPM development data - count response (count of serious infections in exposure period) 
*~~~~~~~~~~~~


use "serinf_full_study.dta", clear
by studyno: egen countevent = sum(serinf)if pyears <= 4
drop if pyears > timevent
drop if timevent <= 0

gen trt1 = 1 if trtment == 1
gen trt2 = 1 if trtment == 2
gen trt4 = 1 if trtment == 4
gen trt1561 = 1 if trtment == 1561
gen trt10060 = 1 if trtment == 10060
recode trt1 .=0
recode trt2 .=0
recode trt4 .=0
recode trt1561 .=0
recode trt10060 .=0

gen smoke0 = 1 if smoke == 0
gen smoke1 = 1 if smoke == 1
gen smoke2 = 1 if smoke == 2
recode smoke0 .=0
recode smoke1 .=0
recode smoke2 .=0

gen roundtime = 1 if (timevent > 0 & timevent <= 1)
replace roundtime = 2 if (timevent > 1 & timevent <= 2)
replace roundtime = 3 if (timevent > 2 & timevent <= 3)
replace roundtime = 4 if (timevent > 3 & timevent <= 4)

save "poisson_model.dta", replace
saveold "R:\BSRBR\Analyses\lucy_bull\Serious_infection\R Datasets\poisson_model.dta", ///
replace version(12)

*~~~~~~~~~~~~
*CPM development data - cloglog (prob at least one per year) 
*~~~~~~~~~~~~

use "baseline_model2.dta", clear

replace event = 0 if (timevent > 4 & event == 1)
replace timevent = 4 if timevent > 4
gen roundtime = 1 if (timevent > 0 & timevent <= 1)
replace roundtime = 2 if (timevent > 1 & timevent <= 2)
replace roundtime = 3 if (timevent > 2 & timevent <= 3)
replace roundtime = 4 if (timevent > 3 & timevent <= 4)
gen steroids = on_steroids_at_baseline

save "cloglog_model.dta", replace
saveold ".\cloglog_model.dta", ///
replace version(12)

log close
