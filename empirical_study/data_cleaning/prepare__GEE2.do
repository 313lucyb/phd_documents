*GEE model predictions in Stata

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

log using "$logdir/Dofile_GEE2.log", text replace

use "long_TDCM2.dta", clear
drop if start > 4
replace serinf = 0 if (end > 4 & serinf == 1)
replace end = 4 if end > 4
gen interval = end - start

gen roundtime = 1 if (end >= 0 & end <= 1)
replace roundtime = 2 if (end > 1 & end <= 2)
replace roundtime = 3 if (end > 2 & end <= 3)
replace roundtime = 4 if (end > 3 & end <= 4)
drop if interval == 0
save "long_GEE2.dta", replace
saveold ".\long_GEE.dta", ///
replace version(12)

log close
