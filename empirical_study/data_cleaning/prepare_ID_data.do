clear all
set more off 
 
*set up global variables containing file path to data directory*
global projdir R:\BSRBR\Analyses\lucy_bull\Serious_infection\
global logdir $projdir/Logfiles
global dodir $projdir/Dofiles
global datadir $projdir/Datasets
 
capture log close 
 
cd "$datadir"
pwd

log using "$logdir/Dofile_ID.log", text replace

*Creating list of IDs including in analysis dataset
use "baseline_serinf_final.dta", clear
*use "serinf_full_study.dta", clear
*sort studyno pyears
*by studyno: keep if _n == 1
keep studyno
rename studyno wantedID

save "study_IDs.dta", replace

use "one_big_file_20181130_for company reports.dta", clear
keep studyno orig_studyno
keep if studyno != orig_studyno
sort studyno
by studyno: gen seq = _n
keep if seq == 1 
gen wantedID = orig_studyno

*Keeping those where orig_studyno is included in the original baseline cohort, only interested in future study IDs to the originals.

merge m:1 wantedID using "study_IDs.dta"
keep if _merge == 3
keep studyno orig_studyno

sort orig_studyno
by orig_studyno: gen seq = _n
count if seq == 2
*62
save "extra_studyno.dta", replace

*Adding the extra study IDs to the wanted list for follow-up information.

use "baseline_serinf_final.dta", clear
*use "serinf_full_study.dta", clear
*sort studyno pyears
*by studyno: keep if _n == 1
keep studyno

*This should not "merge" as such, but add the extra ones to the list
merge 1:1 studyno using "extra_studyno.dta"
keep studyno orig_studyno

*dataset to use in prepare_follow_data do file
save "wanted_IDs.dta", replace

erase study_IDs.dta
erase extra_studyno.dta

log close
