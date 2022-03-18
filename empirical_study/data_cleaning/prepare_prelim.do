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

log using "$logdir/Dofile_Prelim.log", text replace

use "serinf_full.dta", clear
sort studyno pyears
by studyno: egen event = max(serinf)
by studyno: egen age2 = max(age)
by studyno: egen timevent = max(pyears) if event == 0
replace timevent = pyears if serinf == 1
by studyno: egen timevent2 = min(timevent) if event == 1
replace timevent = timevent2 if event == 1
drop timevent2
egen groupid = group(studyno)
save "longdat.dta", replace
saveold ".\longdat.dta", ///
replace version(12)

preserve
sort studyno pyears
by studyno: keep if _n == 1
gen days = timevent*365.25
gen lastdate = date + days
format lastdate %td
save "datind.dta", replace
saveold ".\datind.dta", ///
replace version(12)
restore

preserve
keep if dascore !=.
keep studyno pgen age2 pyears dascore event groupid
save "long_DAS.dta", replace
saveold ".\long_DAS.dta", ///
replace version(12)
restore 

preserve
keep if ovmean !=.
keep studyno pgen age2 pyears ovmean event groupid
save "long_HAQ.dta", replace
saveold ".\long_HAQ.dta", ///
replace version(12)
restore

**Under exposure only
use "longdat.dta", clear
drop event timevent
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
drop on_steroid_at_baseline Oralsteriods

merge 1:1 studyno pyears using "trt_breaks.dta"
sort studyno pyears
drop if _merge == 2
drop _merge

*capped follow-up to exposure only
gen maxfup = pyears+0.25 if (largebreak == 1 | firsttreat != trtment) 
by studyno: egen mfup = min(maxfup)
by studyno: egen maxfup2 = max(pyears)
replace mfup = maxfup2 if mfup ==. 
drop if pyears > mfup
replace trtment = firsttreat

*removed subjects who died before time cap
sort studyno pyears
by studyno: egen event = max(serinf)
by studyno: egen died = max(death)
gen timedeath = pyears if death == 1
by studyno: egen diedtime = max(timedeath)
gen dropp = 1 if diedtime < mfup
egen groupid2 = group(studyno) if dropp == 1
sum groupid2, det
drop if dropp == 1
drop largebreak maxfup maxfup2 mfup timedeath diedtime death groupid2 dropp

*events and censoring
by studyno: egen maxfup = max(pyears)
sort studyno pyears 
by studyno: egen timevent = max(maxfup) if event == 0
replace timevent = pyears if serinf == 1
by studyno: egen maxfup2 = max(timevent) if event == 1
replace maxfup = maxfup2 if event == 1
*by studyno: egen timeventnew = min(timevent) if event == 1
*replace timevent = timeventnew if event == 1
*drop timeventnew
gsort studyno -pyears
by studyno: replace timevent = timevent[_n-1] if (timevent ==. & event == 1)
sort studyno pyears
drop if pyears > maxfup
drop maxfup2

recode serinf 0=.
by studyno: gen serinf2 = serinf[_n-1] if serinf ==.
recode serinf .=0
by studyno: replace serinf2 = serinf2[_n-1] if  serinf2 ==.
rename serinf2 prev_sinf
recode prev_sinf .=0

by studyno: egen numev = sum(serinf)
drop if timevent <=0


save "longdata_exp.dta", replace
saveold ".\longdata_exp.dta", ///
replace version(12)


use "longdata_exp.dta", clear
sort studyno pyears
by studyno: keep if _n == 1
save "data_exp_ind.dta", replace
saveold ".\data_exp_ind.dta", ///
replace version(12)

use "longdata_exp.dta", clear
sort studyno pyears
by studyno: keep if _n == 1
saveold ".\longdata_ind.dta", ///
replace version(12)

use "longdata_exp.dta", clear
keep if dascore !=.
keep studyno pgen age2 pyears dascore event timevent groupid
keep if pyears <= 4
replace event = 0 if (timevent > 4 & event == 1)
replace timevent = 4 if timevent >= 4

save "long_DAS_exp.dta", replace
saveold ".\long_DAS_exp.dta", ///
replace version(12)
keep if pyears <= 1.5
save "long_DAS_exp1.dta", replace
saveold ".\long_DAS_exp1.dta", ///
replace version(12)

use "longdata_exp.dta", clear
keep if ovmean !=.
keep studyno pgen age2 pyears ovmean event timevent groupid
keep if pyears <= 4
replace event = 0 if (timevent > 4 & event == 1)
replace timevent = 4 if timevent >= 4
save "long_HAQ_exp.dta", replace
saveold ".\long_HAQ_exp.dta", ///
replace version(12)
keep if pyears <= 1.5
save "long_HAQ_exp.dta", replace
saveold ".\long_HAQ_exp1.dta", ///
replace version(12)

use "longdata_exp.dta", clear
drop if pyears == 0
sort studyno pyears
by studyno: egen sumdiab = max(diabetes)
by studyno: egen sumren = max(renal)
by studyno: egen sumlung = max(lung)
by studyno: keep if _n == 1

replace event = 0 if (timevent > 4 & event == 1)
replace timevent = 4 if timevent >= 4

save "com_fup.dta", replace
saveold ".\com_fup.dta", ///
replace version(12)

use "longdata_exp.dta", clear
sort studyno pyears
by studyno: egen sumdiab = max(diabetes)
by studyno: egen sumren = max(renal)
by studyno: egen sumlung = max(lung)
by studyno: keep if _n == 1

replace event = 0 if (timevent > 4 & event == 1)
replace timevent = 4 if timevent >= 4

save "com_exp.dta", replace
saveold ".\com_exp.dta", ///
replace version(12)

use "longdata_exp.dta", clear

keep if pyears <=4
replace event = 0 if (timevent > 4 & event == 1)
replace timevent = 4 if timevent >= 4

keep if event == 1 
gen diabtime = pyears if diabetes == 1
gen lungtime = pyears if lung == 1
gen renaltime = pyears if renal == 1
sort studyno pyears
by studyno: egen timediab = min(diabtime)
by studyno: egen timelung = min(lungtime)
by studyno: egen timerenal = min(renaltime)
save "longdata_P.dta", replace

use "longdata_P.dta", clear
drop if (timediab ==. | timediab ==0)
gen distD = (timevent - timediab)
replace distD =. if distD < 0
by studyno: egen distDD = min(distD)
by studyno: keep if _n == 1
keep groupid distD distDD
save "diab_times.dta", replace
saveold ".\diab_times.dta", ///
replace version(12)


use "longdata_P.dta", clear
drop if (timelung ==. | timelung ==0)
gen distL = (timevent - timelung)
replace distL =. if distL < 0
by studyno: egen distLL = min(distL)
by studyno: keep if _n == 1
keep groupid distL distLL
save "lung_times.dta", replace
saveold ".\lung_times.dta", ///
replace version(12)


use "longdata_P.dta", clear
drop if (timerenal ==. | timerenal ==0)
gen distR = (timevent - timerenal)
replace distR =. if distR < 0
by studyno: egen distRR = min(distR)
by studyno: keep if _n == 1
keep groupid distR distRR
save "renal_times.dta", replace
saveold ".\renal_times.dta", ///
replace version(12)

log close
