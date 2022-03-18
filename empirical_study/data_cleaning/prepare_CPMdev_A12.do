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

log using "$logdir/Dofile_BaselineCPM2.log", text replace

*~~~~~~~~~~~~
*CPM development data - survival response (time to first serious infection, with censoring at 12 months (horizon time)) 
*~~~~~~~~~~~~
use "serinf_full_study.dta", clear
recode on_steroid_at_baseline .=0

*Variables of interest
keep studyno groupid pgen pyears age firsttreat trtment on_steroid_at_baseline ovmean dascore bmi died renal lung diabetes smoke previous_dmards disdur serinf mfup

gen maxfup = mfup
keep if pyears <= 1
sort studyno pyears 
total(serinf)
by studyno: egen event = max(serinf)
by studyno: egen timevent = max(maxfup) if event == 0
replace timevent = 1 if (timevent > 1 & event == 0)
replace timevent = pyears if serinf == 1
by studyno: egen timeventnew = min(timevent) if event == 1
replace timevent = timeventnew if event == 1
drop timeventnew mfup

*1 row per subject in survival data
by studyno: gen seq = _n
keep if seq == 1

*egen groupid = group(studyno)
rename on_steroid_at_baseline steroids

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

save "baseline_model.dta", replace
saveold "R:\BSRBR\Analyses\lucy_bull\Serious_infection\R Datasets\baseline_model.dta", ///
replace version(12)
*~~~~~~~~~~~~
*Summarising CPM development data
*~~~~~~~~~~~~

sum bmi, det
bys event: sum bmi, det
bys event: sum age, det
bys event: tab pgen
sum age, det
sum disdur, det
bys event: sum disdur, det
sum dascore, det
bys event: sum dascore, det
sum ovmean, det
bys event: sum ovmean, det
tab smoke
bys event: tab smoke
tab steroids
bys event: tab steroids
sum previous_dmards, det
bys event: sum previous_dmards, det
tab lung
tab diabetes
bys event: tab diabetes
bys event: tab lung
tab renal
bys event: tab renal
tab firsttreat
bys event: tab firsttreat
total(timevent)
total(event)
total(timevent) if event == 0 
total(timevent) if event == 1
*total(event) if event == 0
*total(event) if event == 1
stset timevent, id(studyno) failure(event)
stptime
stptime, by(event)

*~~~~~~~~~~~
*CPM development data - Binary response
*~~~~~~~~~~~

use "baseline_model.dta", clear

*no censoring - logistic regression
drop if (timevent < 1 & event == 0)
drop timevent

save "baseline_binary.dta", replace
saveold "R:\BSRBR\Analyses\lucy_bull\Serious_infection\R Datasets\baseline_binary.dta", ///
replace version(12)

*~~~~~~~~~~~~
*Summarising CPM development data
*~~~~~~~~~~~~

sum bmi, det
bys event: sum bmi, det
bys event: sum age, det
bys event: tab pgen
sum age, det
sum disdur, det
bys event: sum disdur, det
sum dascore, det
bys event: sum dascore, det
sum ovmean, det
bys event: sum ovmean, det
tab smoke
bys event: tab smoke
tab steroids
bys event: tab steroids
sum previous_dmards, det
bys event: sum previous_dmards, det
tab lung
tab diabetes
bys event: tab diabetes
bys event: tab lung
tab renal
bys event: tab renal
tab firsttreat
bys event: tab firsttreat
total(event)
total(event) if pgen == 0
total(event) if pgen == 1

*~~~~~~~~~~~~
*CPM development data - survival response (time to first serious infection, with censoring at 12 months (horizon time)) 
*~~~~~~~~~~~~
use "serinf_full_study.dta", clear
recode on_steroid_at_baseline .=0

*Variables of interest
keep studyno groupid pgen pyears age firsttreat trtment on_steroid_at_baseline ovmean dascore bmi died renal lung diabetes smoke previous_dmards disdur serinf mfup

gen maxfup = mfup
keep if pyears <= 4
sort studyno pyears 
total(serinf)
by studyno: egen event = max(serinf)
by studyno: egen timevent = max(maxfup) if event == 0
replace timevent = 4 if (timevent > 4 & event == 0)
replace timevent = pyears if serinf == 1
by studyno: egen timeventnew = min(timevent) if event == 1
replace timevent = timeventnew if event == 1
drop timeventnew mfup

*1 row per subject in survival data
by studyno: gen seq = _n
keep if seq == 1

*egen groupid = group(studyno)
rename on_steroid_at_baseline steroids

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

save "baseline_model2.dta", replace
saveold ".\baseline_model2.dta", ///
replace version(12)

drop studyno
save "baseline_model2_CSF.dta", replace
saveold ".\baseline_model2_CSF.dta", ///
replace version(12)

*~~~~~~~~~~~~
*Summarising CPM development data
*~~~~~~~~~~~~
use "baseline_model2.dta", clear

sum bmi, det
bys event: sum bmi, det
bys event: sum age, det
bys event: tab pgen
sum age, det
sum disdur, det
bys event: sum disdur, det
sum dascore, det
bys event: sum dascore, det
sum ovmean, det
bys event: sum ovmean, det
tab smoke
bys event: tab smoke
tab steroids
bys event: tab steroids
sum previous_dmards, det
bys event: sum previous_dmards, det
tab lung
tab diabetes
bys event: tab diabetes
bys event: tab lung
tab renal
bys event: tab renal
tab firsttreat
bys event: tab firsttreat
total(timevent)
total(event)
total(timevent) if event == 0 
total(timevent) if event == 1
total(event) if pgen == 0
total(event) if pgen == 1
stset timevent, id(studyno) failure(event)
stptime
stptime, by(event)

log close



