*Creating data for JM model

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

log using "$logdir/Dofile_JM2.log", text replace

use "serinf_full_study.dta", clear
keep studyno groupid pgen pyears age firsttreat trtment on_steroid_at_baseline ovmean dascore bmi renal lung diabetes smoke previous_dmards disdur steroids Oralsteriods serinf mmfup mfup

*LOCF (not strictly necessary)
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode renal 0=.
recode diabetes 0=.
recode lung 0=.
*by studyno: replace bmi = bmi[_n-1] if  bmi ==.
by studyno: replace age = age[_n-1] if age ==. 
by studyno: replace steroids = steroids[_n-1] if steroids ==. 
by studyno: replace renal = renal[_n-1] if renal ==. 
by studyno: replace diabetes = diabetes[_n-1] if diabetes ==. 
by studyno: replace lung = lung[_n-1] if lung ==. 
recode renal  .=0
recode diabetes .=0
recode lung .=0
recode steroids .=0
by studyno: replace previous_dmards = previous_dmards[_n-1] if  previous_dmards ==.
by studyno: replace disdur = disdur[_n-1] if  disdur ==.
by studyno: replace smoke = smoke[_n-1] if smoke ==.
drop on_steroid_at_baseline Oralsteriods
replace trtment = firsttreat

sort studyno pyears 
by studyno: egen event = max(serinf)
by studyno: egen timevent = max(mfup) if event == 0
replace timevent = pyears if serinf == 1
by studyno: egen timeventnew = min(timevent) if event == 1
replace timevent = timeventnew if event == 1
drop timeventnew

drop if pyears > timevent
drop if timevent <=0
count if serinf == 1
*egen groupid = group(studyno)

gen expcap = 1 if timevent > mmfup
recode expcap .=0
gen LMvadtime = timevent if expcap == 0
replace LMvadtime = mmfup if expcap == 1
drop mmfup
replace event = 0 if (event == 1 & timevent > 4)
replace timevent = 4 if timevent > 4

save "long_JM2.dta", replace
saveold ".\long_JM.dta", replace version(12)

drop if ovmean ==. 
by studyno: replace dascore = dascore[_n-1] if  dascore ==.
drop studyno
save "JMlong_haq2.dta", replace
saveold ".\JMlong_haq.dta", replace version(12)

use "long_JM2.dta", clear
drop if dascore ==.
by studyno: replace ovmean = ovmean[_n-1] if  ovmean ==.
drop studyno
save "JMlong_das2.dta", replace
saveold ".\JMlong_das.dta", replace version(12)

use "long_JM2.dta", clear
by studyno: gen seq = _n
keep if seq == 1
drop studyno
saveold ".\JMsurv.dta", replace version(12)

log close
