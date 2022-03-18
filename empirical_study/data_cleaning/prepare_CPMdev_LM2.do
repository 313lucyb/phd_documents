*Creating data for LM model

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

log using "$logdir/Dofile_LM2.log", text replace

use "serinf_full_study.dta", clear
keep studyno groupid  pgen pyears age firsttreat trtment on_steroid_at_baseline ovmean dascore bmi renal lung diabetes smoke previous_dmards disdur steroids Oralsteriods serinf mmfup mfup mintimevent

*LOCF
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode renal 0=.
recode diabetes 0=.
recode lung 0=.
by studyno: replace ovmean = ovmean[_n-1] if  ovmean ==.
by studyno: replace dascore = dascore[_n-1] if  dascore ==.
by studyno: replace bmi = bmi[_n-1] if  bmi ==.
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
by studyno: replace smoke = smoke[_n-1] if  smoke ==.
drop on_steroid_at_baseline Oralsteriods
replace trtment = firsttreat

gen expcap = 1 if mintimevent > mmfup
recode expcap .=0
gen LMvadtime = mintimevent if expcap == 0
replace LMvadtime = mmfup if expcap == 1
drop mmfup mintimevent

*events and censoring
gen maxfup = mfup
sort studyno pyears 
by studyno: egen event = max(serinf)
by studyno: egen timevent = max(maxfup) if event == 0
replace timevent = pyears if serinf == 1
by studyno: egen maxfup2 = max(timevent) if event == 1
replace maxfup = maxfup2 if event == 1
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
count if serinf == 1

gen smoke0 = 1 if smoke == 0
gen smoke1 = 1 if smoke == 1
gen smoke2 = 1 if smoke == 2
recode smoke0 .=0
recode smoke1 .=0
recode smoke2 .=0
*egen groupid = group(studyno)
replace timevent = 4 if (event == 0 & timevent > 4)
save "long_LM2.dta", replace
saveold ".\long_LM.dta", replace version(12)


keep studyno groupid event timevent LMvadtime maxfup pyears age
save "TA_JM_events.dta", replace
saveold ".\TA_JM_events.dta", replace version(12)

log close
 
