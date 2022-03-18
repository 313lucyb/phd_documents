*Creating data for TDCM model

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

log using "$logdir/Dofile_TDCM2.log", text replace

use "serinf_full_study.dta", clear
keep studyno groupid pgen pdays pyears age firsttreat trtment on_steroid_at_baseline ovmean dascore bmi renal lung diabetes smoke previous_dmards disdur steroids Oralsteriods serinf mfup

*LOCF
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode steroids 8=.
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


*identify those with multiple events
by studyno: gen event = serinf[_n+1] if serinf[_n] < serinf[_n+1]
replace event = 1 if (serinf == 1 & pyears == 0)
recode event .=0
gen fupyears = mfup
drop if fupyears <= 0
by studyno: gen delete = 1 if (serinf == 1 & pyears == fupyears)
by studyno: gen start = pyears 
by studyno: gen end = pyears[_n+1]
replace end = pyears + 0.5 if (end ==. & delete !=1) 
drop if delete == 1 
replace serinf = event
drop event delete
*egen groupid = group(studyno)
count if serinf == 1

recode serinf 0=.
by studyno: gen serinf2 = serinf[_n-1] if serinf ==.
recode serinf .=0
by studyno: replace serinf2 = serinf2[_n-1] if  serinf2 ==.
rename serinf2 prev_sinf
recode prev_sinf .=0

gen interval = end - start
drop if interval == 0
drop interval

save "long_TDCM2.dta", replace
saveold ".\long_TDCM.dta", replace version(12)

log close

