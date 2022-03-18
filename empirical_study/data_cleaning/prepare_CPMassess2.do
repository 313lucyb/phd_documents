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

log using "$logdir/Dofile_AssessCPM2.log", text replace

use "serinf_full_study.dta", clear
keep studyno groupid pgen pyears age firsttreat trtment on_steroid_at_baseline ovmean dascore bmi died renal lung diabetes smoke previous_dmards disdur steroids Oralsteriods serinf mmfup mfup mintimevent
sort studyno pyears
gen expcap = 1 if mintimevent > mmfup
recode expcap .=0
gen LMvadtime = mintimevent if expcap == 0
replace LMvadtime = mmfup if expcap == 1
drop mmfup mintimevent
save "temporal_base.dta", replace

*~~~~~~~~~~
*Validation dataset: landmark time = 6 months [SIs allowed, no switching, no treatment breaks > 90 days before LT]. 
*~~~~~~~~~~~

use "temporal_base.dta", clear
keep if LMvadtime > 0.5

*Events and censoring

keep if pyears <= 1.5
sort studyno pyears
gen predwindow = 1 if (pyears > 0.5 & pyears <=1.5)
by studyno: egen event = max(serinf) if predwindow == 1
gen timevent = pyears if (serinf == 1 & event == 1)
by studyno: egen eventnew = max(event) 
recode eventnew .=0
by studyno: replace timevent = mfup if eventnew == 0
replace timevent = 1.5 if (timevent > 1.5 & eventnew == 0)
drop predwindow
*first serious infection in prediction window
by studyno: egen timeventnew = min(timevent) if eventnew == 1
replace timevent = timeventnew if eventnew == 1
drop event timeventnew
rename eventnew event
*timefromLT 
gen timefromLT = timevent - 0.5

*LOCF
sort studyno pyears
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode steroids .=0
recode renal 0=.
recode diabetes 0=.
recode lung 0=.
by studyno: replace ovmean = ovmean[_n-1] if  ovmean ==.
by studyno: replace dascore = dascore[_n-1] if  dascore ==.
by studyno: replace bmi = bmi[_n-1] if  bmi ==.
by studyno: replace age = age[_n-1] if age ==. 
by studyno: replace renal = renal[_n-1] if renal ==. 
by studyno: replace diabetes = diabetes[_n-1] if diabetes ==. 
by studyno: replace lung = lung[_n-1] if lung ==. 
recode renal  .=0
recode diabetes .=0
recode lung .=0
by studyno: replace previous_dmards = previous_dmards[_n-1] if  previous_dmards ==.
by studyno: replace disdur = disdur[_n-1] if  disdur ==.
by studyno: replace smoke = smoke[_n-1] if  smoke ==.
drop on_steroid_at_baseline Oralsteriods
gen trt1 = 1 if trtment == 1
gen trt2 = 1 if trtment == 2
gen trt4 = 1 if trtment == 4
gen trt1561 = 1 if trtment == 1561
gen trt10037 = 1 if trtment == 10037
gen trt10038 = 1 if trtment == 10038
gen trt10060 = 1 if trtment == 10060
recode trt1 .=0
recode trt2 .=0
recode trt4 .=0
recode trt1561 .=0
recode trt10037 .=0
recode trt10038 .=0
recode trt10060 .=0
gen smoke0 = 1 if smoke == 0
gen smoke1 = 1 if smoke == 1
gen smoke2 = 1 if smoke == 2
recode smoke0 .=0
recode smoke1 .=0
recode smoke2 .=0
recode serinf 0=.
by studyno: gen serinf2 = serinf[_n-1] if serinf ==.
recode serinf .=0
by studyno: replace serinf2 = serinf2[_n-1] if  serinf2 ==.
rename serinf2 prev_sinf
recode prev_sinf .=0
recode steroids 8=0
*Taking most recent covariates at prediction time
drop if pyears > 0.5
by studyno: gen seq = _n
by studyno: gen max = _N
keep if seq == max
gen updated = 1 if pyears > 0
count if updated == 1
drop seq max updated
count if event == 1
save "validation_6months.dta", replace
saveold ".\validation_6months.dta", ///
replace version(12)

*no censoring in prediction window
use "validation_6months.dta", clear
drop if (timevent < 1.5 & event == 0)
save "validation_6binary.dta", replace
saveold ".\validation_6binary.dta", ///
replace version(12)

*~~~~~~~~~~
*Validation dataset: landmark time = 12 months [SIs allowed, no switching, no treatment breaks > 90 days before LT]. 
*~~~~~~~~~~~

use "temporal_base.dta", clear
keep if LMvadtime > 1

*Events and censoring

keep if pyears <= 2
sort studyno pyears
gen predwindow = 1 if (pyears > 1 & pyears <=2)
by studyno: egen event = max(serinf) if predwindow == 1
gen timevent = pyears if (serinf == 1 & event == 1)
by studyno: egen eventnew = max(event) 
recode eventnew .=0
by studyno: replace timevent = mfup if eventnew == 0
replace timevent = 2 if (timevent > 2 & eventnew == 0)
drop predwindow
*first serious infection in prediction window
by studyno: egen timeventnew = min(timevent) if eventnew == 1
replace timevent = timeventnew if eventnew == 1
drop event timeventnew
rename eventnew event
*timefromLT 
gen timefromLT = timevent - 1

*LOCF
sort studyno pyears
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode steroids .=0
recode renal 0=.
recode diabetes 0=.
recode lung 0=.
by studyno: replace ovmean = ovmean[_n-1] if  ovmean ==.
by studyno: replace dascore = dascore[_n-1] if  dascore ==.
by studyno: replace bmi = bmi[_n-1] if  bmi ==.
by studyno: replace age = age[_n-1] if age ==. 
by studyno: replace renal = renal[_n-1] if renal ==. 
by studyno: replace diabetes = diabetes[_n-1] if diabetes ==. 
by studyno: replace lung = lung[_n-1] if lung ==. 
recode renal  .=0
recode diabetes .=0
recode lung .=0
by studyno: replace previous_dmards = previous_dmards[_n-1] if  previous_dmards ==.
by studyno: replace disdur = disdur[_n-1] if  disdur ==.
by studyno: replace smoke = smoke[_n-1] if  smoke ==.
drop on_steroid_at_baseline Oralsteriods
gen trt1 = 1 if trtment == 1
gen trt2 = 1 if trtment == 2
gen trt4 = 1 if trtment == 4
gen trt1561 = 1 if trtment == 1561
gen trt10037 = 1 if trtment == 10037
gen trt10038 = 1 if trtment == 10038
gen trt10060 = 1 if trtment == 10060
recode trt1 .=0
recode trt2 .=0
recode trt4 .=0
recode trt1561 .=0
recode trt10037 .=0
recode trt10038 .=0
recode trt10060 .=0
gen smoke0 = 1 if smoke == 0
gen smoke1 = 1 if smoke == 1
gen smoke2 = 1 if smoke == 2
recode smoke0 .=0
recode smoke1 .=0
recode smoke2 .=0
recode serinf 0=.
by studyno: gen serinf2 = serinf[_n-1] if serinf ==.
recode serinf .=0
by studyno: replace serinf2 = serinf2[_n-1] if  serinf2 ==.
rename serinf2 prev_sinf
recode prev_sinf .=0

*Taking most recent covariates at prediction time
drop if pyears > 1
by studyno: gen seq = _n
by studyno: gen max = _N
keep if seq == max
gen updated = 1 if pyears > 0.5
count if updated == 1
drop seq max updated
count if event == 1
save "validation_12months.dta", replace
saveold ".\validation_12months.dta", ///
replace version(12)

*no censoring in prediction window
use "validation_12months.dta", clear
drop if (timevent < 2 & event == 0)
save "validation_12binary.dta", replace
saveold ".\validation_12binary.dta", ///
replace version(12)

*~~~~~~~~~~
*Validation dataset: landmark time = 18 months [SIs allowed, no switching, no treatment breaks > 90 days before LT]. 
*~~~~~~~~~~~

use "temporal_base.dta", clear
keep if LMvadtime > 1.5

*Events and censoring
keep if pyears <= 2.5
sort studyno pyears
gen predwindow = 1 if (pyears > 1.5 & pyears <=2.5)
by studyno: egen event = max(serinf) if predwindow == 1
gen timevent = pyears if (serinf == 1 & event == 1)
by studyno: egen eventnew = max(event) 
recode eventnew .=0
by studyno: replace timevent = mfup if eventnew == 0
replace timevent = 2.5 if (timevent > 2.5 & eventnew == 0)
drop predwindow
*first serious infection in prediction window
by studyno: egen timeventnew = min(timevent) if eventnew == 1
replace timevent = timeventnew if eventnew == 1
drop event timeventnew
rename eventnew event
*timefromLT 
gen timefromLT = timevent - 1.5

*LOCF
sort studyno pyears
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode steroids .=0
recode renal 0=.
recode diabetes 0=.
recode lung 0=.
by studyno: replace ovmean = ovmean[_n-1] if  ovmean ==.
by studyno: replace dascore = dascore[_n-1] if  dascore ==.
by studyno: replace bmi = bmi[_n-1] if  bmi ==.
by studyno: replace age = age[_n-1] if age ==. 
by studyno: replace renal = renal[_n-1] if renal ==. 
by studyno: replace diabetes = diabetes[_n-1] if diabetes ==. 
by studyno: replace lung = lung[_n-1] if lung ==. 
recode renal  .=0
recode diabetes .=0
recode lung .=0
by studyno: replace previous_dmards = previous_dmards[_n-1] if  previous_dmards ==.
by studyno: replace disdur = disdur[_n-1] if  disdur ==.
by studyno: replace smoke = smoke[_n-1] if  smoke ==.
drop on_steroid_at_baseline Oralsteriods
gen trt1 = 1 if trtment == 1
gen trt2 = 1 if trtment == 2
gen trt4 = 1 if trtment == 4
gen trt1561 = 1 if trtment == 1561
gen trt10037 = 1 if trtment == 10037
gen trt10038 = 1 if trtment == 10038
gen trt10060 = 1 if trtment == 10060
recode trt1 .=0
recode trt2 .=0
recode trt4 .=0
recode trt1561 .=0
recode trt10037 .=0
recode trt10038 .=0
recode trt10060 .=0
gen smoke0 = 1 if smoke == 0
gen smoke1 = 1 if smoke == 1
gen smoke2 = 1 if smoke == 2
recode smoke0 .=0
recode smoke1 .=0
recode smoke2 .=0
recode serinf 0=.
by studyno: gen serinf2 = serinf[_n-1] if serinf ==.
recode serinf .=0
by studyno: replace serinf2 = serinf2[_n-1] if  serinf2 ==.
rename serinf2 prev_sinf
recode prev_sinf .=0

*Taking most recent covariates at prediction time
drop if pyears > 1.5
by studyno: gen seq = _n
by studyno: gen max = _N
keep if seq == max
gen updated = 1 if pyears > 1
count if updated == 1
drop seq max updated
count if event == 1
save "validation_18months.dta", replace
saveold ".\validation_18months.dta", ///
replace version(12)

*no censoring in prediction window
use "validation_18months.dta", clear
drop if (timevent < 2.5 & event == 0)
save "validation_18binary.dta", replace
saveold ".\validation_18binary.dta", ///
replace version(12)

*~~~~~~~~~~
*Validation dataset: landmark time = 24 months [SIs allowed, no switching, no treatment breaks > 90 days before LT]. 
*~~~~~~~~~~~

use "temporal_base.dta", clear
keep if LMvadtime > 2

*Events and censoring
keep if pyears <= 3
sort studyno pyears
gen predwindow = 1 if (pyears > 2 & pyears <=3)
by studyno: egen event = max(serinf) if predwindow == 1
gen timevent = pyears if (serinf == 1 & event == 1)
by studyno: egen eventnew = max(event) 
recode eventnew .=0
by studyno: replace timevent = mfup if eventnew == 0
replace timevent = 3 if (timevent > 3 & eventnew == 0)
drop predwindow
*first serious infection in prediction window
by studyno: egen timeventnew = min(timevent) if eventnew == 1
replace timevent = timeventnew if eventnew == 1
drop event timeventnew
rename eventnew event
*timefromLT 
gen timefromLT = timevent - 2

*LOCF
sort studyno pyears
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode steroids .=0
recode renal 0=.
recode diabetes 0=.
recode lung 0=.
by studyno: replace ovmean = ovmean[_n-1] if  ovmean ==.
by studyno: replace dascore = dascore[_n-1] if  dascore ==.
by studyno: replace bmi = bmi[_n-1] if  bmi ==.
by studyno: replace age = age[_n-1] if age ==. 
by studyno: replace renal = renal[_n-1] if renal ==. 
by studyno: replace diabetes = diabetes[_n-1] if diabetes ==. 
by studyno: replace lung = lung[_n-1] if lung ==. 
recode renal  .=0
recode diabetes .=0
recode lung .=0
by studyno: replace previous_dmards = previous_dmards[_n-1] if  previous_dmards ==.
by studyno: replace disdur = disdur[_n-1] if  disdur ==.
by studyno: replace smoke = smoke[_n-1] if  smoke ==.
drop on_steroid_at_baseline Oralsteriods
gen trt1 = 1 if trtment == 1
gen trt2 = 1 if trtment == 2
gen trt4 = 1 if trtment == 4
gen trt1561 = 1 if trtment == 1561
gen trt10037 = 1 if trtment == 10037
gen trt10038 = 1 if trtment == 10038
gen trt10060 = 1 if trtment == 10060
recode trt1 .=0
recode trt2 .=0
recode trt4 .=0
recode trt1561 .=0
recode trt10037 .=0
recode trt10038 .=0
recode trt10060 .=0
gen smoke0 = 1 if smoke == 0
gen smoke1 = 1 if smoke == 1
gen smoke2 = 1 if smoke == 2
recode smoke0 .=0
recode smoke1 .=0
recode smoke2 .=0
recode serinf 0=.
by studyno: gen serinf2 = serinf[_n-1] if serinf ==.
recode serinf .=0
by studyno: replace serinf2 = serinf2[_n-1] if  serinf2 ==.
rename serinf2 prev_sinf
recode prev_sinf .=0

*Taking most recent covariates at prediction time
drop if pyears > 2
by studyno: gen seq = _n
by studyno: gen max = _N
keep if seq == max
gen updated = 1 if pyears > 1.5
count if updated == 1
drop seq max updated
count if event == 1
save "validation_24months.dta", replace
saveold ".\validation_24months.dta", ///
replace version(12)

*no censoring in prediction window
use "validation_24months.dta", clear
drop if (timevent < 3 & event == 0)
save "validation_24binary.dta", replace
saveold ".\validation_24binary.dta", ///
replace version(12)

*~~~~~~~~~~
*Validation dataset: landmark time = 30 months [SIs allowed, no switching, no treatment breaks > 90 days before LT]. 
*~~~~~~~~~~~

use "temporal_base.dta", clear
keep if LMvadtime > 2.5

*Events and censoring
keep if pyears <= 3.5
sort studyno pyears
gen predwindow = 1 if (pyears > 2.5 & pyears <=3.5)
by studyno: egen event = max(serinf) if predwindow == 1
gen timevent = pyears if (serinf == 1 & event == 1)
by studyno: egen eventnew = max(event) 
recode eventnew .=0
by studyno: replace timevent = mfup if eventnew == 0
replace timevent = 3.5 if (timevent > 3.5 & eventnew == 0)
drop predwindow
*first serious infection in prediction window
by studyno: egen timeventnew = min(timevent) if eventnew == 1
replace timevent = timeventnew if eventnew == 1
drop event timeventnew
rename eventnew event
*timefromLT 
gen timefromLT = timevent - 2.5

*LOCF
sort studyno pyears
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode steroids .=0
recode renal 0=.
recode diabetes 0=.
recode lung 0=.
by studyno: replace ovmean = ovmean[_n-1] if  ovmean ==.
by studyno: replace dascore = dascore[_n-1] if  dascore ==.
by studyno: replace bmi = bmi[_n-1] if  bmi ==.
by studyno: replace age = age[_n-1] if age ==. 
by studyno: replace renal = renal[_n-1] if renal ==. 
by studyno: replace diabetes = diabetes[_n-1] if diabetes ==. 
by studyno: replace lung = lung[_n-1] if lung ==. 
recode renal  .=0
recode diabetes .=0
recode lung .=0
by studyno: replace previous_dmards = previous_dmards[_n-1] if  previous_dmards ==.
by studyno: replace disdur = disdur[_n-1] if  disdur ==.
by studyno: replace smoke = smoke[_n-1] if  smoke ==.
drop on_steroid_at_baseline Oralsteriods
gen trt1 = 1 if trtment == 1
gen trt2 = 1 if trtment == 2
gen trt4 = 1 if trtment == 4
gen trt1561 = 1 if trtment == 1561
gen trt10037 = 1 if trtment == 10037
gen trt10038 = 1 if trtment == 10038
gen trt10060 = 1 if trtment == 10060
recode trt1 .=0
recode trt2 .=0
recode trt4 .=0
recode trt1561 .=0
recode trt10037 .=0
recode trt10038 .=0
recode trt10060 .=0
gen smoke0 = 1 if smoke == 0
gen smoke1 = 1 if smoke == 1
gen smoke2 = 1 if smoke == 2
recode smoke0 .=0
recode smoke1 .=0
recode smoke2 .=0
recode serinf 0=.
by studyno: gen serinf2 = serinf[_n-1] if serinf ==.
recode serinf .=0
by studyno: replace serinf2 = serinf2[_n-1] if  serinf2 ==.
rename serinf2 prev_sinf
recode prev_sinf .=0

*Taking most recent covariates at prediction time
drop if pyears > 2.5
by studyno: gen seq = _n
by studyno: gen max = _N
keep if seq == max
gen updated = 1 if pyears > 2
count if updated == 1
drop seq max updated
count if event == 1
save "validation_30months.dta", replace
saveold ".\validation_30months.dta", ///
replace version(12)

*no censoring in prediction window
*CHANGE
use "validation_30months.dta", clear
drop if (timevent < 3.5 & event == 0)
save "validation_30binary.dta", replace
saveold ".\validation_30binary.dta", ///
replace version(12)

*~~~~~~~~~~
*Validation dataset: landmark time = 36 months [SIs allowed, no switching, no treatment breaks > 90 days before LT]. 
*~~~~~~~~~~~

use "temporal_base.dta", clear
keep if LMvadtime > 3

*Events and censoring
keep if pyears <= 4
sort studyno pyears
gen predwindow = 1 if (pyears > 3 & pyears <=4)
by studyno: egen event = max(serinf) if predwindow == 1
gen timevent = pyears if (serinf == 1 & event == 1)
by studyno: egen eventnew = max(event) 
recode eventnew .=0
by studyno: replace timevent = mfup if eventnew == 0
replace timevent = 4 if (timevent > 4 & eventnew == 0)
drop predwindow
*first serious infection in prediction window
by studyno: egen timeventnew = min(timevent) if eventnew == 1
replace timevent = timeventnew if eventnew == 1
drop event timeventnew
rename eventnew event
*timefromLT 
gen timefromLT = timevent - 3

*LOCF
sort studyno pyears
recode steroids 0=.
recode on_steroid_at_baseline 0=.
recode Oralsteriods 0=.
replace steroids = on_steroid_at_baseline if on_steroid_at_baseline !=.
replace steroids = Oralsteriods if Oralsteriods !=.
recode steroids .=0
recode renal 0=.
recode diabetes 0=.
recode lung 0=.
by studyno: replace ovmean = ovmean[_n-1] if  ovmean ==.
by studyno: replace dascore = dascore[_n-1] if  dascore ==.
by studyno: replace bmi = bmi[_n-1] if  bmi ==.
by studyno: replace age = age[_n-1] if age ==. 
by studyno: replace renal = renal[_n-1] if renal ==. 
by studyno: replace diabetes = diabetes[_n-1] if diabetes ==. 
by studyno: replace lung = lung[_n-1] if lung ==. 
recode renal  .=0
recode diabetes .=0
recode lung .=0
by studyno: replace previous_dmards = previous_dmards[_n-1] if  previous_dmards ==.
by studyno: replace disdur = disdur[_n-1] if  disdur ==.
by studyno: replace smoke = smoke[_n-1] if  smoke ==.
drop on_steroid_at_baseline Oralsteriods
gen trt1 = 1 if trtment == 1
gen trt2 = 1 if trtment == 2
gen trt4 = 1 if trtment == 4
gen trt1561 = 1 if trtment == 1561
gen trt10037 = 1 if trtment == 10037
gen trt10038 = 1 if trtment == 10038
gen trt10060 = 1 if trtment == 10060
recode trt1 .=0
recode trt2 .=0
recode trt4 .=0
recode trt1561 .=0
recode trt10037 .=0
recode trt10038 .=0
recode trt10060 .=0
gen smoke0 = 1 if smoke == 0
gen smoke1 = 1 if smoke == 1
gen smoke2 = 1 if smoke == 2
recode smoke0 .=0
recode smoke1 .=0
recode smoke2 .=0
recode serinf 0=.
by studyno: gen serinf2 = serinf[_n-1] if serinf ==.
recode serinf .=0
by studyno: replace serinf2 = serinf2[_n-1] if  serinf2 ==.
rename serinf2 prev_sinf
recode prev_sinf .=0

*Taking most recent covariates at prediction time
drop if pyears > 3
by studyno: gen seq = _n
by studyno: gen max = _N
keep if seq == max
gen updated = 1 if pyears > 2.5
count if updated == 1
drop seq max updated
count if event == 1
save "validation_36months.dta", replace
saveold ".\validation_36months.dta", ///
replace version(12)

*no censoring in prediction window
use "validation_36months.dta", clear
drop if (timevent < 4 & event == 0)
save "validation_36binary.dta", replace
saveold ".\validation_36binary.dta", ///
replace version(12)

*~~~~~~~~
* Creating temporal datasets safe to use on CSF
*~~~~~~~~

use "validation_6months.dta", clear
drop studyno
save "validation_6months_CSF.dta", replace
saveold ".\validation_6months_CSF.dta", ///
replace version(12)

use "validation_12months.dta", clear
drop studyno
save "validation_12months_CSF.dta", replace
saveold ".\validation_12months_CSF.dta", ///
replace version(12)

use "validation_18months.dta", clear
drop studyno
save "validation_18months_CSF.dta", replace
saveold ".\validation_18months_CSF.dta", ///
replace version(12)

use "validation_24months.dta", clear
drop studyno
save "validation_24months_CSF.dta", replace
saveold ".\validation_24months_CSF.dta", ///
replace version(12)

use "validation_30months.dta", clear
drop studyno
save "validation_30months_CSF.dta", replace
saveold ".\validation_30months_CSF.dta", ///
replace version(12)

use "validation_36months.dta", clear
drop studyno
save "validation_36months_CSF.dta", replace
saveold ".\validation_36months_CSF.dta", ///
replace version(12)

log close
