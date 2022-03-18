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

log using "$logdir/Dofile_BigData.log", text replace


**BASELINE FOUNDATION

use "baseline_serinf_final.dta", clear
rename treatstart date
rename pgen gender
drop orig_studyno sumbio switchbiono switcher sumtnf switchtnfno switchedtnf switchertnf past
gen baseline = 1
save "serinf_full.dta", replace

**FOLLOW-UP CLINICAL INFORMATION

use "tbl_fup_20181130.dta", clear
keep studyno Fupno formdate dasdat inputdt weight datend  daswol daesr dacrp dastot daglob Oralsteriods steroids
rename weight weight1
rename dastot dascore
//Correcting mistyped dates (date back to 1900)
replace formdate = dasdat if (studyno == "0020968" & Fupno == 2)
replace formdate = dasdat if (studyno == "0020904" & Fupno == 1)
replace formdate = inputdt if (studyno == "0021540" & Fupno == 1)
save "fup_clinical.dta", replace
//obtaining follow-up information on ONLY the individuals we have in our baseline dataset
use "wanted_IDs.dta", clear
merge 1:m studyno using "fup_clinical.dta"
drop if _merge == 2 | _merge == 1
drop _merge
order studyno Fupno formdate dasdat inputdt weight1 dascore datend daswol daesr dacrp daglob
gen formdate_n = dofc(formdate)
format formdate_n %td
sort studyno formdate_n 
save "fup_clinical.dta", replace
//Add imputation dates
use "tbl_MailingTracking", clear
keep studyno mailingid datesent datereturned datereminder
rename mailingid Fupno
format datesent %td
format datereturned %td 
format datereminder %td
sort studyno Fupno
bysort studyno Fupno: gen seq=_n
drop if seq == 2
drop seq
save "input_dates.dta", replace
use "fup_clinical.dta", clear
merge 1:1 studyno Fupno using "input_dates.dta"
drop if _merge == 2
drop _merge
save "fup_clinical.dta", replace
//Date imputation (Order of priority: FORMDATE DASDATE RETDATE REMDATE DATESENT)
use "fup_clinical.dta", clear
gen dasdate = dofc(dasdat)
format dasdate %td
gen date = formdate_n
format date %td
replace date = dasdate if date ==. 
replace date = datereturned if date ==. 
replace date = datereminder if date ==.
replace date = datesent if date ==. 
drop inputdt dasdat datereturned datereminder datesent formdate_n formdate dasdate
order studyno Fupno date
save "fup_clinical.dta", replace
*Add follow-up clinical data to the baseline data
use "serinf_full.dta", clear
append using "fup_clinical.dta", nolabel generate(filenum)
sort studyno date
tab filenum
drop filenum
order studyno gender date Fupno dascore weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace

**FOLLOW-UP PATIENT INFORMATION (FORMDATE PRETDATE REMDATE DATESENT)

use "tbl_patientbaselinefup_20181130.dta", clear
keep studyno Fupno HAQbasedate ovmean cursmokefu currfrmcmp inputdt pretdate
replace currfrmcmp = HAQbasedate if currfrmcmp ==.
rename currfrmcmp formdate
drop HAQbasedate
gen formdate_n = dofc(formdate)
format formdate_n %td
gen pretdate_n = dofc(pretdate)
format pretdate_n %td
sort studyno Fupno
merge 1:1 studyno Fupno using "input_dates.dta"
drop if _merge == 2
drop formdate pretdate inputdt _merge
sort studyno formdate_n
drop if Fupno == 0
//mis-typed dates in cohort
replace formdate_n = pretdate_n if (studyno == "0019137" & Fupno == 3)
replace formdate_n = pretdate_n if (studyno == "0019076" & Fupno == 3)
replace formdate_n = pretdate_n if (studyno == "0017676" & Fupno == 5)
replace formdate_n = pretdate_n if (studyno == "0019879" & Fupno == 5)
rename formdate_n date
save "fup_patient.dta", replace
//only for individuals we are interested in first time biologic (anti-TNF/biosimilar)
use "wanted_IDs.dta", clear
merge 1:m studyno using "fup_patient.dta"
drop if _merge == 2 | _merge == 1
drop _merge
sort studyno Fupno
recode ovmean 9=.
recode ovmean 8=.
recode cursmokefu 9=.
recode cursmokefu 8=. 
drop if (ovmean ==. & cursmokefu ==.)
replace date = pretdate_n if date ==.
replace date = datereturned if date ==.
replace date = datereminder if date ==.
replace date = datesent if date ==. 
drop pretdate_n datereturned datereminder datesent
save "fup_patient.dta", replace
*Add patient reported follow-up to longitudinal cohort
use "serinf_full.dta", clear
append using "fup_patient.dta", nolabel generate(filenum)
sort studyno date
tab filenum
drop filenum
order studyno gender date Fupno ovmean cursmokefu dascore pribio weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace

** ADDING HAQ MEASUREMENTS TAKEN >6 MONTHS AFTER TREATMENT STARTED

use "serinf_full.dta", clear
append using "first_haq.dta", nolabel generate(filenum)
sort studyno date
tab filenum
drop filenum
order studyno gender date Fupno ovmean cursmokefu dascore pribio weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace

**FOLLOW-UP COMORBIDITIES/ SERIOUS INFECTIONS / DEATHS

use "adverse_events_20181130.dta", clear
keep studyno date serinf death diedate meddraPT meddraHLT meddraHLGT inactive  tbinf
keep if (serinf == 1 | tbinf == 1 | death == 1 | meddraHLT == 242 | meddraHLGT == 66 | meddraHLGT == 242 | meddraHLGT == 191 | (meddraHLGT == 34 & meddraPT != 3594) | (meddraHLT == 705 & meddraPT != 268))
drop if inactive == 1
drop inactive
save "fup_adverse.dta", replace
use "wanted_IDs.dta", clear
merge 1:m studyno using "fup_adverse.dta"
save "fup_adverse.dta", replace
//-------------------------
keep if _merge == 2
keep studyno
save "no_ae_fup.dta", replace
//---------------------------
use "fup_adverse.dta", clear
drop if _merge == 1 | _merge == 2
drop _merge
drop diedate
save "fup_adverse.dta", replace 
use "fup_adverse.dta", clear
keep if (serinf == 1 | tbinf == 1)
save "serinf_events.dta", replace
use "fup_adverse.dta", clear
keep if death == 1
drop meddraPT serinf meddraHLT meddraHLGT
save "death_events.dta", replace
use "fup_adverse.dta", clear
keep if ((meddraHLGT == 34 & meddraPT != 3594) | (meddraHLT == 705 & meddraPT != 268))
gen lung = 1
drop serinf death
save "LD_events.dta", replace
use "fup_adverse.dta", clear
keep if (meddraHLGT == 242 | meddraHLGT == 191) 
gen renal = 1
drop serinf death
save "renal_events.dta", replace
use "fup_adverse.dta", clear
keep if (meddraHLT == 242 | meddraHLGT == 66)
gen diabetes = 1
drop death serinf
save "diabetes_events.dta", replace
//add to longitudinal dataset
use "serinf_full.dta", clear
append using "serinf_events.dta", nolabel generate(filenum)
sort studyno date
drop filenum
order studyno gender date Fupno serinf ovmean cursmokefu dascore pribio weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace
use "serinf_full.dta", clear
append using "death_events.dta", nolabel generate(filenum)
sort studyno date
drop filenum
order studyno gender date Fupno death serinf ovmean cursmokefu dascore pribio weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace
use "serinf_full.dta", clear
rename lung_baseline lung
append using "LD_events.dta", nolabel generate(filenum)
sort studyno date
drop filenum
order studyno gender date Fupno death lung serinf ovmean cursmokefu dascore pribio weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace
use "serinf_full.dta", clear
append using "renal_events.dta", nolabel generate(filenum)
sort studyno date
drop filenum
order studyno gender date Fupno death renal lung serinf ovmean cursmokefu dascore pribio weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace
use "serinf_full.dta", clear
append using "diabetes_events.dta", nolabel generate(filenum)
sort studyno date
drop filenum
order studyno gender date Fupno death diabetes renal lung serinf ovmean cursmokefu dascore pribio weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace


**FOLLOW-UP TREATMENT SWITCHING

use "cleaned_spells_20181130.dta", clear
format start %td
format stop %td
keep studyno start stop bionum
keep studyno start bionum
save "fup_treatments.dta", replace
use "wanted_IDs.dta", clear
merge 1:m studyno using "fup_treatments.dta"
drop if _merge == 1 | _merge == 2
drop _merge
sort studyno start
rename start date
rename bionum trtment
save "fup_treatments.dta", replace
use "serinf_full.dta", clear
append using "fup_treatments.dta", nolabel generate(filenum)
sort studyno date
tab filenum
drop filenum
order studyno gender date Fupno trtment ovmean cursmokefu dascore pribio weight1 datend daswol daesr dacrp daglob
save "serinf_full.dta", replace


**CLEANING AND FOLLOW-UP DAYS/YEARS

use "serinf_full.dta", clear
by studyno: egen pgen = max(gender)
by studyno: egen firsttreat = max(pribio)
order studyno pgen date Fupno trtment firsttreat
sort studyno date
bys studyno: gen pdays = date[_n] - date[1]
gen pyears = pdays/365.25
order studyno pgen date pyears Fupno trtment firsttreat
drop pribio gender
save "serinf_full.dta", replace
sort studyno date
replace trtment = firsttreat if ((trtment ==. | trtment == 9) & pyears == 0)
replace smoke = cursmokefu if smoke ==.
drop cursmokefu
recode renal .=0
recode diabetes .=0
recode diabetes 2=0
recode diabetes 9=0
recode serinf .=0
recode lung .=0
recode death .=0
recode disdur 9999=.
recode dascore 99.99=.
recode ovmean 9=.
recode acrrpos 8=.
recode acrrpos 9=.
recode acrrpos 2=.
recode smoke 8=.
recode smoke 9=.
recode smoke 6=.
recode renal 2=0
recode renal 9=.
recode datend 88=.
recode datend 99=.
recode daswol 88=.
recode daswol 99=.
recode daesr 999=.
recode daesr 888=.
recode dacrp 888=.
recode dacrp 999=.
recode dacrp 999.9=.
recode daglob 888=.
recode daglob 999=.
recode dascore 999=.
recode dascore 888=.
recode weight1 9999=.
recode weight1 8888=.
replace weight1 =. if weight1 > 300
recode steroids 9=.
recode Oralsteriods 8=.
recode Oralsteriods 9=.
recode Oralsteriods 88=.
recode Oralsteriods 99=.
replace datend =. if datend > 28
replace daswol =. if daswol > 28
replace daglob =. if daglob > 100
replace dascore =. if dascore > 10
lab var studyno "Study ID number"
lab var pgen "Patient Gender (female = 1)"
lab var trtment "Treatment switches"
lab var firsttreat "Firstline treatment"
lab var lung "Chronic lung disease onset"
drop HAQdate HAQnew
lab var age "Age patient entered the cohort"
drop haq previous_bio
drop if date ==. 
drop pyears pdays
save "serinf_full.dta", replace

*MERGING NEW STUDY IDS WITH OLD STUDY IDS

use "serinf_full.dta", clear
*Solve missing orig_studyno
replace orig_studyno = studyno if orig_studyno == ""
sort orig_studyno date
*Merge data for the same orig_studyno
gen updatt = 1 if studyno != orig_studyno
replace studyno = orig_studyno if updatt == 1
sort studyno date
drop updatt orig_studyno

**CHECK || bringing any follow-up data collected before baseline to determined baseline point (assuming no more than one?)

sort studyno date
gen basedate = date if baseline == 1
by studyno: egen datebase = max(basedate)
format datebase %td
drop basedate
by studyno: replace date = datebase if date < datebase
*62 changes
save "serinf_full.dta", replace

**MERGING DATA ON SAME DATE

use "serinf_full.dta", clear
sort studyno date
by studyno: egen gender = max(pgen)
by studyno: egen pribio = max(firsttreat)
replace pgen = gender
replace firsttreat = pribio
drop gender pribio Fupno
order studyno pgen date firsttreat trtment ovmean dascore datend daswol daesr dacrp daglob height1 heightmetres bmi serinf tbinf death renal lung asthma copd diabetes Oralsteriods smoke
drop ra firstbioB 
drop dmtherapy baseline 

foreach var of varlist ovmean-bmi{
by studyno date: egen temp`var' = mean(`var')
}
foreach var of varlist serinf-Oralsteriods{
by studyno date: egen temp`var' = max(`var')
}
by studyno date: egen tempsmoke = min(smoke)
by studyno date: egen temptrtment = min(trtment)
foreach var of varlist weight1-disdur{
by studyno date: egen temp`var' = mean(`var')
}
foreach var of varlist on_steroid_at_baseline-steroids{
by studyno date: egen temp`var' = max(`var')
}
foreach var of varlist trtment-steroids{
replace `var' = temp`var'
}

bysort studyno: carryforward trtment, replace
bysort studyno: carryforward firsttreat, replace

*by studyno date: gen countt = _n
*keep if ((count[_n] == 1 & count[_n+1] == 2) | (count[_n] == 2))

by studyno date: keep if _n == 1

foreach var of varlist ovmean-steroids{
drop temp`var'
}
sort studyno date
save "serinf_full.dta", replace

use "serinf_full.dta", clear
sort studyno date
gen event = 1 if serinf == 1 | tbinf == 1
replace event = 0 if event ==.
drop serinf tbinf
rename event serinf

sort studyno date

bys studyno: gen pdays = date[_n] - date[1]
gen pyears = pdays/365.25
order studyno pgen pyears trtment firsttreat
sort studyno pyears
save "serinf_full.dta", replace

erase fup_clinical.dta
erase input_dates.dta
erase fup_patient.dta
erase fup_adverse.dta
erase death_events.dta
erase LD_events.dta
erase serinf_events.dta
erase renal_events.dta
erase diabetes_events.dta
erase fup_treatments.dta

*by studyno: keep if _n == 1

*Exposed-base longitudinal dataset

use "serinf_full.dta", clear
sort studyno pyears
recode serinf .=0
keep studyno pgen pyears pdays age firsttreat trtment on_steroid_at_baseline ovmean dascore bmi weight1 heightmetres death renal lung diabetes smoke previous_dmards disdur steroids Oralsteriods serinf
by studyno: egen heightm = max(heightmetres)
replace bmi = weight1 / (heightm^2) if (heightm !=. & weight1 !=.)
egen groupid = group(studyno)
sum groupid, det
drop groupid

*incorporate info about trtment switch/breaks
merge 1:1 studyno pyears using "trt_breaks.dta"
sort studyno pyears
drop if _merge == 2
drop _merge

*restrict to exposed follow-up
gen switch = 1 if largebreak == 1 
gen switch2 = 1 if firsttreat != trtment
replace switch = 0 if switch ==.
replace switch2 = 0 if switch2 ==.
tab switch
tab switch2
gen maxfup = pyears+0.25 if (switch == 1 | switch2 ==1) 
sum maxfup, det
*expands these values for the subjects where the exposure has been capped
by studyno: egen mfup = min(maxfup)

*Keeping note of exposure cap time for temporal assessmnet datasets (mmfup)
gen mmfup = mfup - 0.25 if mfup !=. 

by studyno: egen maxfup2 = max(pyears)
replace mfup = maxfup2 if mfup ==. 
replace mmfup = mfup if mmfup ==. 
gen dropp = 1 if pyears > mfup
egen groupid = group(studyno) if dropp == 1
sum groupid, det
drop if dropp == 1
drop groupid dropp

*motivating data plots
preserve
egen groupid = group(studyno)
keep if groupid < 50
*set scheme plottig
xtline weight1, t(pyears) i(groupid) overlay legend(off) scheme(s1mono) xtitle("Time since baseline (years)") ytitle("Weight (kg)")
xtline ovmean if pyears < 4, t(pyears) i(groupid) overlay legend(off) scheme(s1mono) xtitle("Time since baseline (years)")
xtline dascore, t(pyears) i(groupid) overlay legend(off) scheme(s1mono) xtitle("Time since baseline (years)")
recode steroids .=0
xtline steroids, t(pyears) i(groupid) overlay legend(off) scheme(s1mono) xtitle("Time since baseline (years)") ylab(0 1) ytitle("On steroids")
*Restrict to 4 years of follow-up
restore

sort studyno pyears
drop if pyears > 4
*Events and censoring
by studyno: egen event = max(serinf)
by studyno: gen timevent2 = pyears if serinf == 1
by studyno: replace timevent2 = timevent2[_n-1] if  timevent2 ==.
gsort studyno -pyears
by studyno: replace timevent2 = timevent2[_n-1] if  timevent2 ==.
sort studyno pyears
bysort studyno: gen timevent = timevent2 if event == 1
replace timevent = mfup if event == 0
replace timevent = 4 if timevent > 4
drop if pyears > timevent

*removed those who died before event/censoring time = data error
sort studyno pyears
by studyno: egen died = max(death)
gen timedeath = pyears if death == 1
by studyno: egen diedtime = max(timedeath)
gen dropp = 1 if diedtime < mfup
egen groupid = group(studyno) if dropp == 1
sum groupid, det
drop if dropp == 1
drop largebreak maxfup maxfup2 timedeath diedtime death groupid dropp

recode steroids 8=0
*restrict to complete-cases at baseline
by studyno: egen mintimevent = min(timevent)
drop if mintimevent <= 0
*long treatment break from baseline
drop if mintimevent == 0.25 & event == 0

egen groupid = group(studyno)
sum groupid, det
by studyno: egen numev = sum(serinf)
replace trtment = firsttreat if trtment != firsttreat
drop if (firsttreat == 10038 | firsttreat == 10037)
save "serinf_full_study_noncomplete.dta", replace
saveold "R:\BSRBR\Analyses\lucy_bull\Serious_infection\R Datasets\long_noncomplete.dta", ///
replace version(12)

use "serinf_full_study_noncomplete.dta", clear

*Complete case from baseline
drop groupid
gen miss = 1 if ((firsttreat ==. | ovmean ==. | dascore ==. | bmi ==. | disdur ==. | previous_dmards ==. | smoke ==. | age ==. ) & (pyears == 0))
by studyno: egen drop = max(miss)
egen groupid = group(studyno) if drop == 1
sum groupid, det
drop if drop == 1
drop drop miss groupid

egen groupid = group(studyno)
sum groupid, det
save "serinf_full_study.dta", replace
saveold "R:\BSRBR\Analyses\lucy_bull\Serious_infection\R Datasets\long_complete.dta", ///
replace version(12)



*MAKING TREATMENT BREAK/SWITCH DATA





**Supplementary dataset = identifying when large breaks of treatment start
use "cleaned_spells_20181130.dta", clear
format start %td
format stop %td
keep studyno start stop bionum
sort studyno start
save "fup_treatments.dta", replace

use "wanted_IDs.dta", clear
merge 1:m studyno using "fup_treatments.dta"

*drop those with no treatment info in cleaned_spells, and those who are not included in baseline dataset
drop if _merge == 1 | _merge == 2
drop _merge
sort studyno start
by studyno: gen breaktime = stop-start if bionum == 9
gen largebreak = 1 if (breaktime >= 90 & breaktime !=.)

replace studyno = orig_studyno if ((studyno != orig_studyno) & (orig_studyno != ""))

drop orig_studyno

sort studyno start
by studyno: gen pdays = start[_n]-start[1]
by studyno: gen pyears = pdays/365.25
keep if largebreak == 1
keep studyno pyears largebreak
order studyno pyears largebreak



save "trt_breaks.dta", replace











****[[ Data Issue 2 ]]
**QUANTIFYING DATE ISSUE - SCROLL DOWN FOR DATASET CREATION

use "tbl_fup_20181130.dta", clear
keep studyno Fupno formdate dasdat inputdt weight datend  daswol daesr dacrp dastot daglob Oralsteriods whichbio steroids
save "fup_serinf_app.dta", replace
use "baseline_serinf_final_two.dta", clear
keep studyno
merge 1:m studyno using "fup_serinf_app.dta"
drop if _merge == 2 | _merge == 1
drop _merge
save "fup_serinf_app.dta", replace 
replace datend =. if datend > 28
replace daswol =. if daswol > 28
recode daesr 999=.
recode daesr 888=.
recode dacrp 999.9=.
recode dacrp 888=.
recode dacrp 999=.
replace daglob =. if daglob > 100
replace dastot =. if dastot > 10
gen formdate_n = dofc(formdate)
format formdate_n %td
sort studyno formdate_n
foreach var of varlist datend-dastot{
by studyno formdate_n: egen std`var' = sd(`var')
sum std`var', det
}

foreach var of varlist datend-dastot{
count if std`var' > 0 & std`var' !=.
}
save "fup_serinf_app.dta", replace 

use "fup_serinf_app.dta", clear 
keep if stddastot > 0 & stddastot !=.
by studyno formdate_n: keep if _n == 1
gen count = 1
total(count)

use "fup_serinf_app.dta", clear 
keep if stddatend > 0 & stddatend !=.
by studyno formdate_n: keep if _n == 1
gen count = 1
total(count)

use "fup_serinf_app.dta", clear 
keep if stddaswol > 0 & stddaswol !=.
by studyno formdate_n: keep if _n == 1
gen count = 1
total(count)

use "fup_serinf_app.dta", clear 
keep if stddaesr > 0 & stddaesr !=.
by studyno formdate_n: keep if _n == 1
gen count = 1
total(count)

use "fup_serinf_app.dta", clear 
keep if stddacrp > 0 & stddacrp !=.
by studyno formdate_n: keep if _n == 1
gen count = 1
total(count)

use "fup_serinf_app.dta", clear 
keep if stddaglob > 0 & stddaglob !=.
by studyno formdate_n: keep if _n == 1
gen count = 1
total(count)

count if formdate_n ==.

log close
