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

log using "$logdir/Dofile_BaselineData.log", text replace

use "tbl_baseline_20181130.dta", clear
keep studyno registdt inputdt biodate todday todmon todyear age height1 weight1 pgen pribio firstbio disdur dascore datend daswol daesr dacrp daglob acrrpos renal copd asthma syspulm newdmard casecont biologic indication diabetes haq dmtherapy ra smoke 

*[[ Exclusion Criterion 1 - Flow Diagram 1 ]]
keep if indication == 1

*[[ Exclusion Criterion 2 - Flow Diagram 1 ]]
keep if firstbio == 1 & (pribio == 1 | pribio == 2 | pribio == 4 | pribio == 1392 | pribio == 1561 | pribio == 10037 | pribio == 10038 | pribio == 10060 | pribio == 10061)

rename firstbio firstbioB
sort studyno
save "tbl_serinf_baseline.dta",replace

**One big file - Treatment history check, and high dose of steroids calculated
use "one_big_file_20181130_for company reports.dta", clear

**Assumptions made about whether individuals are on a high steroid dose or not (Taken from Diederik's analysis on the RABBIT score) 
gen highGCdose=1 if on_steroid_at_baseline==1 & diabetes==0 & smoke==0
recode  highGCdose .=0 if on_steroid_at_baseline==0| diabetes==1 | smoke==1| smoke==2

**Total number of anti-TNFs or anti-TNFs biosimilars and number of treatment switches for individuals

gen sumbio = had1+had2+had3+had4+had5+had1390+had1391+had1392+had1394+had1561+had10037+had10038+had10060+had10061+had10159
label variable sumbio `"Number of biologics/biosimilars received"' 
gen switchbiono=sumbio-1
recode switchbiono -1=0
label variable switchbiono `"How many biologic switches has the patient had?"'
gen switchedbio=1 if switchbiono >0 & switchbiono !=.
recode switchedbio .=0
label variable switchedbio `"Has the patient ever switched a biologic?"' 
bys studyno: egen switcher=max(switchedbio)
label variable switcher `"Patients who have ever switched a biologic"'
gen sumtnf=had1
replace sumtnf=sumtnf+1 if had2==1
replace sumtnf=sumtnf+1 if had4==1
replace sumtnf=sumtnf+1 if had1392==1
replace sumtnf=sumtnf+1 if had1561==1
replace sumtnf=sumtnf+1 if had10037==1
replace sumtnf=sumtnf+1 if had10038==1
replace sumtnf=sumtnf+1 if had10060==1
replace sumtnf=sumtnf+1 if had10061==1
label variable sumtnf `"Number of anti-TNF/anti-TNF biosimilars received"'
gen switchtnfno=sumtnf-1 
recode switchtnfno -1=0
label variable switchtnfno `"How many anti-TNF switches has the patient had?"'
 gen switchedtnf=1 if switchtnfno >0 & switchtnfno !=.
recode switchedtnf .=0
label variable switchedtnf `"Has ever switched from TNF to another anti-TNF"'
bys studyno: egen switchertnf=max(switchedtnf)
label variable switchertnf `"Ever switched to another anti-TNF"'

**Treatment start date, steroid treatment, treatment and treatment past to check this is their first anti-TNF/biosimilar treatment

keep studyno start highGCdose firstbio pribio previous_bio indication sumbio switchbiono switcher sumtnf switchtnfno switchedtnf switchertnf on_steroid_at_baseline orig_studyno previous_dmards
keep if firstbio == 1 & (pribio == 1 | pribio == 2 | pribio == 4 | pribio == 1392 | pribio == 1561 | pribio == 10037 | pribio == 10038 | pribio == 10060 | pribio == 10061)
keep if indication == 1
sort studyno start
by studyno: keep if _n==1 
save "baseline_big_file.dta", replace

**Extracting baseline HAQ data

use "tbl_patientbaselinefup_20181130.dta", clear
sort studyno Fupno 
keep studyno Fupno ovmean HAQbasedate inputdt
rename inputdt inputdtHAQ
keep if Fupno == 0
save "baseline_haq.dta", replace

**Merging all three datasets together

use "tbl_serinf_baseline.dta", clear
merge 1:1 studyno using "baseline_haq.dta"
drop if _merge == 2 //Study ID numbers in the baseline dataset are those we are interested in (first biologic, anti-tnf/biosimilar)
drop _merge
merge 1:1 studyno using "baseline_big_file.dta"
drop if _merge == 2 //same as above
drop _merge

**Observing whether we need to investigate patient history as well as update predictions over time

gen past = 1 if studyno != orig_studyno
recode past .=0

**Asthma and COPD at baseline are collected (different for comorbidites during follow-up)

gen lung_baseline = 1 if (asthma ==1 | copd ==1 | syspulm ==1)
recode lung_baseline .=0
order studyno Fupno start registdt todday todmon todyear inputdt HAQbasedate inputdtHAQ biodate pgen age height1 weight1 pribio previous_dmards disdur dascore ovmean on_steroid_at_baseline acrrpos datend daswol daesr dacrp daglob firstbio firstbioB

lab var start "Baseline treatment start date"
lab var biodate "Date started on biologic"
lab var registdt "Registration date - baseline data date (except HAQ)"
lab var orig_studyno "first study id in BSRBR"
lab var previous_dmards "number of dmards before follow-up date"
lab var highGCdose "on estimated high steroid dose"
lab var Fupno "Follow-up number"
lab var inputdt "Date baseline data was entered into database"
lab var inputdtHAQ "Date baseline HAQ was entered into database"
lab var ovmean "HAQ score"
lab var pribio "Resitered biologic drug ID from baseline dataset"
lab var todday "Day of registration date"
lab var todmon "Month of registration date"
lab var todyear "Year of registration date"
lab var pgen "Patient gender, 1 = female"
lab var age "Patient age at follow-up date"
lab var ra "Patient has RA"
lab var disdur "Disease duration (years)"
lab var biologic "Starting biologic treatment"
lab var asthma "Patient has asthma"
lab var copd "Patient diagnosed with Chronic Obstructive Pulmonary Disease"
lab var renal "Patient diagnosed with Chronic Kidney Disease"
lab var diabetes "Patient has diabetes"
lab var dmtherapy "How is patient's diabetes treated, 1 insulin, 2 tablet controlled, 3 diet controlled, 9 missing"
lab var smoke "Smoking status 0 current smoker, 1, past smoker, 2 never smoked"
lab var weight1 "Patient weight"
lab var height1 "Patient height"
lab var dascore "DAS28 score"
lab var acrrpos "Rheumatoid factor positive?"
lab var datend "28-Tender joint count"
lab var daswol "28 - Swollen joint count"
lab var daesr "ESR count"
lab var dacrp "CRP count" 
lab var haq "Patient has a HAQ recording (yes = 1)"
lab var past "We have past information on a previous studyno for this patient"
lab var lung_baseline "Patient has COPD or asthma"
lab var daglob "Physician's global assessment of disease activity"
lab var firstbioB "is this their first biologic? as recorded in baseline dataset"
lab var HAQbasedate "Date baseline HAQ score was measured"
lab var syspulm "Patient ever had pulmonary fibrosis?"

**Cleaning up missingness
**Assuming missing for "Not required" and for "don't know", for comorbidities, "don't know" assumes patient does not have the condition

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
recode daglob 888=.
recode daglob 999=.
recode dascore 999=.
recode dascore 888=.
recode todday 99=.
recode todday 88=.
recode todmon 99=.
recode todmon 88=.
recode todyear 9999=.
recode todyear 8888=.
replace ovmean = 3 if ovmean > 3
recode diabetes 2=.
recode diabetes 9=.

**Adjusting height and weight where they have been entered incorrectly- same assumptions as project last year

replace height1 = . if height1== 9999 | height1 == 8888 | height1 == 888 | height1 == 999
replace weight1 = . if weight1 == 9999 | weight1 == 8888
replace weight1 = . if (weight1 > 300)
replace height1 = .  if height1 < 100 & weight1 < 100
gen tochange = 0
replace tochange = 1 if height1 < 100 & weight1 > 100
gen weight2 = weight1
replace weight1 = height1 if (tochange == 1)
replace height1 = weight2 if (tochange == 1)
drop weight2 tochange
gen heightmetres = height1 / 100 if (height1 != .)
gen bmi = weight1 / (heightmetres^2) if weight1 != . & height1 != .
replace bmi =. if bmi > 100
save "baseline_serinf_final.dta", replace


**Solving baseline date problem - baseline = treatment start date

use "baseline_serinf_final.dta", clear
gen registdate = dofc(registdt)
format registdate %td
gen treatstart = dofc(biodate)
format treatstart %td
gen HAQdate = dofc(HAQbasedate)
format HAQdate %td
replace HAQdate = registdate if (HAQdate ==. & ovmean !=.)

gen error = treatstart - registdate
lab var error "Time(days) discrepancy between study start and treatment start"
hist error, bin(10) saving($projdir\Output\hist1_baseline, replace)
di "{stata graph use hist1_baseline:Display hist1_baseline graph}"
gen error_two = treatstart - HAQdate
lab var error_two "Time(days) discrepancy between baseline HAQ and treatment start"
hist error_two, bin(10) saving($projdir\Output\hist2_baseline, replace)
di "{stata graph use hist2_baseline:Display hist2_baseline graph}"
scatter error error_two, msize(vsmall) xtitle("HAQ measurement time(days before treatment started)") ytitle("Registration time (days before treatment started)") saving($projdir\Output\scatter1_baseline, replace)
di "{stata graph use scatter1_baseline:Display scatter1_baseline graph}"

*[[ Exclusion Criterion 3 - Flow Chart 1 ]] 

drop if (error > 180 | error < -180)
scatter error error_two, msize(vsmall) xtitle("HAQ measurement time(days before treatment started)") ytitle("Registration time (days before treatment started)") saving($projdir\Output\scatter2_baseline, replace)
di "{stata graph use scatter2_baseline:Display scatter2_baseline graph}"
gen HAQnew = ovmean if (error_two > 180 | error_two < -180)
replace ovmean =. if (error_two > 180 | error_two < -180)

preserve
keep studyno HAQdate HAQnew error_two
keep if HAQnew !=.
keep if error_two < 0
drop error_two
save "first_haq.dta", replace
restore

drop HAQnew HAQdate indication newdmard casecont error error_two registdate biodate inputdtHAQ HAQbasedate inputdt registdt start todday todmon todyear
recode Fupno .=0
drop firstbio
order studyno Fupno treatstart
replace dascore = . if dascore > 10 //2changes

*[[ Exclusion criterion 4 ]]
drop if sumbio > 1 | sumtnf > 1

*[[ Exclusion criterion 5 ]]
drop if pribio == 10037 | pribio ==  10038

save "baseline_serinf_final.dta", replace

**Erase temporary datasets
erase tbl_serinf_baseline.dta 
erase baseline_big_file.dta 
erase baseline_haq.dta

//baseline data (full dataset)

sum bmi, det
bys pgen: sum bmi, det
bys pgen: sum age, det
tab pgen
sum age, det
sum disdur, det
bys pgen: sum disdur, det
sum dascore, det
bys pgen: sum dascore, det
sum ovmean, det
bys pgen: sum ovmean, det
bys pgen: tab acrrpos
tab smoke
bys pgen: tab smoke
tab on_steroid_at_baseline
bys pgen: tab on_steroid_at_baseline
tab highGCdose
bys pgen: tab highGCdose
sum previous_dmards, det
bys pgen: sum previous_dmards, det
tab lung_baseline
tab diabetes
bys pgen: tab diabetes
bys pgen: tab lung_baseline
tab renal
bys pgen: tab renal
tab past
bys pgen: tab past
sum sumbio, det
bys pgen: sum sumbio, det
bys pgen: tab switchedtnf
bys pgen: tab switcher
bys pgen: sum switchbiono, det
bys pgen: sum switchtnfno, det
sum sumtnf, det
bys pgen: sum sumtnf, det

log close
