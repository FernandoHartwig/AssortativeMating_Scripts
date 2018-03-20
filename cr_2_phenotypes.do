//Neil Davies 29/01/18
//This cleans the ALSPAC phenotype data
/* What we need:
1. Height
2. Education

for mum, dad, offspring

PCs, birth year and month
*/

//Partner's height and weight
use aln paw009 paw002 using "pa_3c.dta",clear
replace paw002=. if paw002<0
replace paw009=. if paw009<0
replace paw009=paw009/39.3701
//Convert height to m
save "workingdata/partnerheight",replace

//Partner's education
use "pb_3b.dta",clear

//Clean years of full time education
//Use Okbay method for defining years of education
/*
1 == College or University degree 				ISCED=5 (20 years)
5 == NVQ or HND or HNC or equivalent			ISCED=5 (19 years)
6 == Other prof. qual. eg: nursing, teaching 	ISCED=4 (15 years)
2 == A levels/AS levels or equivalent			ISCED=3 (13 years)
3 == O levels/GCSEs or equivalent 				ISCED=2 (10 years)
4 == CSEs or equivalent							ISCED=2 (10 years)
7 == None of the above 							ISCED=1 (7 years)
8 == Prefer not to answer"						EXCLUDED

"What is the highest level of schooling you have completed?
1) None      ISCED 0
2) CSE only  ISCED 2
3) O-levels  ISCED 3
4) A-levels  ISCED 4
5) University degree and vocational training ISCED 5	
*/
//Father's education
gen father_eduyears=.

replace father_eduyears=20 if pb321==1
replace father_eduyears=19 if (pb319==1|pb318==1) & father_eduyear!=20
replace father_eduyears=15 if (pb320==1|pb315==1|pb316==1) & father_eduyears==.
replace father_eduyears=13 if pb312==1& father_eduyears==.
replace father_eduyears=10 if (pb310==1|pb311==1)& father_eduyears==.
replace father_eduyears=7 if pb322==1 & father_eduyears==.

//Replace eduyears using the recoded variable if missing:
replace father_eduyears=10 if (pb325a==1| pb325a==2) & father_eduyears==.

//Father's report of mother's education
gen f_mother_eduyears=.
replace f_mother_eduyears=20 if pb337==1
replace f_mother_eduyears=19 if (pb334==1|pb335==1) & f_mother_eduyears!=20
replace f_mother_eduyears=15 if (pb336==1|pb331==1|pb332==1) & f_mother_eduyears==.
replace f_mother_eduyears=13 if pb328==1& f_mother_eduyears==.
replace f_mother_eduyears=10 if (pb326==1|pb327==1)& f_mother_eduyears==.
replace f_mother_eduyears=7 if pb338==1 & f_mother_eduyears==.

//Replace eduyears using the recoded variable if missing:
replace f_mother_eduyears=10 if (pb342a==1| pb342a==2) & f_mother_eduyears==.

keep aln f_mother_eduyears father_eduyears
compress
save "workingdata/partnereduc",replace

//Get mothers height
use aln dw021 dw002 using "d_3c.dta",clear
replace dw021=. if dw021<0
replace dw002=. if dw002<0
compress
save "workingdata/motherheight",replace

//Get mothers education
use "c_7d.dta",clear

//Mother's education
gen mother_eduyears=20 if c641==1
replace mother_eduyears=19 if (c638==1|c639==1) & mother_eduyears!=20
replace mother_eduyears=15 if (c635==1|c636==1|c640==1) & mother_eduyears==.
replace mother_eduyears=13 if c632==1& mother_eduyears==.
replace mother_eduyears=10 if (c631==1|c630==1)& mother_eduyears==.
replace mother_eduyears=7 if c642==1 & mother_eduyears==.

//Replace eduyears using the recoded variable if missing:
replace mother_eduyears=10 if (c645a==1|c645a==2) & mother_eduyears==.

//Mother's report of father's education
gen m_father_eduyears=.
replace m_father_eduyears=20 if c661==1
replace m_father_eduyears=19 if (c658==1|c659==1) & m_father_eduyears!=20
replace m_father_eduyears=15 if (c655==1|c656==1|c660==1) & m_father_eduyears==.
replace m_father_eduyears=13 if c652==1& m_father_eduyears==.
replace m_father_eduyears=10 if (c651==1|c650==1)& m_father_eduyears==.
replace m_father_eduyears=7 if c662==1 & m_father_eduyears==.

//Replace eduyears using the recoded variable if missing:
replace m_father_eduyears=10 if (c666a==1| c666a==2) & m_father_eduyears==.

keep aln mother_eduyear m_father_eduyears
compress
save "workingdata/mothereduc",replace

joinby aln using "workingdata/partnereduc",unmatched(both)
drop _m
joinby aln using "workingdata/partnerheight",unmatched(both)
drop _m
joinby aln using "workingdata/motherheight",unmatched(both)
drop _m

//Gen combined education measures

gen c_mother_eduyear=mother_eduyears
replace c_mother_eduyear=f_mother_eduyears if c_mother_eduyear==.
gen c_father_eduyear=father_eduyears
replace c_father_eduyear=m_father_eduyears if c_father_eduyear==.
compress
save "workingdata/combined",replace

//Get the kids education and height
//Height
use aln qlet f7ms010 using "f07_3d.dta",clear
replace f7ms010=. if f7ms010<0
compress
save "workingdata/childheight",replace

joinby aln using "workingdata/combined",unmatched(both)
compress
save "workingdata/phenotypes",replace


