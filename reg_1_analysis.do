//Neil Davies 13/02/18
//This runs the regressions for the assortative mating paper.

use "20180208_B2347_Davies_1.dta",clear

rename f7ms010  offspring_height
rename paw009 father_height
rename dw021 mother_height
rename average_points_raw offspring_educ

replace father_eduyears =m_father_eduyears if father_eduyears ==.
replace mother_eduyears =f_mother_eduyears if mother_eduyears ==.

//Summary of data
sum offspring_height offspring_educ mother_height mother_eduyears father_height father_eduyears o_height_score o_ea2_score m_height_score m_ea2_score f_height_score   f_ea2_score if offspring_height!=. & offspring_educ!=., sep(0)

//Phenotypic correlations offspring, mother, father for EA and height
pwcorr  mother_height father_height offspring_height mother_eduyears father_eduyears offspring_educ if offspring_height!=. & offspring_educ!=., star(0.05)   obs

//Genotypic correlations offspring, mother, father for EA and height
pwcorr m_height_score f_height_score o_height_score m_ea2_score f_ea2_score o_ea2_score    if offspring_height!=. & offspring_educ!=.   , star(0.05)      obs

//An example observational association height on EA, unadjusted, adjusted for parental genotype
reg offspring_educ offspring_height,ro
reg offspring_educ offspring_height m_height_score f_height_score,ro 

//An example observational association height, unadjusted and adjusted for parental genotype 
reg offspring_height offspring_educ,ro 
reg offspring_height offspring_educ m_ea2_score f_ea2_score,ro 

//An example MR of height on EA, unadjusted, adjusted for parental genotype
ivreg2 offspring_educ (offspring_height=o_height_score),ro endog(offspring_height)
ivreg2 offspring_educ (offspring_height=o_height_score) m_height_score f_height_score,ro endog(offspring_height)

ivreg2 offspring_educ (offspring_height mother_height father_height =o_height_score  m_height_score f_height_score),ro endog(offspring_height)

