//Neil Davies 26/01/18
//This cleans the ALSPAC fathers dosage data


cap prog drop import_snps
prog def import_snps
args file


import delimited "`file'.", delimiter(space, collapse) encoding(ISO-8859-1) clear

local chr=v[2]
local pos=v[4]
local effect=v[5]
local other=v[6]

drop if _n<7 
gen id=round(_n/3+1/6)
bys id: gen n=_n
drop if n==3

destring v, replace
replace v=v*2 if n==1 

bys id: egen X`chr'_`pos'_`effect'_`other'=total(v)
drop if n==2
drop v

drop n

save "`file'_clean.dta",replace

end



fs split_height_father*_trans

fs 
foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}

//Merge all of the files back togeather
use split_height_father_snps_01.genaaaa_trans_clean.dta,clear
fs *clean.dta
foreach i in `r(files)'{
	local j=`j'+1
	di "`j'"
	di "`i'"
	joinby id using `i' 
	}
save "snps",replace

//Merge in the rsids
use "rawdata/extracted_snps/split/snps.dta",clear
joinby id using "snplists/combined_height_educ",unmatched(master)

gen chr=word(subinstr(chr_pos,":"," ",1),1)
gen pos=word(subinstr(chr_pos,":"," ",1),2)


gen effect=""
gen other=""

local rsid="X"
while rsid!=""{
	local i=`i'+1
	local rsid=rsid[`i']
	local chr=chr[`i']
	local pos=pos[`i']

	if length("`chr'")==1{
		local chr="0`chr'"
		}
		
	ds X`chr'_`pos'
	local effect=word(subinstr("`r(varlist)'","_"," ",5),3)
	local other=word(subinstr("`r(varlist)'","_"," ",5),4)

	
	di "CHR=`chr'"
	di "POS=`pos'"
	di "rsid=`rsid'"
	di "effect=`effect'"
	di "other=`other'"
	
	replace effect ="`effect'" in `i'
	replace other ="`other'" in `i'
	
	rename X`chr'_`pos' `rsid'_`effect'_`other'
	}
keep id  rs*
compress
save "workingdata/fathers_snpdata",replace

//847 SNPs extracted
//Merge in the IDs from the sample file

import delimited "rawdata/data.sample", delimiter(space) clear 
drop in 1
keep id_1 
gen n=_n
save "workingdata/ids",replace

use "workingdata/fathers_snpdata",clear
gen n=_n
joinby n using  "workingdata/ids",unmatched(master)
drop _m
order id_1
rename id_1 aln
drop id n
compress
save "workingdata/fathers_snpdata",replace

//*****************************
//MOTHERs 
//*****************************


cap prog drop import_snps
prog def import_snps
args file

import delimited "`file'.", delimiter(space, collapse) encoding(ISO-8859-1) clear

local rsid=v[2]
local effect=v[4]
local other=v[5]

drop if _n<6
gen id=round(_n/3+1/6)
bys id: gen n=_n
drop if n==3

destring v, replace
replace v=v*2 if n==1 

bys id: egen `rsid'_`effect'_`other'=total(v)
drop if n==2
drop v

drop n

save "`file'_clean.dta",replace

end

cd "extracted_snps/split"

fs split_height_mother*_trans

foreach i in `r(files)' {
	di "`i'"
	import_snps `i'
	}

//Merge all of the files back togeather
use split_height_mother_snps_01.genaaaa_trans_clean.dta,clear
fs split_height_mo*clean.dta
foreach i in `r(files)'{
	local j=`j'+1
	di "`j'"
	di "`i'"
	joinby id using `i' 
	}

compress
save "workingdata/snps_mother",replace

//Merge in the IDs from the sample file

import delimited "rawdata/data.sample_mothers.", delimiter(space) clear 
drop in 1
keep id_1 
gen n=_n
save "workingdata/mothers_ids",replace

//849 SNPs extracted
//Merge in the IDs from the sample file

use "workingdata/snps_mother",clear
gen n=_n
joinby n using  "workingdata/mothers_ids",unmatched(master)
drop _m
order id_1
rename id_1 aln
drop id n
compress
save "workingdata/mothers_snpdata",replace

append using "workingdata/fathers_snpdata"

save "workingdata/combined_snpdata",replace

//Constructing allele scores
use "workingdata/combined_snpdata",clear
gen n=_n

//Next we need to extract the coefficients from Okbay and Wood
//And construct the scores

joinby n using "snplists/woodcoefficients",unmatched(master) update
drop _m
joinby n using "snplists/okbaycoefficients",unmatched(master) update
drop _m

destring height_eaf,replace

//Construct the scores:
//Harmonize the effect alleles
gen height_als_effect=""
gen height_als_other=""
gen height_als_eaf=.

forvalues i=6(1)702{
	local rsid=height_rsid[`i']
	cap:ds `rsid'
	if _rc==111{
		di "SNP `rsid' not found"
		}
	else{
		local snpname="`r(varlist)'"
		replace height_als_effect=word(subinstr("`snpname'","_"," ",5),2) in `i'
		replace height_als_other=word(subinstr("`snpname'","_"," ",5),3) in `i'
		tabstat `rsid',save
		replace height_als_eaf=el(r(StatTotal),1,1)/2 in `i'
		}
	}

//Harmonize height alleles
gen flip=1 if height_als_effect==height_other_allele & height_als_other==height_effect_allele

replace height_beta=height_beta*-1 if flip==1
replace height_eaf=1-height_eaf if flip==1
replace height_effect_allele=height_als_effect if flip==1
replace height_other_allele=height_als_other if flip==1

twoway scatter height_eaf height_als_eaf 

//Harmonize education alleles
gen ea2_als_effect=""
gen ea2_als_other=""
gen ea2_als_eaf=.

forvalues i=1(1)162{
	local rsid=ea2_rsid[`i']
	cap:ds `rsid'
	if _rc==111{
		di "SNP `rsid' not found"
		}
	else{
		local snpname="`r(varlist)'"
		replace ea2_als_effect=word(subinstr("`snpname'","_"," ",5),2) in `i'
		replace ea2_als_other=word(subinstr("`snpname'","_"," ",5),3) in `i'
		tabstat `rsid',save
		replace ea2_als_eaf=el(r(StatTotal),1,1)/2 in `i'
		}
	}

replace flip=.

replace flip=1  if ea2_als_effect==upper(ea2_other_allele) & ea2_als_other==upper(ea2_effect_allele)

twoway scatter ea2_eaf ea2_als_eaf 

replace ea2_beta=ea2_beta*-1 if flip==1
replace ea2_eaf=1-ea2_eaf if flip==1
replace ea2_effect_allele=ea2_als_effect if flip==1
replace ea2_other_allele=ea2_als_other if flip==1

twoway scatter ea2_eaf ea2_als_eaf 

//Create scores
gen ea2_score=0
gen height_score=0
forvalues i=6(1)702{
	local rsid=height_rsid[`i']
	cap:ds `rsid'
	if _rc==111{
		di "SNP `rsid' not found"
		}
	else{
		di "`rsid' `beta' `i'"
		local beta=height_beta[`i']
		replace height_score=`rsid'*`beta'+height_score if `rsid'!=.
		}
	}

forvalues i=1(1)162{
	local rsid=ea2_rsid[`i']
	cap:ds `rsid'
	if _rc==111{
		di "SNP `rsid' not found"
		}
	else{
		di "`rsid' `beta' `i'"
		local beta=ea2_beta[`i']
		replace ea2_score=`rsid'*`beta'+ea2_score if `rsid'!=.
		}
	}
compress	
save "workingdata/snpdata_scores",replace
use "workingdata/snpdata_scores",clear

preserve 
keep if qlet=="F"
keep aln ea2_score height_score rs*
replace aln=substr(aln,1,5)
save "workingdata/snpdata_scores_f",replace
restore
  	
preserve 
keep if qlet=="M"
keep aln ea2_score height_score rs*
replace aln=substr(aln,1,5)
save "workingdata/snpdata_scores_m",replace
restore

preserve 
keep if qlet=="A"
keep aln ea2_score height_score rs*
replace aln=substr(aln,1,5)
save "workingdata/snpdata_scores_o",replace
restore


use "workingdata/snpdata_scores_o",clear
renpfix rs O_rs
rename height_score O_height_score
rename ea2_score O_ea2_score
joinby aln using "workingdata/snpdata_scores_f",unmatched(none)

renpfix rs F_rs
rename height_score F_height_score
rename ea2_score F_ea2_score

joinby aln using "workingdata/snpdata_scores_m",unmatched(none)

renpfix rs M_rs
rename height_score M_height_score
rename ea2_score M_ea2_score

//Check for mismatches
ds O_*
foreach i in `r(varlist)'{
	local rsid=word(subinstr("`i'","_"," ",5),2)
	di "`rsid'"
	gen diff_`rsid'=F_`rsid'/2+ M_`rsid'/2-O_`rsid'
	}

egen mean_diff=rowmean(diff_*)
sum mean_diff
compress
save "workingdata/snpdata_scores",replace

