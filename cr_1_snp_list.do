//Neil Davies 23/02/18
//This converts the rsids from the SNP list into CHR:POS so that the SNPs can be extracted from the ALSPAC father's data



forvalues i=10(1)22{
	if `i'<10{
		local i="0`i'"
		}
	import delimited "data_chr`i'.snp-stats.",clear
	save "workingdata/fathers_data_cdr`i'_snp-stats",replace
	}

//Select the 697 height variants and the 162 education variants (859 in total):
import delimited "okbay_162.txt", clear 
rename v1 rsid
save "workingdata/okbay_162",replace

import delimited "wood_697.txt", clear 
rename v1 rsid
save "workingdata/wood_697",replace

//Merge with SNP-stats data
forvalues i=1(1)22{
	if `i'<10{
		local i="0`i'"
		}
	use "workingdata/fathers_data_cdr`i'_snp-stats",clear
	rename rsid chr_pos
	rename snpid rsid
	joinby rsid using "workingdata/okbay_162", unmatched(both)
	tab _m
	gen okbay_snp=(_m==3)
	drop _m
	joinby rsid using "workingdata/wood_697", unmatched(both)
	tab _m
	gen wood_snp=(_m==3)
	drop _m
	keep if okbay_snp==1|wood_snp==1
	keep chr_pos rsid
	compress
	save "workingdata/combined_`i'",replace
	}
use "workingdata/combined_22", clear
forvalues i=1(1)21{
	if `i'<10{
		local i="0`i'"
		}		
	append using "workingdata/combined_`i'",
	}
	
joinby rsid using "workingdata/okbay_162", unmatched(both)
tab _m
gen okbay_snp=(_m==3|_m==2)
gen matched=(_m==3)
drop _m
joinby rsid using "workingdata/wood_697", unmatched(both)
tab _m
gen wood_snp=(_m==3|_m==2)
replace matched=1 if _m==3
drop _m

//We have matched 847 SNPs, 12 are missing
keep if chr_pos !=""
compress
gen id=_n
save "snplists/combined_height_educ",replace

//Import Okbay coefficients
import delimited "snplists/EduYears_excl_23andMe_ALSPAC_sumstats.txt", clear 

rename markername rsid
compress
save "snplists/EduYears_excl_23andMe_ALSPAC_sumstats",replace

import delimited "snplists/okbay_162.txt", clear 
rename v1 rsid

joinby rsid using "snplists/EduYears_excl_23andMe_ALSPAC_sumstats",unmatched(master)

drop _m chr pos pval n
gen n=_n
ds effect_allele-se
foreach i in `r(varlist)'{
	rename `i' ea2_`i'
	}
save "snplists/okbaycoefficients",replace

//Import Wood coefficients
//copy from supplementary materials
rename snp rsid
rename var4 height_effect_allele
rename var5 height_other_allele
rename var6 height_eaf
rename var7 height_beta
rename var8 height_se
drop chr position var9 var10
gen n=_n
compress
save "snplists/woodcoefficients",replace
