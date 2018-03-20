//Neil Davies 30/01/18
//This merges the genotype and phenotype data

use "workingdata/phenotypes",clear

tostring aln,replace
joinby aln using "workingdata/snpdata_scores",unmatched(using) _merge(XX)
compress
save "workingdata/formerging",replace
