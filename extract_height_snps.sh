#!/bin/bash


output_dir="height_father_snps_${chr}.gen"
snp_list="alspac_recode.txt"

echo "HERE"

genfile="data_chr${chr}.bgen"
samplefile="data.sample"

# qctool -g example.bgen -og subsetted.gen -incl-rsids rsids_to_include.txt
# The bgen files do not have RSids only chromsome and position.
# Will extract the CHR:POS from the SNP stats file and repeat.




# for chr in `seq 17 22`;
# do
#   echo $chr
echo "Running gtools extraction"
echo "${gtool}  -g ${genfile} -incl-rsids ${snp_list} -og ${output_dir}"
module add apps/qctool-2.0 

qctool -g ${genfile} -og ${output_dir} -incl-rsids ${snp_list} 
# done




