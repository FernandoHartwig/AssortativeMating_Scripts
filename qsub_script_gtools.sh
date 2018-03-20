#!/bin/bash

#This queues up all of the gtools commands to extract all fo the height and BMI snps

for chr in {01..22};
do
  qsub extract_height_snps.sh  -v chr=$chr
  sleep 1
done



