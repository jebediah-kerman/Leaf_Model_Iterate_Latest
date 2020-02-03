#! /bin/bash
# pipe_sepals.sh Nsample
# 

Nsample=400 # number of points on the half leaf

scale_type=${PWD##*/}

# Scale=ALL_DATA # Mingyuan # Lilan # says which scale to use

# go in each genotype folder, and then in each type folder

unset genlist
unset typelist
for genotype in */
do
	echo $genotype
	cd $genotype
	genlist+="$genotype"
	genotype=${genotype%?}
	cd ..
done




echo PCA analysis
echo =========
# PCA analysis ...
# ------------------------

mkdir ../ImagesPCA_$scale_type
../SepalContour/fda_PCA_Sepals.R $Nsample $genlist $scale_type
