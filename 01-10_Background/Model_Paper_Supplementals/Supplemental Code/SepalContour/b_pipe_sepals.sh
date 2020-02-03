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
count=0
unset Rejlist
for file in *.tif #*.tif #*.jpeg *.tiff
do
echo $file;	count=`expr $count + 1`
filetemp=$(basename $(echo $file) .tif)
all="$type-$genotype-$filetemp"
echo $all
# check if the file is to be rejected
if grep -q $all ../Rejected.txt;then echo $count;Rejlist+="-$count";fi
done
		# registering the data...
		# ------------------------
	
echo registering data
echo =========
echo $Nsample
echo $Rejlist

mkdir registration_on_common_landmarks$Nsample
../../SepalContour/fda_registerSepals.R $Nsample $Rejlist
	cd ..
done
