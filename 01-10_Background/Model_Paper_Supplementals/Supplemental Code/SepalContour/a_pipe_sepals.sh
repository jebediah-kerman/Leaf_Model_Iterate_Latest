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
	logfile='pipe_sepals.log'
	echo '--------------------------' >> $logfile
	date >> $logfile
	echo '--------------------------' >> $logfile

	# write the contours to R, measure length, width, etc...
	# -----------------------
	rm ContoursData_sample$Nsample/*
	
	mkdir registered_curves
count=0
unset Rejlist
for file in *.jpeg #*.tif #*.jpeg *.tiff
do
echo $file;
../../SepalContour/detect_2Dcontour-scaled.py $file $Nsample v
	cd ..
done
done
