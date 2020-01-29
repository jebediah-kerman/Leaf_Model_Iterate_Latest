#!/bin/bash

# Cleaning garbage
rm -r Output_Summary
mkdir Output_Summary
mkdir Output_Summary/Mesh
mkdir Output_Summary/Density

# ... faire pareil pour toutes les variables
# ! nombres entiers ! il faudra faire la conversion dans le .txt
# D is the total run number calculated as the number of rows in params.txt minus 1
# D is also used as the seed number in main.cpp
D=0
Dmax=$(($(cat params.txt | wc -l)-1))
while read -r i j
do

	# if D is zero (reading the headline), create variable names
	if [ ${D} -eq 0 ]
	then
		iName=${i}
		jName=${j}
		D=$((1+$D))
		continue
	fi

	# Cleaning garbage
	rm -r ${D}_${i}_${j}

	# Creating new directories
    mkdir ${D}_${i}_${j}
    cd ${D}_${i}_${j}
    mkdir Plot
	mkdir Plot/Anis
	mkdir Data
	mkdir Plot/Displ
	mkdir Plot/Elast
	mkdir Plot/Matrix
	mkdir Plot/Mesh
	mkdir Plot/RotMat
	mkdir Plot/Stress
	mkdir Plot/Density
	mkdir Sepal

	# copie des fichiers sauves dans la source
	cp ../Source/*.cpp Sepal/
	echo "                                                 "
	echo "                                                 "
	echo "                                                 "
	echo "                                                 "
	echo "                                                 "
	echo "                                                 "
	echo "                                                 "
	echo "Generation du code"
	echo "Valeur de Repetition    $D"
	echo "                                                 "

	# ecrire dans le fichier
	#echo "include 'iostream'" > Sepal/main.cpp
	#echo "include 'cfloat'" >> Sepal/main.cpp
	#echo "using namespace std;" >> Sepal/main.cpp
	echo "int simnumber=$D;" > Sepal/main.cpp
	echo "int seed=$D;" >> Sepal/main.cpp

	# Parameters to look at
	echo "real ${iName}=${i};" >> Sepal/main.cpp
	echo "real ${jName}=${j};" >> Sepal/main.cpp

	# Paste the end
	cat ../Source/End.cpp >> Sepal/main.cpp

	# executer le fichier
	cd Sepal
	FreeFem++-nw main.cpp

	# Summarize output by picking a mid and final graph in Mesh
	cd ../Plot/Mesh
	mid=$(($(ls | sort -n | wc -l)/2))
	ls | sort -n > tmp1
	tail -n 1 tmp1 > tmp2
	head -n $mid tmp1 | tail -n 1 >> tmp2
	while read name
	do
		cp $name ../../../Output_Summary/Mesh/
	done < tmp2

	# Also pick Density
	cd ../Density
	mid=$(($(ls | sort -n | wc -l)/2))
	ls | sort -n > tmp1
	tail -n 1 tmp1 > tmp2
	head -n $mid tmp1 | tail -n 1 >> tmp2
	while read name
	do
		cp $name ../../../Output_Summary/Density/
	done < tmp2

	# Iteration
    cd ../../..
    D=$((1+$D))

done < params.txt



