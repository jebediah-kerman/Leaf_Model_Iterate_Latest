#!/bin/bash

# Cleaning garbage
rm -r Output_Summary
mkdir Output_Summary
mkdir Output_Summary/Mesh

# ... faire pareil pour toutes les variables
# ! nombres entiers ! il faudra faire la conversion dans le .txt
# D is also used as the seed number in main.cpp
D=1
while [ $D -le 2 ]
do
	# Cleaning garbage
	rm -r $D

	# Creating new directories
    mkdir $D
    cd $D
    mkdir Plot
	mkdir Plot/Anis
	mkdir Data
	mkdir Plot/Displ
	mkdir Plot/Elast
	mkdir Plot/Matrix
	mkdir Plot/Mesh
	mkdir Plot/RotMat
	mkdir Plot/Stress
	mkdir Plot/Analyses
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
	cat ../Source/End.cpp >> Sepal/main.cpp

	# executer le fichier
	cd Sepal
	FreeFem++-nw main.cpp

	# Summarize output by picking a mid and final graph
	cd ../Plot/Mesh
	mid=$(($(ls | sort -n | wc -l)/2))
	ls | sort -n > tmp1
	tail -n 1 tmp1 > tmp2
	head -n $mid tmp1 | tail -n 1 >> tmp2
	while read name
	do
		cp $name ../../../Output_Summary/Mesh/
	done < tmp2

	# Iteration
    cd ../../..
    D=$((1+$D))

done



