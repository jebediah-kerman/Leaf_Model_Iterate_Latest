#!/bin/bash

# Cleaning garbage
#rm -r Output_Summary
#mkdir Output_Summary
#mkdir Output_Summary/Elast
#mkdir Output_Summary/Density
#mkdir Output_Summary/GrowthRate

# ... faire pareil pour toutes les variables
# ! nombres entiers ! il faudra faire la conversion dans le .txt
# D is the total run number calculated as the number of rows in params.txt minus 1
# D is also used as the seed number in main.cpp
D=100
Dmax=$(($(cat params.txt | wc -l)-1))

# Header of good_output.txt
#echo SimNumber$'\t'i$'\t'j$'\t'k$'\t'Crit_AFInit$'\t'Crit_HWRatio > Output_Summary/good_output.txt



while read -r i j k
do

	# if D is zero (reading the headline), create variable names
	if [ ${D} -eq 100 ]
	then
	iName=${i}
		jName=${j}
		kName=${k}
		D=$((1+$D))
		continue
	fi

	# Three replications
	for rep in 1 2 3
	do

		# Cleaning garbage
		rm -r ${i}_${j}_${k}_${D}

		# Creating new directories
		mkdir ${i}_${j}_${k}_${D}
		cd ${i}_${j}_${k}_${D}
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
		mkdir Plot/DiffusionConst
		mkdir Plot/GrowthRate
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
		echo "real ${kName}=${k};" >> Sepal/main.cpp

		# Paste the end
		cat ../Source/End.cpp >> Sepal/main.cpp

		# Paste additional outputs
		echo fferr\ \<\<\ \"${D}\\t${i}\\t${j}\\t${k}\\t\"\ \<\<\ CritAFInit\ \<\<\ \"\\t\"\ \<\<\ CritHWRatio\ \<\<\ endl\; >> Sepal/main.cpp

		# executer le fichier
		cd Sepal
		FreeFem++-nw main.cpp

		# Summarize output by picking a mid and final graph in Elasticity
		cd ../Plot/Elast
		mid=$(($(ls | sort -n | wc -l)/2))
		ls | sort -n > tmp1
		tail -n 1 tmp1 > tmp2
		## head -n $mid tmp1 | tail -n 1 >> tmp2	# Output mid-point?
		while read name
		do
			cp $name ../../../Output_Summary/Elast/${i}_${j}_${k}_${D}_${name}
		done < tmp2

		# Also Density
		cd ../Density
		mid=$(($(ls | sort -n | wc -l)/2))
		ls | sort -n > tmp1
		tail -n 1 tmp1 > tmp2
		## head -n $mid tmp1 | tail -n 1 >> tmp2	# Output mid-point?
		while read name
		do
				cp $name ../../../Output_Summary/Density/${i}_${j}_${k}_${D}_${name}
		done < tmp2

		# Also GrowthRate
		cd ../GrowthRate
		mid=$(($(ls | sort -n | wc -l)/2))
		ls | sort -n > tmp1
		tail -n 1 tmp1 > tmp2
		## head -n $mid tmp1 | tail -n 1 >> tmp2	# Output mid-point?
		while read name
		do
			cp $name ../../../Output_Summary/GrowthRate/${i}_${j}_${k}_${D}_${name}
		done < tmp2

		# Iteration
	    cd ../../..
	    D=$((1+$D))

	done
done < params.txt



