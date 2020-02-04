# Create heading
echo frontArrHeightIni$'\t'fAspeed$'\t'fADiffusionFactor > params.txt

# Create parameter list
for i in 2 4 6
do
	for j in 0.5 1.0
	do
		for k in 0.0 0.5 1.0
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done

for i in 8 10
do
	for j in 0.05 0.1 0.2
	do
		for k in 0.0 0.5 1.0
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done

for i in 8 10
do
	for j in 0.5 1.0
	do
		for k in 0.0 0.5 1.0
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done