# Create heading
echo frontArrHeightIni$'\t'fAspeed$'\t'fAElastFactor > params.txt

# Create parameter list
for i in 5 6 7
do
	for j in 0.05 0.1 0.2
	do
		for k in 1.5 2. 3.
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
