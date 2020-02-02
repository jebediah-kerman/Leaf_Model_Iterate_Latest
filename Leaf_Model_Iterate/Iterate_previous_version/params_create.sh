# Create heading
echo frontArrHeightIni$'\t'fAspeed$'\t'fAElastFactor > params.txt

# Create parameter list
for i in 2 4 6
do
	for j in 0.025 0.05 0.1
	do
		for k in 1. 2. 10.
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
