# Create heading
echo fAElastFactor$'\t'dRho$'\t'junk > params.txt

# Create parameter list
for i in 2. 3. 5.
do
	for j in 2 3 4
	do
		for k in 1
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
