# Create heading
echo Rhz$'\t'dRho$'\t'RelEl > params.txt

# Create parameter list
for i in 0.02 0.03 0.05
do
	for j in 1.5 2 3
	do
		for k in 7 10 15
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done