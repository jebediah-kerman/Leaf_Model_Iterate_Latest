# Create heading
echo Rhz$'\t'dRho > params.txt

# Create parameter list
for j in 0.01 0.03 0.1 0.3
do
	for k in 2. 4. 6. 99.
	do
		echo ${j}$'\t'${k} >> params.txt
	done
done
