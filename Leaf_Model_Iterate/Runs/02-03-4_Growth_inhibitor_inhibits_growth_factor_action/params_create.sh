# Create heading
echo RhzInh$'\t'dRhoInh$'\t'RelElInh > params.txt

# Create parameter list
for i in 0.05 0.1
do
	for j in 3. 20.
	do
		for k in 1. 5. 10.
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
