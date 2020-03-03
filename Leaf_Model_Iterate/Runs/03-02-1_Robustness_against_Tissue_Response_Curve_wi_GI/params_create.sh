# Create heading
echo Rhz$'\t'junk2$'\t'junk3 > params.txt

# Create parameter list
for i in 0.01 0.02 0.03 0.04 0.05
do
	for j in 1.
	do
		for k in 1.
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
