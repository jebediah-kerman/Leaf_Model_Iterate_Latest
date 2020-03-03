# Create heading
echo DInh$'\t'GI2transitionalGFconc$'\t'junk3 > params.txt

# Create parameter list
for i in 0.1 0.2 0.4
do
	for j in 0.0003 0.0001 0.00003 0.00001
	do
		for k in 1.
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
