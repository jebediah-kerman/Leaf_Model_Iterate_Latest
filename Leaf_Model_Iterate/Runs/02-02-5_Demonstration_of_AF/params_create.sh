# Create heading
echo frontArrHeightIni$'\t'junk1$'\t'junk2 > params.txt

# Create parameter list
for i in 4. 5. 6. 7. 8.
do
	for j in 1
	do
		for k in 1
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
