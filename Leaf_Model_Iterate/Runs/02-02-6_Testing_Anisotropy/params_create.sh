# Create heading
echo Anis$'\t'junk1$'\t'junk2 > params.txt

# Create parameter list
for i in 0. 0.1 0.2 0.3 0.4
do
	for j in 1
	do
		for k in 1
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
