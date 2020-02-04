# Create heading
echo AnisStart$'\t'AnisEnd$'\t'junk > params.txt

# Create parameter list
for i in 0.8 0.6 0.4
do
	for j in 0.2 0.1 0.0
	do
		for k in 1
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
