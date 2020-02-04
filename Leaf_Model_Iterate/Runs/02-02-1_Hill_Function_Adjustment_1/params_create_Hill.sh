# Create heading
echo Rhz$'\t'dRho$'\t'RelEl > params.txt

# Create parameter list
for i in 0.01 0.03 0.1
do
	for j in 2 4 20
	do
		for k in 5. 10. 50.
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
echo 0.1$'\t'20$'\t'1. >> params.txt