# Create heading
echo RhzInh$'\t'dRhoInh$'\t'RelElInh > params.txt

# Create parameter list
for i in 0.05 0.1 0.2
do
	for j in 3. 6. 20.
	do
		for k in 5. 10. 20.
		do
			echo ${i}$'\t'${j}$'\t'${k} >> params.txt
		done
	done
done
