#!/usr/bin/python
from sys import argv
import csv

file_list=argv[1:len(argv)-1]
a=len(file_list)
readers = []
for fileI in file_list:				# creates a list of lists of rows (rows are strings!)
	f=open(fileI,'r').readlines()
	readers.append(f)
dictionaries = [dict() for i in range(a)] #creates a list of 'a' empty dictionaries


for i,dictI in enumerate(dictionaries):
	for j, rowJ in enumerate(readers[i]):
			new_row=rowJ.rstrip('\n').split(',')
			if j != 0 and len(new_row) == 2:
				dictI[new_row[0]]=new_row[1]

dict_keys = [dictI.keys() for dictI in dictionaries]

for i,keysI in enumerate(dict_keys[1:]):	#i goes from 0 and dict_keys goes from 1 so we need to always add 1 to i
	for n in keysI:
		dictI = dictionaries[i+1]	
		dict_prev = dictionaries[i]
		z = dictI[n]
		try:
			if z in dict_keys[i]:
				dictI[n] = dict_prev[z]
			else:
				del dictI[n]
		except KeyError:
			del dictI[n]


'''
for i,file in enumerate(file_list):
	print 'this is file number ',i,' : ',file
'''

output=open(argv[-1],'w')        #argv[6] should be a file that contains labels of T5 with parent labels from T0
writer=csv.writer(output,delimiter=',')
headline=('Label','Parent label')
writer.writerow(headline)
for item in dictionaries[-1]:
	row=(item,dictionaries[-1][item])
	writer.writerow(row)
#this loop for each item in dictionary creates a list containing key and value and then writes this list as one row in csv file
output.close()


