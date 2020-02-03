#!/usr/bin/python
# -*- coding: utf-8 -*-

##################################################
##################################################
##                                              ##
##                 Freefem++ script             ##
##                                              ##
##  author: Mathilde Dumond                     ##
##                                              ##
##  Supplemental Information of the paper:      ##
##  Variable cell growth yields reproducible    ##
##  organ development through spatiotemporal    ##
##  averaging                                   ##
##                                              ##
##  Hong et al.                                 ##
##  Dev Cell 2016                               ##
##                                              ##
##################################################
##################################################


from __future__ import print_function
import sys
import os
import matplotlib.pyplot as plt
import Image


from class_mesh import *


def inter(a,b):
	return b,a

def rearrange(liste, ind):
	liste2 = liste[ind:]
	liste2.extend(liste[:ind])
	return liste2

def ReadFile(filename):
	NewTime = mesh()
	NewTime.readFile(filename)
	return NewTime

def tofloatSt(string, line):
	filename = "Output_names.txt"
	st = string.replace(' ','')
	if st == '':
		return None
	else:
		try:
			st = float(st)
			return st
		except ValueError:
			print("Oops, that was not a valid number: "+ st + "\nAt line: "+str(line)+" in the file "+filename)
			raise

def ReadAllNames():
	'''
	Reads the file Output_names.txt
	'''
	filename = "Output_names.txt"
	f = open(filename,'r')
	csvin = csv.reader(f, delimiter=';')
	listAllOutput = []
	listAllTitles = []
	listAllTypes = []
	listAllColors = []
	listAllScales = []
	f.next() # header line
	i=1
	for row in csvin:
		listAllOutput.append(row[0].replace(' ',''))
		listAllTitles.append(row[1].replace(' ',''))
		listAllTypes.append(row[2].replace(' ',''))
		listAllColors.append(row[3].replace(' ',''))
		listAllScales.append([tofloatSt(row[4], i), tofloatSt(row[5], i)])
		i+=1
	f.close()
	return [listAllOutput, listAllTitles, listAllTypes ,listAllColors, listAllScales]

def cleanListStrings(l):
	'''
	list of strings -> list of strings
	removes all blank spaces from each string in the list of strings
	'''
	k = []
	for i in l:
		k.extend(i.replace('\n','').split(' '))
	k = filter(lambda a: a != '', k) # remove the empty strings
	return k

def ReadToAnalyze():
	'''
	Reads the file Input_names and gets the info
	'''
	filename = "Input_names.txt"
	filein = open(filename,'r')
	ToDoGraphs = []
	ToDoOvTimeGraphs = []
	line = filein.next()
	while line[0] == '#':
		line = filein.next()
	while not '# End ToDoGraphs' in line: # read until this line = all the info for the sepal graphs
		# the empty lines will go into the list, but will be removed in the clean function
		line = line.split(',')
		ToDoGraphs.extend(line)
		line = filein.next()
	while line[0] == '#':
		line = filein.next()
	while not '# End ToDoOvTimeGraphs' in line: # read the lines for the over time graphs
		# the empty lines will go into the list, but will be removed in the clean function
		line = line.split(',')
		ToDoOvTimeGraphs.extend(line)
		line = filein.next()
	filein.close()
	ToDoGraphs = cleanListStrings(ToDoGraphs)
	ToDoOvTimeGraphs = cleanListStrings(ToDoOvTimeGraphs)
	return ToDoGraphs, ToDoOvTimeGraphs

def OneTimePoint(NewTime, currFile, sevSim=False):
	'''
	Main function for the analysis of one data point
	'''
	listAll = ReadAllNames()
	listOutput, _ = ReadToAnalyze()
	if listOutput == []:
		return
	NewTime.setSubList() # if you want to display a subset of all the points - the list is suuposed to give equally distributed points - could be improved
	print('Analyzing the '+NewTime.getPrefix()+' file.')
	numb = 1
	for i in listOutput:
		# print(i)
		sys.stdout.flush()
		# find if 'D', find if 'LOCVAR'
		name, ind, title, typ, col, limScale = findSpec(i, listAll)
		# different treatment for each kind (typ) of data
		if typ == 'M': # Mesh
			TrIndex, Index, [X, Y] = eval('NewTime.get'+i+'()')
			NewTime.plotMesh(TrIndex, Index, [X, Y],title)
		elif typ == 'S': # Scalar
			Index, Value = NewTime.getScalar(name)
			NewTime.plotScal(Index, Value, title, sevSim, colormap = col, colorrange = limScale)
		elif typ == 'D': #'distribution'
			cF = open(currFile, 'a')
			cF.write(NewTime.getPrefix()+'\n\n')
			cF.close()
			Index, Value = NewTime.getScalar(name)
			NewTime.plotDistr(Value, title, currFile)
		elif typ == 'T': # scalar defined on triangles (and not vertex)
			TrIndex, Index, [X, Y], Z = eval('NewTime.get'+name+'()')
			NewTime.plotScalTr(TrIndex, Index, [X, Y], Z, title, sevSim=sevSim, col=col, limScale=limScale)
		elif typ == 'O':
			pass
		# display
		if numb/len(listOutput) == 1:
			print(str(numb*100/len(listOutput))+'%'+' done!')
		else:
			print(str(numb*100/len(listOutput))+'%'+' done!', end='\r')
			sys.stdout.flush()
		numb+=1

def findSpec(name, listAll):
	'''
	Returns whether there is a 'D' and/or a 'LOCVAR' substring in the string
	and the canonical name
	If you want to add other specifications do it here
	'''
	distr = False
	if name[-1] == 'D':
		distr = True
		name = name[:-1]
	ind = listAll[0].index(name)
	# find the link in listAll (= Output_names.txt)
	if distr:
		title = listAll[1][ind]+'_dist'
		typ = 'D'
	else:
		title = listAll[1][ind]
		typ = listAll[2][ind]
	col = listAll[3][ind]
	limScale = listAll[4][ind]
	if name == 'TrContour':
		name = 'One'
	return name, ind, title, typ, col, limScale

def OnePoint(name, currFile):
	NewTime = ReadFile(name)
	OneTimePoint(NewTime, currFile, sevSim=False)

def AllTimePoints(simname, currFile):
	'''
	Goal of this function = plot all the point to easily make a video using ImageJ
	'''
	listOutput, _ = ReadToAnalyze()
	simname = simname+'/'
	nbTotIt = len([name for name in os.listdir(simname+'Data/')])
	iterations = [i for i in range(nbTotIt)]
	iterations = [str(i) for i in iterations]
	iterations = [(5-len(i))*'0'+i for i in iterations]
	iterFiles = []
	# First, make the list of all files to be treated
	for i in iterations:
		for j in os.listdir(simname+'Data/'):
			if i in j:
				iterFiles.append(j)
				break
	iterFiles = [simname+'Data/'+i for i in iterFiles]
	# Then, treat the files
	numb = 1
	for i in iterFiles:
		# print(i)
		NewTime = ReadFile(i)
		# cF = open(currFile, 'a')
		# cF.write(NewTime.getPrefix()+'\n\n')
		# cF.close()
		OneTimePoint(NewTime, currFile, listOutput)
		# display
		if numb/len(listOutput) == 1:
			print(str(numb*100/len(iterFiles))+'%'+' done!')
		else:
			print(str(numb*100/len(iterFiles))+'%'+' done!', end='\r')
			sys.stdout.flush()

def SevTimePoints(simname, currFile, nbPts=[3,3]):
	'''
	Goal of this function: get a panel of different time (equally distributed) points showing the evolution of the simulation
	If you want more than 9 steps, please change the nbPts default argument
	'''
	listOutput, _ = ReadToAnalyze()
	if listOutput == []:
		return
	simname = simname+'/'
	# total number of iterations
	nbTotIt = len([name for name in os.listdir(simname+'Data/')])
	# divided by the number of steps
	nbIt = nbPts[0]*nbPts[1]
	step = nbTotIt/(nbIt-1)
	if step*(nbIt-1)<(nbTotIt-2):
		step +=1
	iterations = [min(i*step, (nbTotIt-2)) for i in range(nbIt)]
	iterations = [str(i) for i in iterations]
	iterations = [(5-len(i))*'0'+i for i in iterations]
	iterFiles = []
	for i in iterations:
		for j in os.listdir(simname+'Data/'):
			if i in j:
				iterFiles.append(j)
				break
	iterFiles = [simname+'Data/'+i for i in iterFiles]
	# analysis of each point
	for i in iterFiles:
		# print(i)
		NewTime = ReadFile(i)
		# cF = open(currFile, 'a')
		# cF.write(NewTime.getPrefix()+'\n\n')
		# cF.close()
		OneTimePoint(NewTime, currFile, listOutput)
	# Merge of the graphs done independently
	MergeGraphs(iterations, simname, nbPts)

def LastTimePoint(simname, currFile, sevSim = False):
	'''
	Analysis of only the last time point
	'''
	# Find which file
	simname = simname+'/'
	nbTotIt = len([name for name in os.listdir(simname+'Data/')])
	j = os.listdir(simname+'Data/')
	for tfile in j:
		if ('00'+str(nbTotIt-1)) in tfile:
			lastFile = tfile
	lastFile = simname+'Data/'+lastFile
	# cF = open(currFile, 'a')
	# cF.write(lastFile.split('/')[-1]+'\n\n')
	# cF.close()
	NewTime = ReadFile(lastFile)
	OneTimePoint(NewTime, currFile, sevSim=sevSim)
	if sevSim:
		NewTime.saveContour()
	return NewTime.getLength(), NewTime.getWidth(), NewTime.getArea(), NewTime.getPerim()

def MergeGraphs(iterations, simname, nbPts):
	'''
	Finds all png files corresponding to the same data, merges them all together and deletes the original files
	Deals with the iteration numbers so that if you did another step it is not erased
	'''
	listOutput, _ = ReadToAnalyze()
	# First list of what is there
	f = []
	for (dirpath, dirnames, filenames) in os.walk(simname+'Plot/Analyses/'):
		f.extend(filenames)
		break
	f.sort()
	listAll = ReadAllNames()
	# sort the files
	allFiles = SortByType(f, iterations, simname)
	# merge
	for l in allFiles.keys():
		MergeAGraph(allFiles[l], nbPts, l, simname)

def SortByType(listFiles, iterations, simname):
	'''
	from a list of graphes, gives back a dictionnary with as key the title of the graph, and as entry the ordered list of all png files
	@argument: list of files where the step number and the name are separated by '___'
	'''
	listFilesO = listFiles[:]
	dic = {}
	for f in listFilesO:
		f2=f.split('.')[0] # removes the .png
		f2=f2.split('___') # separator between the iteration number (f2[0]) and the graph type (f2[1])
		it = f2[0].split('_')[-1]
		if it in iterations:
			if f2[1] in dic.keys():
				dic[f2[1]].append(simname+'Plot/Analyses/'+f)
			else:
				dic[f2[1]] = [simname+'Plot/Analyses/'+f]
	return dic

def MergeAGraph(listFiles, nbPts, title, simname):
	'''
	From a list of png files, makes a big png file with all of them
	'''
	# print(listFiles)
	firstIm = Image.open(listFiles[0])
	# dimension of the file (same for every file)
	dim = firstIm.size
	# big image initialization
	blankIm = Image.new("RGB",(dim[0]*nbPts[0],dim[1]*nbPts[1]))
	fileIndex = 0
	# print(listFiles)
	# open each file
	for j in range(nbPts[1]):
		for i in range(nbPts[0]):
			tempIm = Image.open(listFiles[fileIndex])
			# find where to copy the new image
			coord1 = dim[0]*i
			coord2 = dim[1]*j
			# stitch
			blankIm.paste(tempIm, (coord1,coord2))
			# delete the original file
			os.remove(listFiles[fileIndex])
			# increase index
			fileIndex += 1
	blankIm.save(simname +'Plot/Analyses/'+title+'.png')

def InitializeDataFile(NewTime):
	'''
	Initalizes the file in which are stored the over time data (xmin & co) - each simulation has one file
	'''
	name = str(NewTime.getOutputFolder())
	name = name.split('/')
	outputFold = str(NewTime.getOutputFolder())+'../../'+name[-4]
	currOTFile = open(outputFold+'_overTimeData.txt','w')
	currOTFile.write('t; xmin; xmax; ymin; ymax\n')
	currOTFile.close()
	return outputFold+'_overTimeData.txt'

def GraphsOverTime(simname, currFile):
	'''
	Main function for the graphs 'over time' - to have the evolution of a quantity over time
	'''
	print('Doing over time graphs.')
	_, listOutputG = ReadToAnalyze()
	if listOutputG == []:
		return
	listAll = ReadAllNames()
	# list of all time steps
	f = []
	for (dirpath, dirnames, filenames) in os.walk(simname+'/Data/'):
		f.extend(filenames)
		break
	f.sort()
	f = f[:-1] # the last file is removed because if the simulation crashed this file cannot be used
	AllData = [[] for n in listOutputG] # one sublist for each kind of data
	Step = []
	ini = True
	numb = 1
	maxi = len(f)
	for i in f:
		if numb/maxi == 1:
			print(str(numb*100/maxi)+'%'+' done!')
		else:
			print(str(numb*100/maxi)+'%'+' done!', end='\r')
			sys.stdout.flush()
		numb+=1
		NewTime = mesh()
		NewTime.readFile(simname+'/Data/'+i)
		for j in range(len(listOutputG)):
			if listOutputG[j]=='Data':
				if ini: # initialization = header of the file
					ini = False
					currOTFileName = InitializeDataFile(NewTime)
				DataOverTime(currOTFileName,i , NewTime)
			else:
				ind = listAll[0].index(listOutputG[j])
				title = listAll[1][ind]
				typ = listAll[2][ind]
			if typ == 'O':
				Value = eval('NewTime.get'+listOutputG[j]+'()')
			else:
				Value = NewTime.getDistrParam(listOutputG[j])
			# Value is either just a float (length), or a list [mean, sd] and added in the AllData structure
			AllData[j].append(Value)
			s = getStepfromName(i)
			if j == 0:
				Step.append(s)
	# Doing the plots
	for i in range(len(listOutputG)):
		if listOutputG[i]=='Data':
			pass
		else:
			ind = listAll[0].index(listOutputG[i])
			title = listAll[1][ind]
			color = listAll[3][ind]
			PlotTime(simname, Step, AllData[i], title, col = color)

def DataOverTime(currOTFileName,t, NewTime):
	xmin, xmax, ymin, ymax = NewTime.getxyminmax()
	currOTFile = open(currOTFileName, 'a')
	currOTFile.write(t+'; '+str(round(xmin,2))+'; '+str(round(xmax,2))+'; '+str(round(ymin,2))+'; '+str(round(ymax,2))+'\n')
	currOTFile.close()

def SevSim(foldername, currFile, overTime=False):
	'''
	Main function dealing with several simulations at once
	For each simulation, launches the main function of one simulation
	'''
	listSim = next(os.walk(foldername))[1]
	listSim.remove('Source') # 'Source' is the name of the template of the simulation code
	# The GoodOuput file allows to select the simulations to be analyzed (for example on whether the simulation crashed)
	GoodOutput = []
	numb = 1
	# if there is no such file, it won't crash, and analyze all simulations in the folder
	try:
		filename = foldername + '/good_output.txt'
		csvin = csv.reader(open(filename,'r'), delimiter=' ')
		for row in csvin:
			GoodOutput.append(str(row[0]))
		newlistSim = []
		for i in GoodOutput:
			newlistSim.append(i)
		listSim = newlistSim
	finally:
		listSim = [foldername+'/'+i for i in listSim]
		Id = []
		Lengths = []
		Widths = []
		Areas = []
		Perimeters = []
		for sim in listSim:
			print ('Analyzing simulation '+numb+ '/'+len(listSim))
			# cF = open(currFile, 'a')
			# cF.write(sim+'\n\n')
			# cF.close()
			# print(sim)
			l, w, a, p = LastTimePoint(sim, currFile, listOutput=[], sevSim=True)
			if overTime:
				GraphsOverTime(sim, currFile)
			Id.append(sim.split('/')[-1])
			Lengths.append(l)
			Widths.append(w)
			Areas.append(a)
			Perimeters.append(p)
			numb+=1
	# histograms showing the distribution of data over all the simulations
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
	nbins = 20
	ax1.hist(Lengths, nbins, histtype='step', cumulative=True)
	ax1.set_title('Lengths')
	ax2.hist(Widths, nbins, histtype='step', cumulative=True)
	ax2.set_title('Widths')
	ax3.hist(Areas, nbins, histtype='step', cumulative=True)
	ax3.set_title('Areas')
	ax4.hist(Perimeters, nbins, histtype='step', cumulative=True)
	ax4.set_title('Perimeters')
	fig.savefig(foldername + '/Variability_Histograms.png')
	plt.close(fig)
	# writing of the data for each simulation
	cF = open(currFile, 'a')
	cF.write('Id	Length	Width	Area	Perimeter\n')
	for p in range(len(Lengths)):
		cF.write(Id[p]+'	'+str(Lengths[p])+'	'+str(Widths[p])+'	'+str(Areas[p])+'	'+str(Perimeters[p])+'\n')
	cF.write('\n\n\n\n')
	cF.write ('Length:	'+ str(np.mean(np.array(Lengths)))+ '	'+str(np.std(np.array(Lengths))) +'\n')
	cF.write ('Width:	'+ str(np.mean(np.array(Widths)))+ '	'+str(np.std(np.array(Widths))) +'\n')
	cF.write ('Area:	'+ str(np.mean(np.array(Areas)))+ '	'+str(np.std(np.array(Areas))) +'\n')
	cF.write ('Perimeter:	'+ str(np.mean(np.array(Perimeters)))+ '	'+str(np.std(np.array(Perimeters))) )
	cF.close()

def PlotTime(simname, Step, Data, savename, col='black'):
	# Data can either be a float (length, width, etc) or a couple [mean, sd] when a distribution is displayed (Elast)
	if not(isinstance(Data[0],float)):
		# [mean, sd] case
		Data = zip(*Data)
		Mean = np.array(Data[0])
		Sd = np.array(Data[1])
	else:
		# only one data case
		Mean = Data
		Sd = []
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	try: # the color is asked in the Output_names.txt file. If the color is not recognized, then balck is used
		col = plt.colors.to_rgb(col)
	except:
		col = 'black'
	ax.plot(Step, Mean, color=col)
	# if there is a standard deviation specified
	if not(Sd==[]):
		ax.plot(Step, Mean-Sd,'--')
		ax.plot(Step, Mean+Sd,'--')
	# Save the plot
	fig.savefig(simname+'/Plot/Analyses/'+savename+'.png')
	plt.close(fig)

def getStepfromName(filename):
	fspl = filename.split('__')
	return int(fspl[-2])

def main(name, opt):
	'''
	Main function. From the input argument, decides if there is a single time step, a simulation or a set of simulations
	'''
	# is the name the name of a file or a simulation or a set of simulations?
	if name[-4:] == '.txt': # Only one time step
		l = name.split('/')
		n = l[-1] # name of the file (without the path)
		l = l[:-2]
		l = '/'.join(l) # path
		currFile = l+'/Log_'+n+'.txt'
		cF = open(currFile,'a')
		cF.write(n+'\n\n\n')
		cF.close()
		# we only have one file
		OnePoint(name, currFile)
	else: # Folder
		l = name.split('/')
		currFile = name+'/Log.txt'
		cF = open(currFile,'a')
		cF.write(l[-1]+'\n\n\n')
		cF.close()
		# -> is it one simulation or many?
		if 'Data' in os.listdir(name): # every simulation has a Data folder = this case is one simulation
			print('Analyzing the simulation '+l[-1]+'...')
			# Dealing with options
			if 'l' in opt:
				LastTimePoint(name, currFile) # l = only the last point
			elif 'a' in opt: # a = all time points
				AllTimePoints(name, currFile)
			else:
				SevTimePoints(name, currFile) # default = several time points
			if 't' in opt:
				GraphsOverTime(name, currFile)
		else: # this is one set of simulations
			print('Analyzing the simulations in the folder '+l[-1]+'...')
			# several simulations: default = last time point for each simulation
			SevSim(name, currFile, overTime='t' in opt)	# t = with the graphs over time

if __name__ == "__main__":
	if len(sys.argv) > 2:
		name = sys.argv[2]
		opt = sys.argv[1]
	else:
		name = sys.argv[1]
		opt = ''
	main(name, opt)
