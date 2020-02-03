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


from class_point import *
import csv
import random
import pylab
from tvtk.api import tvtk
from mayavi.api import OffScreenEngine
import matplotlib.pyplot as plt



class mesh:
	def __init__(self):
		self.listPoints = [] # list containing point objects
		self.prefix = None # Prefix of the mesh (= step in the simulation)
		self.outputFolder = None # Where to save plots
		self.subList = [] # smaller list of points
		self.listIndexes = [] # list of only the indexes of the points in the listPoints list

	def getPrefix(self):
		return self.prefix

	def setPrefix(self, s):
		self.prefix = s

	def getOutputFolder(self):
		return self.outputFolder

	def setOutputFolder(self,s):
		self.outputFolder = s

	def getList(self):
		return self.listPoints

	def setListIndexes(self):
		for i in self.getList():
			self.listIndexes.append(i.getIndex())

	def getListIndexes(self):
		if self.listIndexes == []:
			self.setListIndexes()
		return self.listIndexes

	def getnbTotP(self):
		return len(self.getList())

	def getNbTotTr(self): # total number of triangles
		TrMax = 0
		for Point in self.getList():
			if Point.getnbTr()>TrMax:
				TrMax = Point.getnbTr()
		return TrMax+1

	def addPoint(self, NewPoint):
		self.listPoints.append(NewPoint)

	def readFile(self, filename):
		# reads the txt file containing the mesh info
		csvin = csv.reader(open(filename,'r'), delimiter=' ')
		csvin.next()
		for row in csvin:
			NewPoint = point()
			NewPoint.initiate(row)
			self.addPoint(NewPoint)
		self.findFilenameInfo(filename)

	def getScalar(self, scalarName, subList = False):
		'''
		Returns a list of the value of scalarName for each point in the mesh - usable by plotScalar
		'''
		data = []
		Index = []
		if subList:
			pointList = self.getSubList()
		else:
			pointList = self.getList()
		for Point in pointList:
			if not (Point.getIndex() in Index): # takes care of the repetitions
				Index.append(Point.getIndex())
				data.append(eval('Point.get'+scalarName+'()'))
		return Index, data

	def getScalarByTriangle(self, TrIndex, Index, IndexS, Scalar):
		'''
		TrIndex, Index = output of getTriangles
		IndexS, Scalar = output of a function that outputs a scalar
		'''
		TrMax = self.getNbTotTr()
		Z = np.zeros((TrMax,3)) # list of triplets
		for i in range(len(TrIndex)): # for each triangle
			for j in range(len(Index[0])): # For each point of the triangle
				# Ind = index of the vertex
				Ind = Index[i,j]
				# Get where in the Scalar array this index is
				IndS = IndexS.index(Ind)
				# Get the Scalar value
				Sc = Scalar[IndS]
				# Store
				Z[i,j] = Sc
		return Z

	def getCoords(self):
		X=[]
		Y=[]
		Index = []
		for Point in self.getList():
			if not(Point.getIndex() in Index):
				Index.append(Point.getIndex())
				X.append(Point.getCoords()[0])
				Y.append(Point.getCoords()[1])
		return Index, [X, Y]

	def getTriangles(self):
		'''
		Returns a list of the triangles for the plot, the list of the vertexes in each traingle and 2 (x and y) lists of triplets
		'''
		TrMax = self.getNbTotTr()
		X = np.zeros((TrMax,3))
		Y = np.zeros((TrMax,3))
		Index = np.zeros((TrMax,3))-1 # initialization at -1 because there is a vertex nb 0
		TrIndex = np.zeros((TrMax))-1
		for Point in self.getList():
			if not(Point.getnbTr() in TrIndex): # that triangle has not been started yet
				# find where in the array I need to go (first space with a -1)
				i = np.where(TrIndex==-1)[0][0]
				TrIndex[i]=Point.getnbTr()
				X[i,0] = Point.getCoords()[0]
				Y[i,0] = Point.getCoords()[1]
				Index[i,0] = Point.getIndex()
			else:
				# The triangle has already been started
				# find the triangle
				I =np.where(TrIndex==Point.getnbTr())[0][0]
				# and where to fill in this triangle
				i = np.where(Index[I,:]==-1)[0][0]
				X[I,i] = Point.getCoords()[0]
				Y[I,i] = Point.getCoords()[1]
				Index[I,i] = Point.getIndex()
		return TrIndex, Index, [X, Y]

	def getBigTriangles(self):
		'''
		Returns a list of the triangles after growth (ie +ux, uy) for
		the plot (2 (x and y) lists of triplets)
		'''
		TrMax = self.getNbTotTr()
		X = np.zeros((TrMax,3))
		Y = np.zeros((TrMax,3))
		Index = np.zeros((TrMax,3))-1 # initialization at -1 because there is a vertex nb 0
		TrIndex = np.zeros((TrMax))-1
		for Point in self.getList():
			if not(Point.getnbTr() in TrIndex): # that triangle has not been started yet
				# find where in the array I need to go (first space with a -1)
				i = np.where(TrIndex==-1)[0][0]
				TrIndex[i]=Point.getnbTr()
				x = Point.getCoords()[0]
				y = Point.getCoords()[1]
				ux = Point.getDispl()[0]
				uy = Point.getDispl()[1]
				X[i,0] = x+ux
				Y[i,0] = y+uy
				Index[i,0] = Point.getIndex()
			else:# The triangle has already been started
				# find the triangle
				I =np.where(TrIndex==Point.getnbTr())[0][0]
				# and where to fill in this triangle
				i = np.where(Index[I,:]==-1)[0][0]
				x = Point.getCoords()[0]
				y = Point.getCoords()[1]
				ux = Point.getDispl()[0]
				uy = Point.getDispl()[1]
				X[i,0] = x+ux
				Y[i,0] = y+uy
				Index[I,i] = Point.getIndex()
		return TrIndex, Index, [X, Y]

	def getGrTrArea(self):
		'''
		returns the input for the display of areal growth
		TrIndex, Index, [X,Y], Z where Z is the areal growth (Aire2/Aire1)
		X, Y = One dimensional arrays - one value for each vertex
		'''
		TrMax = self.getNbTotTr()
		VertMax = self.getnbTotP()
		X = np.zeros((VertMax))
		Y = np.zeros((VertMax))
		Z = np.zeros((TrMax))-1
		Index = np.zeros((TrMax,3))-1
		TrIndex = np.zeros((TrMax))-1
		TrIndexS, IndexS, [XS, YS] = self.getTriangles()
		TrIndexB, IndexB, [XB, YB] = self.getBigTriangles()
		for i in range(len(TrIndexS)): # for each triangle
			TrIndex[i] = TrIndexS[i]
			Index[i] = IndexS[i]
			for j in range(len(IndexS[i])): # storage of the small triangles coordiantes for the plot
				X[IndexS[i][j]] = XS[i][j]
				Y[IndexS[i][j]] = YS[i][j]
			AS = TrArea(XS[i],YS[i])
			# find the good triangle in TrIndexB
			ind = np.where(TrIndexB==TrIndex[i])[0][0]
			AB = TrArea(XB[ind],YB[ind])
			Z[i] = AB/AS
		return TrIndex, Index, [X, Y], Z

	def getEdges(self):
		'''
		Returns a list of all the edges:
		returns a list of couples of points (index of points) and a list of triangles in which these couples of points are
		'''
		TrIndex, Index, [X,Y] = self.getTriangles()
		Edges = []
		TrEInd = []
		for i in range(len(TrIndex)): # for each triangle
			# get the 3 edges
			e = Index[i]
			# arrange the vertices indexes in the 3 edges for the unicity ([v1, v2] is the same edge as [v2, v1])
			EdgTr = [[min(e[0],e[1]),max(e[0],e[1])],[min(e[0],e[2]),max(e[0],e[2])],[min(e[2],e[1]),max(e[2],e[1])]]
			for e in EdgTr:
				# for each couple, check if already in the list Edges
				if e in Edges: # yes -> add the new triangle to the already sorted edge
					ind = Edges.index(e)
					TrEInd[ind].append(TrIndex[i])
				else: # no -> new entry in the Edges list
					Edges.append(e)
					TrEInd.append([TrIndex[i]])
		return Edges, TrEInd

	def saveContour(self):
		'''
		Main function to save the contour data (for the contour variability analysis)
		'''
		Index, nbTr = self.getBoundaryInd() # from the points that are already boundary
		# list of indexes and triangle id -> each triangle number is there twice (two vertexes), and
		# eahc vertex index is there twice (belongs to 2 triangles)
		for i in range(2,len(Index)):
			# 2 first points are ok - find either the next vertex in the same triangle, or the next triangle for the same vertex
			# series of permutations of the original lists to have stg like
			# [v1, v2, v2, v3, v3, v4,...]
			# [t1, t1, t2, t2, t3, t3,...]
			if nbTr[i-1]==nbTr[i-2]: # want to find the next triangle for the same vertex
				j=i
				while Index[i]!=Index[i-1]:
					Index[i],Index[j] = inter(Index[i],Index[j])
					nbTr[i], nbTr[j] = inter(nbTr[i],nbTr[j])
					j+=1
			else: # want to find the next vertex for the same triangle
				j=i
				while nbTr[i]!=nbTr[i-1]:
					Index[i],Index[j] = inter(Index[i],Index[j])
					nbTr[i], nbTr[j] = inter(nbTr[i],nbTr[j])
					j+=1
		# simplify the vertex list to remove repetitions
		Index = [Index[2*i] for i in range(len(Index)/2)]
		# get coordinates
		IndexC, [XC, YC] = self.getCoords()
		X = []
		Y = []
		for p in Index:
			X.append(XC[IndexC.index(p)])
			Y.append(YC[IndexC.index(p)])
		# find the lowest point
		d = 10.
		ind = None
		for p in range(len(Index)):
			if (X[p]**2+Y[p]**2)<d:
				d = (X[p]**2+Y[p]**2)
				ind = p
		# rearrange
		Index = rearrange(Index, ind)
		X = rearrange(X, ind)
		Y = rearrange(Y, ind)
		# check anticlockwise
		if (X[1]-X[0]<0):	# correct if not
			Index = Index[::-1]
			X = X[::-1]
			Y = Y[::-1]
		# save
		name = str(self.getOutputFolder())
		name = name.split('/')
		outputFold = str(self.getOutputFolder())+'../../../'+name[-4]
		cFF = open(outputFold+'_Contourxy.vx','w')
		cFF.write('# 6\n')
		cFF.write('# 2\n')
		cFF.write('# '+str(len(Index))+'\n')
		for p in range(len(Index)):
			cFF.write(str(X[p])+' '+str(Y[p])+'\n')
		cFF.close()

	def getBoundaryInd(self):
		'''
		Returns a list of vertexes on the edges and their triangles
		(an edge is on the border if it is not shared with another triangle)
		'''
		Edges, TrEInd = self.getEdges()
		Index = []
		nbTr = []
		for i in range(len(Edges)):
			if len(TrEInd[i]) == 1: # this edges is in one triangle = on the boundary
				Index.extend(Edges[i])
				nbTr.extend(TrEInd[i])
				nbTr.extend(TrEInd[i])
				# nbTr.extend twice because 2 point for each triangle
		return Index, nbTr

	def getBoundaryCoords(self):
		'''
		returns imbricated lists centered on the triangle
		TrB = [t1, t2, t3,...]
		IndB = [[i1,i2,i3],[i4,i5,i6],...]
		XB = [[xi1,xi2,xi3],[xi4,xi5,xi6],...]
		'''
		IndB, nbTr = self.getBoundaryInd()
		# recup X, Y
		Indall, [X, Y] = self.getCoords()
		XB = []
		YB = []
		IB = []
		TrB = []
		for i in range(len(IndB)):
			indb = IndB[i]
			ind = Indall.index(indb)
			if not(nbTr[i] in TrB):
				TrB.append(nbTr[i])
				XB.append([X[ind]])
				YB.append([Y[ind]])
				IB.append([indb])
			else:
				indt = TrB.index(nbTr[i])
				XB[indt].append(X[ind])
				YB[indt].append(Y[ind])
				IB[indt].append(indb)
		return TrB, IB, [XB, YB]

	def findFilenameInfo(self, filename):
		filenameCut = filename.split('/')
		filenameEnd = filenameCut[-1]
		self.setPrefix(filenameEnd[:-8])
		ind = filenameCut.index('Data')
		sep = '/'
		filenameS = filenameCut[:ind]
		filenameS.extend(['Plot','Analyses'])
		outputFold = sep.join(filenameS)
		self.setOutputFolder(outputFold+'/')

	def plotMesh(self, nbTr, Index, Coords, savename = 'mesh'): # Coords = [X, Y]
		'''
		Plots a mesh
		'''
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		self.plotaxMesh(ax, nbTr, Index, Coords)
		title = self.getPrefix()
		title = title.split('__')
		title = int(title[1])
		ax.set_title(str(title))
		# save
		fig.savefig(str(self.getOutputFolder())+str(self.getPrefix())+'_'+savename+'.png')
		plt.close(fig)

	def plotaxMesh(self, ax, nbTr, Index, Coords): # Coords = [X, Y]
		n = len(Coords[0][0])
		if n==2: # for an edge
			X = np.zeros(n)
			Y = np.zeros(n)
		else: # for a triangle (to close the triangle, need the repetition of the first point at the end)
			X = np.zeros(n+1)
			Y = np.zeros(n+1)
		for t in range(len(nbTr)):
			if len(Coords[0][t])==n:
				X[0:n] = Coords[0][t]
				Y[0:n] = Coords[1][t]
			if n!=2:
				X[n] = Coords[0][t][0]
				Y[n] = Coords[1][t][0]
			ax.plot(X, Y, color = 'k')
		ax.axis('equal')

	def convertToPolydata(self, IndexS):
		'''
		Function used for plotScalar - use Mayavi library
		IndexS and Scalar must be numpy arrays
		'''
		nbP = len(IndexS)
		nbT = self.getNbTotTr()
		xyzPoints = np.zeros((nbP,3))
		trianglesI = np.zeros((nbT, 3))-1
		TrIndex = np.zeros(nbT)-1
		for Point in self.getList():
			# find where in the xyzPoints I need to put the xy data, ie the position of this point in IndexS
			i = np.where(IndexS==Point.getIndex())[0][0]
			xyzPoints[i,0] = Point.getCoords()[0]
			xyzPoints[i,1] = Point.getCoords()[1]
			# is this triangle already in the trianglesI list?
			if not(Point.getnbTr() in TrIndex):
				# find where in the array I need to go (first space with a -1)
				j = np.where(TrIndex==-1)[0][0]
				TrIndex[j] = Point.getnbTr()
				trianglesI[j,0] = i
			else:
				j = np.where(TrIndex == Point.getnbTr())[0][0]
				# find where to fill this triangle
				k = np.where(trianglesI[j,:]==-1)[0][0]
				trianglesI[j,k] = i
		return xyzPoints, trianglesI

	def plotScal(self, IndexS, Scalar, savename, square=True, sevSim=False, colormap = 'jet', colorrange=[None,None]): # 'Greys'
		'''
		Input: Output of a fonction that outputs a scalar (with the associated index of the vertices)
		plots a mesh filled with colors corresponding to a scalar value
		Possible colormaps: (to add in the Output_names.txt file)
		'Accent' or 'Blues' or 'BrBG' or 'BuGn' or 'BuPu' or 'Dark2' or 'GnBu' or 'Greens' or 'Greys' or 'OrRd' or 'Oranges' or 'PRGn' or 'Paired' or 'Pastel1' or 'Pastel2' or 'PiYG' or 'PuBu' or 'PuBuGn' or 'PuOr' or 'PuRd' or 'Purples' or 'RdBu' or 'RdGy' or 'RdPu' or 'RdYlBu' or 'RdYlGn' or 'Reds' or 'Set1' or 'Set2' or 'Set3' or 'Spectral' or 'YlGn' or 'YlGnBu' or 'YlOrBr' or 'YlOrRd' or 'autumn' or 'binary' or 'black-white' or 'blue-red' or 'bone' or 'cool' or 'copper' or 'file' or 'flag' or 'gist_earth' or 'gist_gray' or 'gist_heat' or 'gist_ncar' or 'gist_rainbow' or 'gist_stern' or 'gist_yarg' or 'gray' or 'hot' or 'hsv' or 'jet' or 'pink' or 'prism' or 'spectral' or 'spring' or 'summer' or 'winter'
		'''
		from mayavi.sources.vtk_data_source import VTKDataSource
		from mayavi.modules.surface import Surface
		import mayavi.mlab
		# size of a transparent shape to have the same scale on all images
		xmS = -14.
		xMS = 14.
		ymS = 20.
		yMS = 26.
		square = [[xmS,xmS,xMS,xMS],[ymS,yMS,yMS,ymS]]
		xyzSquare = [[xmS,ymS,0],[xmS,yMS,0],[xMS,yMS,0],[xMS,ymS,0]] # because Mayavi is in 3D
		trSquare = [[0,1,2],[0,2,3]]
		IndexS = np.array(IndexS)
		Scalar = np.array(Scalar)
		meshsq = tvtk.PolyData(points=xyzSquare, polys=trSquare)
		meshsq.point_data.scalars = [0,0]
		meshsq.point_data.scalars.name = 'square'
		xyzPoints, trianglesI = self.convertToPolydata(IndexS)
		# The TVTK dataset.
		mesh = tvtk.PolyData(points=xyzPoints, polys=trianglesI)
		mesh.point_data.scalars = Scalar
		mesh.point_data.scalars.name = savename
		# Create the MayaVi offscreen engine and start it.
		e = OffScreenEngine()
		e.start()
		# Create a new scene.
		win = e.new_scene()
		src = VTKDataSource(data = mesh)
		src2 = VTKDataSource(data = meshsq)
		e.add_source(src2)
		sSq = Surface()
		e.add_module(sSq)
		e.add_source(src)
		s = Surface()
		e.add_module(s)
		# to plot a transparent square, I get the lut table for the square (= where mayavi looks for the colors) and I change the last column (alpha channel) to always transparent (=0)
		lut = sSq.module_manager.scalar_lut_manager.lut.table.to_array()
		lut[:, -1] = np.zeros(256)
		# and then give back the new table to the manager
		sSq.module_manager.scalar_lut_manager.lut.table = lut
		s.module_manager.scalar_lut_manager.lut_mode = colormap
		# change here if you want to remove the scale
		s.module_manager.scalar_lut_manager.show_scalar_bar = True
		s.module_manager.scalar_lut_manager.show_legend = True
		#~ s.module_manager.scalar_lut_manager.show_scalar_bar = False
		#~ s.module_manager.scalar_lut_manager.show_legend = False
		if colormap == 'Greys':
			s.scene.background = (1.0, 1.0, 1.0) # set background to white
		title = self.getPrefix()
		title = title.split('__')
		title = int(title[1])
		win.scene.camera.distance = 10.
		win.scene.z_plus_view()
		#~ mayavi.mlab.title(str(title)) # Uncomment this if you want the number of the step as a title
		if not not(colorrange[1]):
			s.module_manager.scalar_lut_manager.data_range = np.array([colorrange[0],  colorrange[1]])
		if sevSim:
			name = str(self.getOutputFolder())
			name = name.split('/')
			outputFold = str(self.getOutputFolder())+'../../../'+name[-4]
			win.scene.save(outputFold+'_'+savename+'.png', size=(1200, 600))
		else:
			win.scene.save(str(self.getOutputFolder())+str(self.getPrefix())+'_'+savename+'.png', size=(800, 800))
		win.scene.close()
		e.stop()

	def getSubList(self):
		return self.subList

	def setSubList(self):
		'''
		Returns a list of points with approx 1 over 3 points (1 point per triangle kept)
		'''
		listAll = []
		vertInd = []
		trInd = []
		subVertInd = []
		removedVertInd = []
		# first get the list of all points
		for Point in self.getList():
			if not(Point.getIndex() in vertInd):
				vertInd.append(Point.getIndex())
				trInd.append(Point.getnbTr())
				listAll.append(Point)
		indivTr = np.unique(np.array(trInd))
		for tr in indivTr: # for each triangle
			# find the vertices of the triangle
			c = []
			for i in range(len(trInd)):
				if trInd[i]==tr:
					c.append(i)
			# keep only one
			a = random.randrange(len(c))
			subVertInd.append(c[a])
		self.subList = []
		for i in subVertInd:
			self.subList.append(listAll[i])

	def plotScalTr(self, TrIndex, Index, Coords, Z, savename, sevSim=False, col = '', limScale = [None, None]): #Coords = [X,Y]
		fig = pylab.figure()
		ax = fig.add_subplot(1,1,1)
		if not(sevSim):
			TrB, IB, [XB, YB] = self.getBoundaryCoords()
			CoordsB = [XB, YB]
			self.plotaxMesh(ax, TrB, IB, CoordsB)
		#  !! If you want another colormap, change here
		if col == '' or col == 'jet':
			colmap = plt.cm.rainbow
		else:
			print 'Please modify the code in plotScalTr function to take you colormap into account'
			raise
		# Plot
		pylab.tripcolor(Coords[0], Coords[1], Index, vmin = limScale[0], vmax = limScale[1], facecolors=Z, edgecolors='none', cmap=colmap)
		if not(sevSim):
			title = self.getPrefix()
			title = title.split('__')
			title = int(title[1])
			ax.set_title(str(title))
		ax.axis('equal')
		# Axes range
		ax.set_xlim(-20.,20.)
		ax.set_ylim(-5.,25.)
		if sevSim:
			ax.axis('off')
		if not(sevSim):
			pylab.colorbar()
		# Save
		if sevSim:
			name = str(self.getOutputFolder())
			name = name.split('/')
			outputFold = str(self.getOutputFolder())+name[-4]
			fig.savefig(outputFold+'_'+savename+'.png',dpi=400)
		else:
			fig.savefig(str(self.getOutputFolder())+str(self.getPrefix())+'_'+savename+'.png',dpi=400)
		pylab.close(fig)

	def plotDistr(self, Scalar, savename, currFile):
		'''
		Plots a histogram
		'''
		S = np.array(Scalar)
		# save the distribution properties in the Log.txt file
		cF = open(currFile,'a')
		cF.write(savename+'\n')
		cF.write('Distribution properties:    mean;    sd\n')
		cF.write(str(np.mean(S))+"; "+str(np.std(S))+'\n\n\n')
		cF.close()
		fig = pylab.figure()
		ax = fig.add_subplot(1,1,1)
		ax.hist(Scalar, 50)
		# add a title
		title = self.getPrefix()
		title = title.split('__')
		title = int(title[1])
		ax.set_title(str(title))
		fig.savefig(str(self.getOutputFolder())+str(self.getPrefix())+'_'+savename+'.png')
		pylab.close(fig)

	def getLength(self):
		'''
		Returns the difference between highest and lowest point
		'''
		ymin = None
		ymax = None
		for Point in self.getList():
			y = Point.getCoords()[1]
			if (not ymin) or (y<ymin):
				ymin = y
			if (not ymax) or (y > ymax):
				ymax = y
		return ymax - ymin

	def getWidth(self):
		'''
		Returns the difference between the most on the left and the most on the right
		'''
		xmin = None
		xmax = None
		for Point in self.getList():
			x = Point.getCoords()[0]
			if (not xmin) or (x<xmin):
				xmin = x
			if (not xmax) or (x > xmax):
				xmax = x
		return xmax - xmin

	def getxyminmax(self):
		xmin = None
		xmax = None
		ymin = None
		ymax = None
		for Point in self.getList():
			y = Point.getCoords()[1]
			if (not ymin) or (y<ymin):
				ymin = y
			if (not ymax) or (y > ymax):
				ymax = y
			x = Point.getCoords()[0]
			if (not xmin) or (x<xmin):
				xmin = x
			if (not xmax) or (x > xmax):
				xmax = x
		return xmin, xmax, ymin, ymax

	def getPerim(self):
		'''
		Returns the sum of the length of the edges in the boundary
		'''
		Perim = 0
		TrB, IB, [XB, YB] = self.getBoundaryCoords()
		for i in range(len(XB)):
			x = XB[i]
			y = YB[i]
			Perim += np.sqrt((x[0]-x[1])**2+(y[0]-y[1])**2)
		return Perim

	def getArea(self):
		'''
		Returns the sum of the areas of all the triangles
		'''
		Area = 0
		TrIndex, Index, [X, Y] = self.getTriangles()
		for i in range(len(X)):
			x1 = X[i][1]-X[i][0]
			x2 = X[i][2]-X[i][0]
			y1 = Y[i][1]-Y[i][0]
			y2 = Y[i][2]-Y[i][0]
			Area += abs(x1*y2-x2*y1)/2.0
		return Area

	def getDistrParam(self, string):
		_, Scalar = self.getScalar(string, subList = False)
		Scalar = np.array(Scalar)
		return [np.mean(Scalar), np.std(Scalar)]

	def findVert(self, ind):
		# finds the Point object from its index
		listI = self.getListIndexes()
		listV = self.getList()
		return listV[listI.index(ind)]

def TrArea(X, Y):
        '''
        Returns the area of the triangle of coordinates X = [x1,x2,x3], Y = [y1,y2,y3]
        '''
        v1 = [0,0]
        v2 = [0,0]
        v1[0] = X[1]-X[0]
        v1[1] = Y[1]-Y[0]
        v2[0] = X[2]-X[0]
        v2[1] = Y[2]-Y[0]
        return abs(v1[0]*v2[1]-v1[1]*v2[0])/2.0
