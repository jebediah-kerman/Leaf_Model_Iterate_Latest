####################################################################
#                                                                  #
#                                                                  #
#                     ANALYSIS OF MODEL OUTPUT                     #
#                                                                  #
#  author: Mathilde Dumond                                         #
#                                                                  #
#  Supplemental Information of the paper:                          #
#  Variable cell growth yields reproducible organ development      #
#  through spatiotemporal averaging                                #
#                                                                  #
#  Hong et al.                                                     #
#  Dev Cell 2016                                                   #
#                                                                  #
####################################################################



usage:
./main.py options filename_or_foldername
options can be any combination (string of several
                               characters) or none of:
         'a' if we want all data points of a simulation
                   to be analyzed (this can take a lot of time)
         'l' if we want only the last time point of a
                   simulation to be analyzed
         't' if we want a graph of a quantity over time
                   (see second frame)

filename_or_foldername is one of:
         text file output of the model (in the Data folder)
         folder containing one simulation
         folder containing a set of simulations
the options have no meaning if a text file is as input
and won't be addressed
the 'a' and 'l' options will not be used if a set of
simulations is given. for each simulation only the last
time point is analyzed.
if the 't' option is activated then overtime will be
computed for all simulations (! can be very time consuming)

the options 'a', 'l' and 't' are all receivable for a single
simulation, but if 'a' and 'l' are specified ('al'), only
'l' will be executed.

Change the file 'Input_names.txt' to specify which graphs are
to be plotted

An output file is written containing mean, etc info
named "Log.txt"


###################################################################

List of the functions:

## main.py
  - inter(a,b)
    Purpose: invert two objects
    Arguments: a and b: any kind of objects
    Outputs: a and b in the opposite order

  - rearrange(liste, ind)
    Purpose: change the first point of a list (but not the order)
    Arguments: liste: list: list that we want to change
               ind: index of the new first point
    Output: list: new list beginning at the ind index

  - ReadFile(filename)
    Purpose: creates a new mesh and initializes it from a file
    Argument: filename: string: name of the file where the information is
    Output: mesh: mesh object following the specifications in filename

  - tofloatSt(string, line)
    Purpose: When reading the "Output_names.txt" file, checks that the inputs in the Min and Max columns are floats
    Arguments: string: string in the Min or Max column for a given line
               line: line from where the string is taken (for the error output)
    Output: float: float of the corresponding number.
    Note: Raises an error if not a number

  - ReadAllNames()
    Purpose: Reads the file Output_names.txt
    Output: array of strings: array contining info from the Output_names file

  - cleanListStrings(l)
    Purpose: Removes all blank spaces from each string in the list of strings
    Argument: l: list of strings
    Output: list of strings: the same list with no space and no empty value

  - ReadToAnalyze()
    Purpose: Reads the file Input_names and gets the info
    Output: [[strings],[strings]] list of the two kind of output we want (one time point based, or overtime)

  - OneTimePoint(NewTime, currFile, sevSim=False)
    Purpose: Main function for the analysis of one data point
    Arguments: NewTime: mesh: mesh to be analyzed
               currFile: string: name of the Log file
               sevSim: bool: whether several simulations are being analyzed at once

  - findSpec(name, listAll)
    Purpose: Returns whether there is a 'D' and/or a 'LOCVAR' substring in the string and the canonical name
    Arguments: name: name of the wished output (from Input_names.txt)
               listAll: list of all possible outputs (from ReadAllNames())
    Outputs: name: string: short name of the characteristic
             ind: int: where in the listAll is the info
             title: string: long name of the characteristic
             typ: string: type of the characteristic (Scatar, Mesh, etc)
             col: string: color for the plot
             limScale: [float, float] color range for the plot
    Note: If you want to add other specifications do it here

  -  OnePoint(name, currFile)
    Purpose: Launches the OneTimePoint function from the name of the file containing the mesh info
    Arguments: name: string: name of the file containing mesh information
               currFile: string: name of the Log.txt file

  - AllTimePoints(simname, currFile)
    Purpose: Plot all the point to easily make a video using ImageJ
    Arguments: simname: string: path to the simulation to analyse
               currFile: string: name of the Log.txt file

  - SevTimePoints(simname, currFile, nbPts=[3,3])
    Purpose: Get a panel of different time (equally distributed) points showing the evolution of the simulation
    Arguments: simname: string: path to the simulation to analyse
               currFile: string: name of the Log.txt file
               nbPts: [int, int]: shape of the panel - nb of graphs in each direction
    Note: If you want more or less than 9 steps, please change the nbPts default argument

- LastTimePoint(simname, currFile, sevSim = False)
  Purpose: Analysis of only the last time point of a simulation
  Arguments: simname: string: path to the simulation to analyse
             currFile: string: name of the Log.txt file
             sevSim: bool: whether several simulations are being analyzed at once
  Outputs: float: length of the mesh (vertical)
           float: width of the mesh (horizontal)
           float: area of the mesh
           float: perimeter of the border of the mesh

  - MergeGraphs(iterations, simname, nbPts)
    Purpose: Finds all png files corresponding to the same data, merges them all together and deletes the original files. Deals with the iteration numbers so that if you did another step it is not erased
    Arguments: iterations: list of strings: list of the iterations to be in the merged figure
               simname: string: path to the simulation to analyse
               nbPts: [int, int]: shape of the panel - nb of graphs in each direction

  - SortByType(listFiles, iterations, simname)
    Purpose: From a list of filenames, gives back a dictionnary with as key the title of the kind of graph, and as entry the ordered list of all png files
    Argument: listFiles: list of strings: list of filenames where the step number and the name are separated by '___'
              iterations: list of strings: list of the iterations to be in the merged figure
              simanme: string: path to the simulation to analyse
    Output: dictionary: with as keys the title of the graphs (string), and values lists of the png filenames associated with this kind of graph

  - MergeAGraph(listFiles, nbPts, title, simname)
    Purpose: From a list of png files, makes a big png file with all of them
    Arguments: listFiles: list of strings: list of filenames
               nbPts: [int, int]: shape of the panel - nb of graphs in each direction
               title: string: title of the file
               simname: string: path to the simulation to analyse

  - InitializeDataFile(NewTime)
    Purpose: Initalizes the file in which are stored the over time data (xmin & co) - each simulation has one file
    Argument: NewTime: mesh: one mesh of the simulation (the first one in practice but not mandatory for this function)
    Output: string: name and path of the text file

  - GraphsOverTime(simname, currFile)
    Purpose: Main function for the graphs 'over time' - to have the evolution of a quantity over time
    Arguments: simname: string: path to the simulation to analyse
               currFile: string: name of the Log.txt file

  - DataOverTime(currOTFileName,t, NewTime)
    Purpose: Writes in a file the xmin, xmax, ymin and ymax info
    Arguments: currOTFileName: string: name of the txt file where the overtime info is written (initialized in InitializeDataFile)
               t: string: step of the simulation (time information)
               NewTime: mesh: one mesh of the simulation

  - SevSim(foldername, currFile, overTime=False)
    Purpose: Main function dealing with several simulations at once. For each simulation, launches the main function of one simulation
    Arguments: foldername: string: name and path of the folder containing all the simulations (the folder should contain only folders of simulations except for the 'Source' folder)
               currFile: string: name of the Log.txt file
               overTime: bool: whether the overtime plots and data saving will be done

  - PlotTime(simname, Step, Data, savename, col='black')
    Purpose: Plot of the overtime data (second list in Input_names.txt)
    Arguments: simname: string: path to the simulation to analyse
               Step: list of float: list of steps, time course
               Data: list of floats or list of couples: Data to be plotted (can be list of data, or list of [Mean, Sd])
               savename: string: name of the file to be saved
               col: string: color of the plot. Default is 'black'

  - getStepfromName(filename)
    Purpose: Returns the step of a file (depends on the format of the string)
    Argument: filename: string: name of the file containing mesh information
    Output: int: step number

  - main(name, opt)
    Purpose: Main organizer - calls the functions depending on the argument
    Arguments: name: string: path to a folder containing a set of simulations, or path to a single simulation or a txt file of only one step of a simulation
               opt: string: containing the options. For details see beginnning of this file

##################################

class_mesh.py

  - class mesh:

    - __init__(self)
      Purpose: Initialization of the object
      Attributes:
      self.listPoints = [] # list containing point objects
      self.prefix = None # Prefix of the mesh (= step in the simulation)
      self.outputFolder = None # Where to save plots
      self.subList = [] # smaller list of points
	  self.listIndexes = [] # list of only the indexes of the points in the listPoints list

	- getPrefix(self)
      Purpose: get the attribute: prefix

	- setPrefix(self, s)
	  Purpose: set the attribute: prefix
      Argument: s: string: containing simulation index and step in the simulation

	- getOutputFolder(self)
	  Purpose: get the attribute: outputFolder (where to save the plots)

	- setOutputFolder(self,s)
      Purpose: set the attribute: outputFolder (where to save the plots)
      Argument: s: string: path to the folder

	- getList(self)
      Purpose: get the attibute: listPoints (list of point objects)

	- setListIndexes(self)
      Purpose: set the attribute: listIndexes (list of indexes of the points in the same order as listPoints)

	- getListIndexes(self)
      Purpose: get the attribute: listIndexes (list of indexes of the points in the same order as listPoints)

	- getnbTotP(self)
      Purpose: get the number total of points in the mesh
      Output: int: nb of points in the mesh

	- getNbTotTr(self)
      Purpose: get the total number of triangles
      Output: int: total number of triangles

	- addPoint(self, NewPoint)
      Purpose: adds a new point to the list of points
      Argument: NewPoint: point: point to be added

	- readFile(self, filename)
	  Purpose: reads the txt file containing the mesh info, sets the listPoints, prefix, and outputFolder attributes
      Argument: filename: string: name of the file containing the mesh information

	- getScalar(self, scalarName, subList = False)
      Purpose: Returns a list of the value of scalarName for each point in the mesh - usable by plotScalar
      Arguments: scalarName: string: short name of a characteristic which is 'S' (from Output_names.txt)
                 subList: bool: whether we store all the points or only a subset
      Outputs: list of int: list of indexes of the vertices (unique)
               list of float: value of the scalar for each vertex

	- getScalarByTriangle(self, TrIndex, Index, IndexS, Scalar)
      Purpose: Organizes data to be readable by plotScalTr
      Arguments: TrIndex: list of int: list of triangle indexes. output of getTriangles. ex: [t1, t2, t3,...]
                 Index: list of triplets of int: list containing the point indexes for each triangle. output of getTriangles. ex: [[p1, p2, p3],[p2, p3, p4], [p7, p8, p9],...]
                 IndexS: list of int: list of point indexes. output of getScalar. ex: [p1, p2, p3, p4,...]
                 Scalar: list of float: values of scalar for each index. output of getScalar. ex: [s1, s2, s3, s4, ...]
      Output: list of triplets of float: list containing the scalar values  for each triangle. ex: [[s1, s2, s3],[s2, s3, s4], [s7, s8, s9],...]

	- getCoords(self)
      Purpose: Returns the coordinates of the vertices
      Outputs: list of int: list of point indexes. ex: [p1, p2, p3, p4,...]
               list of 2 lists of floats: X, Y coordinates of the vertices. ex: [[x1, x2, x3, x4,...], [y1, y2, y3, y4, ...]]

    - getTriangles(self)
      Purpose: Returns a list of the triangles for the plot, the list of the vertexes in each triangle and 2 (x and y) lists of triplets
      Outputs: list of int: list of triangle indexes. output of getTriangles. ex: [t1, t2, t3,...]
               list of triplets of int: list containing the point indexes for each triangle. output of getTriangles. ex: [[p1, p2, p3],[p2, p3, p4], [p7, p8, p9],...]
               list of 2 lists of triplets of floats: X, Y coordinates of the vertices. ex: [[[x1, x2, x3],[x2, x3, x4] ,...], [[y1, y2, y3], [y2, y3, y4] , ...]]

	- getBigTriangles(self)
      Purpose: Returns a list of the triangles after growth (ie +ux, uy) for the plot and 2 (x and y) lists of triplets
      Outputs: list of int: list of triangle indexes. output of getTriangles. ex: [t1, t2, t3,...]
               list of triplets of int: list containing the point indexes for each triangle. output of getTriangles. ex: [[p1, p2, p3],[p2, p3, p4], [p7, p8, p9],...]
               list of 2 lists of triplets of floats: X, Y coordinates of the vertices with the deformation. ex: [[[x1, x2, x3],[x2, x3, x4] ,...], [[y1, y2, y3], [y2, y3, y4] , ...]]

	- getGrTrArea(self)
      Purpose: Returns the input for the display of areal growth
      Outputs: list of int: list of triangle indexes. output of getTriangles. ex: [t1, t2, t3,...]
               list of triplets of int: list containing the point indexes for each triangle. output of getTriangles. ex: [[p1, p2, p3],[p2, p3, p4], [p7, p8, p9],...]
               list of 2 lists of triplets of floats: X, Y coordinates of the vertices with the deformation. ex: [[[x1, x2, x3],[x2, x3, x4] ,...], [[y1, y2, y3], [y2, y3, y4] , ...]]
               list of float: list of triangle values. ex: [z1, z2, z3,...]

	- getEdges(self)
      Purpose: Returns an list of all the edges: returns a list of couples of points (index of points) and a list of triangles in which these couples of points are
      Outputs: list of couples of int: list containing the edges of the mesh. ex: [[p1, p2],[p2, p3],[p4,p5],...]
               list of int: list of the triangles in which each edge is. ex: [t1, t1, t1, t2, t2, ...]

	- saveContour(self)
      Purpose: Main function to save the contour data (for the contour variability analysis) of name: outputFold+'_Contourxy.vx'

	- getBoundaryInd(self)
      Purpose: Returns a list of vertexes on the border of the mesh and their triangles (an edge is on the border if it is not shared with another triangle)
      Outputs: list of int: list of indexes of the vertices that are on the border of the mesh
               list of int: list of indexes of the triangles corresponding to each vertex in the first list

	- getBoundaryCoords(self)
      Purpose: Returns the coordinates of the vertices that are on the border of the mesh
      Outputs: list of int: list of indexes of triangles on the border. ex: [t1, t2, t3,...]
               list of couples of int: list containing the indexes of the vertices of the edges on the border of the mesh. ex: [[p1, p2],[p2, p3],[p4,p5],...]
               list of 2 lists of couples of floats: X, Y coordinates of the vertices on the border. ex: [[[x1, x2],[x2, x3] ,...], [[y1, y2], [y2, y3] , ...]]

	- findFilenameInfo(self, filename)
      Purpose: Extract useful information from the name of the file. Fills the prefix and outputFold attributes
      Argument: filename: string: name of the file defining the mesh

	- plotMesh(self, nbTr, Index, Coords, savename = 'mesh')
      Purpose: Does the plot of the whole mesh (all the triangles)
      Arguments: nbTr: int: total number of triangles
                 Index: list of couples of triplets of int: list containing the indexes of the vertices of the edges on the border of the mesh (couples), or on the triangles of the mesh (triplets)
                 Coords: list of 2 lists of couples of triplets of floats: each sublist contains the coordinates of the vertices (X and Y resp.) of the edges on the border of the mesh (couples), or on the triangles of the mesh (triplets)
                 savename: string: name of the saved file

	- plotaxMesh(self, ax, nbTr, Index, Coords)
      Purpose: Helper of the plotMesh function. Plots a subpart of the mesh
      Arguments: ax: figure (matplotlib object) on which to plot
                 nbTr: int: total number of triangles
                 Index: list of couples of triplets of int: list containing the indexes of the vertices of the edges on the border of the mesh (couples), or on the triangles of the mesh (triplets)
                 Coords: list of 2 lists of couples of triplets of floats: each sublist contains the coordinates of the vertices (X and Y resp.) of the edges on the border of the mesh (couples), or on the triangles of the mesh (triplets)

	- convertToPolydata(self, IndexS)
      Purpose: Helper of the plotScalar function. Converts data to readable format for plotScalar function
      Arguments: IndexS: list of int: list of point indexes. ex: [p1, p2, p3, p4,...]
	  Outputs: list of triplets: xyz coordinates of each vertex. z=0 in this case. ex: [[x1,y1,0],[x2,y2,0],[x3,y3,0],...]
               list of triplets: indexes of vertices in each triangle (in the order of the previous list). ex: [[p1,p2,p3],[p4,p5,p6],...]

	- plotScal(self, IndexS, Scalar, savename, square=True, sevSim=False, colormap = 'jet', colorrange=[None,None]): # 'Greys'
      Purpose: Does the plot of a scalar value on the mesh based on vertices
      Arguments: IndexS: list of int: list of point indexes. ex: [p1, p2, p3, p4,...]
                 Scalar: list of float: list of value for each point. ex: [s1, s2, s3, s4,...]
                 savename: string: name of the saved file
                 square: bool: whether a transparent square will be used to scale all the plots (for them to have the same scale)
                 sevSim: bool: whether the program is currently analyzing several simulations - or just one
                 colormap: string: kind of colorscale to use. Details in the 'Notes' part
                 colorrange: [float, float]: min and max of the color range. Default is [None, None] and in this case adapts to the data
      Notes: Possible colormaps: (to add in the Output_names.txt file): 'Accent' or 'Blues' or 'BrBG' or 'BuGn' or 'BuPu' or 'Dark2' or 'GnBu' or 'Greens' or 'Greys' or 'OrRd' or 'Oranges' or 'PRGn' or 'Paired' or 'Pastel1' or 'Pastel2' or 'PiYG' or 'PuBu' or 'PuBuGn' or 'PuOr' or 'PuRd' or 'Purples' or 'RdBu' or 'RdGy' or 'RdPu' or 'RdYlBu' or 'RdYlGn' or 'Reds' or 'Set1' or 'Set2' or 'Set3' or 'Spectral' or 'YlGn' or 'YlGnBu' or 'YlOrBr' or 'YlOrRd' or 'autumn' or 'binary' or 'black-white' or 'blue-red' or 'bone' or 'cool' or 'copper' or 'file' or 'flag' or 'gist_earth' or 'gist_gray' or 'gist_heat' or 'gist_ncar' or 'gist_rainbow' or 'gist_stern' or 'gist_yarg' or 'gray' or 'hot' or 'hsv' or 'jet' or 'pink' or 'prism' or 'spectral' or 'spring' or 'summer' or 'winter'

	- getSubList(self)
	  Purpose: get the attribute: subList

	- setSubList(self)
      Purpose: set the attribute: subList, by approx picking randomly a vertex per triangle to keep

	- plotScalTr(self, TrIndex, Index, Coords, Z, savename, sevSim=False, col = '', limScale = [None, None])
      Purpose: Does the plot of a scalar value on the mesh based on triangles
      Arguments: TrIndex: list of int: list of triangle indexes. ex: [t1, t2, t3,...]
                 Index: list of triplets of int: list containing the point indexes for each triangle. ex: [[p1, p2, p3],[p2, p3, p4], [p7, p8, p9],...]
                 Coords: list of 2 lists of triplets of floats: X, Y coordinates of the vertices with the deformation. ex: [[[x1, x2, x3],[x2, x3, x4] ,...], [[y1, y2, y3], [y2, y3, y4] , ...]]
                 Z: list of float: list of triangle values. ex: [z1, z2, z3,...]
                 savename: string: name of the saved file
                 sevSim: bool: whether the program is currently analyzing several simulations - or just one
                 col: string: kind of colorscale to use.
                 limScale: [float, float]: min and max of the color range. Default is [None, None] and in this case adapts to the data

	- plotDistr(self, Scalar, savename, currFile)
      Purpose: Plots a histogram of a scalar
	  Arguments: Scalar: list of float: list of all scalar values over the mesh, and saves mean and sd in currFile
                 savename: string: name of the saved file
                 currFile: string: name of the Log.txt file

	- getLength(self)
      Purpose: Returns the difference between highest and lowest point
      Output: float: difference between highest and lowest point

	- getWidth(self)
      Purpose: Returns the difference between the most on the left and the most on the right
      Output: float:difference between the most on the left and the most on the right

	- getxyminmax(self)
      Purpose: Returns xmin, xmax, ymin and ymax
      Outputs: float: xmin
               float: xmax
               float: ymin
               float: ymax

	- getPerim(self)
      Purpose: Returns the sum of the length of the edges in the boundary
      Output: float: sum of the length of the edges in the boundary of the mesh

	- getArea(self)
      Purpose: Returns the sum of the areas of all the triangles
      Output: float: sum of the areas of all the triangles

	- getDistrParam(self, string)
      Purpose: Returns mean and standard deviation of a scalar
      Argument: string: name of a scalar ('S' in Output_names.txt)
      Output: [float, float]: [mean, standard deviation] for this scalar

	- findVert(self, ind)
      Purpose: Returns a point object when its index is given
      Argument: ind: int: index of the point
      Output: point: corresponding point

  - TrArea(X, Y)
    Purpose: Returns the area of a triangle
    Arguments: X: triplet of float: x-coordinates of the triangle
               Y: triplet of float: y-coordinates of the triangle (in the same order as X)
    Output: float: area of the triangle

##################################

class_point.py

  - class point
    - __init__(self)
      Purpose: Initialization of the object
      Attributes:
      self.Index_Vertex = None # index of the vertex
      self.nbTriangle = None # index of the triangle in which the vertex is
      self.nbInTr = None # index of the vertex in the triangle (in {0,1,2})
      self.coords = [None, None] # coordinates of the vertex ([x,y])
      self.displ = [None,None] # displacement of the vertex computed in the simulation [ux, uy]
      self.Elast = None # elasticity/young's modulus at this position
      self.RotatedMatrix = [None, None, None, None, None, None] # values in the rotated matrix: [AR, BR, IR, GR, ER, FR]; real matrix: [AR, BR, GR; BR, ER, FR; GR, FR, IR]
      self.HookeMatrix = [None, None, None, None] # values in the Hooke's matrix: [A1, A2, B12, C3]; real matrix: [A1, B, 0; B, A2, 0; 0,0, C3]
      self.Deformation = [None, None, None, None] # defromation matrix: [dxux, dxuy, dyux, dyuy]

    - initiate(self, row)
      Purpose: Initialization of a new point from a row containing all the info
      Argument: row: list of string: row from the source file, separated by ' '

    - getIndex(self)
      Purpose: get the attibute: Index_Vertex (index of the vertex)

    - setIndex(self, i)
      Purpose: set the attibute: Index_Vertex (index of the vertex)
      Argument: i: int: index

    - getnbTr(self)
      Purpose: get the attibute: nbTriangle (index of the triangle)

    - setnbTr(self, i)
      Purpose: set the attibute: nbTriangle (index of the triangle)
      Argument: i: int: index of the triangle

    - getnbInTr(self)
      Purpose: get the attibute: nbInTr (index of the vertex in the triangle)

    - setnbInTr(self, i)
      Purpose: set the attibute: nbInTr (index of the vertex in the triangle)
      Argument: i: int in {0,1,2}: index of the vertex in the triangle

    - getCoords(self)
      Purpose: get the attibute: coords ([x,y] coordinates of the vertex)

    - setCoords (self, c)
      Purpose: set the attibute: coords ([x,y] coordinates of the vertex)
      Argument: c: [float, float]: [x,y] coordinates of the vertex

    - getDispl(self)
      Purpose: get the attibute: displ ([ux,uy] coordinates of the displacement at this vertex)

    - setDispl(self, u)
      Purpose: set the attibute: displ ([ux,uy] coordinates of the displacement at this vertex)
      Argument: u: [float, float]: [ux,uy] coordinates of the displacement at this vertex

    - getElast(self)
      Purpose: get the attibute: Elast (value of the elasticity at this vertex)

    - setElast(self, Y)
      Purpose: set the attibute: Elast (value of the elasticity at this vertex)
      Argument: Y: float: value of the elasticity at this vertex

    - getRotMat(self):
      Purpose: get the attibute: RotatedMatrix (rotated elasticity matrix at this vertex)

    - setRotMat(self, mat)
      Purpose: set the attibute: RotatedMatrix (rotated elasticity matrix at this vertex)
      Argument: mat: list of 6 float: rotated elasticity matrix at this vertex

    - getHookeMat(self)
      Purpose: get the attibute: HookeMatrix (elasticity matrix at this vertex)

    - setHookeMat(self, mat)
      Purpose: set the attibute: HookeMatrix (elasticity matrix at this vertex)
      Argument: mat: list of 4 float: elasticity matrix at this vertex

    - getDeformation(self)
      Purpose: get the attibute: Deformation (deformation matrix at this vertex)

    - setDeformation(self, a)
      Purpose: set the attibute: Deformation (deformation matrix at this vertex)
      Argument: a: list of 4 float: deformation matrix at this vertex

    - getGrArea(self)
      Purpose: Returns a vertex-based value of growth (rate) from the deformation: 1+dux/dx+duy/dy
      Output: float: growth (rate) at this vertex

    - getOne(self)
      Purpose: Returns 1 for some plots (only contour for example)
      Output: int: 1

  - Diag(Matrix)
    Purpose: Returns the eigen values of a 2x2 matrix
    Argument: list of 4 float [a,b,c,d]                : matrix [a,b;c,d]
           OR list of 3 float [a,b,c]                  : matrix [a,c;c,b] (symetric matrix)
           OR list of 2 lists of 2 floats [[a,b],[c,d]]: matrix [a,b;c,d]
    Output: list of 2 float: [first eigenvalue, second eigenvalue]

  - Sym(Matrix)
    Purpose: Returns a symmetric 2x2 matrix, where the not on the diagonal terms are meaned (see the deformation formalism)
    Argument: list of  4 float [a,b,c,d]: matrix [a,b;c,d]
    Output: list of  4 float [a,b,b,d]: matrix [a,b;b,d]

  - tofloat(liste)
    Purpose: Converts to floats strings in a list
    Argument: list of strings
    Output: list of floats corresponding to the strings
    Note: there is no check here because this is used to read a file automatically generated, with only numbers
