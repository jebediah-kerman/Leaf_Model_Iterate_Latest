#! /usr/bin/python
import mahotas
import numpy as np
from scipy import ndimage as nd
from scipy.misc import imsave

import mahotas.polygon
import matplotlib.pyplot as plt
import pylab
import os
import sys
from sys import argv
import csv #MATHILDE
#import os.path #MATHILDE

def inverse(I):
	inv = 255*np.ones_like(I)-I
	return inv

def define_structure_element2(d) :
	"""
	defines a spherical structuring element of radius d
	"""
	struct = nd.generate_binary_structure(2, 1)
	struct.astype(int)
	struct = nd.iterate_structure(struct, d).astype(int)
	return struct


def define_structure_element(d):
	s= 2*d+1
	ind_values = range(s)
	imsize = (s,s)
	initial=np.zeros(imsize)
	threshold=1
	c=(d,d)
	for x in ind_values:
		for y in ind_values:
			norm2=(x-c[0])**2+(y-c[1])**2
			if norm2<=d**2 :
				initial[x,y]=threshold
	return initial


def write_vector2vx(y,filename) :
	FILE = open(filename,"w")
	FILE.write('# 6\n')
	FILE.write('# 2\n')
	s = y.shape
	FILE.write('# '+str(s[0])+'\n')
	for i in range(0,s[0]) :
		FILE.write(str(y[i,0])+' '+str(y[i,1])+'\n')
	FILE.close()
	return


def read_contour(filename) :
	"""
	Reads the coordinates of the point-list defining a 2D curve.

	**Parameters**

	filename : string
		name of the file containing the coordinates of the points 
		defining a curve. The file has to contain the 2D cartesian coordinates
		of the points in two columns, one column for each dimension.

	**Returns**

	y : 2-column array
		defines a curve as an ordered list of points, 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates

	"""
	x0 = []
	for line in file(filename):
		if line[0] == "#" :
			line = line.lstrip('# ')
			line = line.rstrip('\n')
		else :
			line = line.rstrip('\n')
			line_list = [float(x) for x in line.split(' ')]
			x0.append(line_list)
	y = np.array(x0)
	return y


def fig_contour(vect,titre) :
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.axis('equal')
	ax.plot(vect[:,0], vect[:,1],'bo')
	ax.plot(vect[:,0], vect[:,1],'b-')
	ax.plot(vect[0,0],vect[0,1],'ro')
	ax.set_title(titre)
	plt.grid(True)
	plt.show()
	return


def resample(y,N_sample):
	# k is the degree of the spline
	# s is a smoothing factor and we do not want to smooth at this stage
	tck,u=interpolate.splprep([y[:,0],y[:,1]],k=2, s=0)
	# evaluate the fitted curve on N_sample points
	u_sample = np.linspace(u.min(), u.max(), N_sample)
	# the resampled
	y_resampled=(np.array(interpolate.splev(u_sample,tck))).transpose()
	return y_resampled


def surface2vector(surf, dist):
	"""
	from a binary image creates a list of points on a grid with step "dist"
	"""
	s = surf.shape
	print s
	y = []
	for i in range(0,s[0],dist):
		for j in range (0,s[1],dist):
			if surf[i,j]>0 :
				y.append([i,j])
	return np.array(y)


def princomp(A):
	""" performs principal components analysis 
	 (PCA) on the n-by-p data matrix A
	 Rows of A correspond to observations, columns to variables. 

	Returns :  
	coeff :
	is a p-by-p matrix, each column containing coefficients 
	for one principal component.
	score : 
	the principal component scores; that is, the representation 
	of A in the principal component space. Rows of SCORE 
	correspond to observations, columns to components.
	latent : 
	a vector containing the eigenvalues 
	of the covariance matrix of A.
	"""
	# computing eigenvalues and eigenvectors of covariance matrix
	M = (A-np.mean(A.T,axis=1)).T # subtract the mean (along columns)
	[latent,coeff] = np.linalg.eig(np.cov(M)) # attention:not always sorted
	score = np.dot(coeff.T,M) # projection of the data in the new space
	return coeff,score,latent
 
# Inputs and parameters
#-----------------------

filename = "1.JPG"
filename = argv[1]
fileroot = filename.split('.')[0]

Npoints = 400 # number of points on a half resampled leaf
Npoints = int(argv[2])

ScaleInfo = argv[4]

#MATHILDE
#scale
print('reading scale...')
scale = 370.0 #so if there is no info, it still works
#reading the file
if os.path.isfile('../../Scales.txt'):
	with open('../../Scales.txt', 'rb') as csvfile:
		scalereader = csv.reader(csvfile, delimiter = ' ', quotechar='|')
		for row in scalereader:
			#finding the right line
			if row[0] in ScaleInfo:
				scale = float(row[1])
#ENDMATHILDE

dist=int(scale/75)+1 # resampling distance in contourdetection

# Output directories
#--------------------

imagedir = 'Images/'
os.system("mkdir -p "+imagedir)

contourdir = 'Contours'+str(Npoints)+'/'
os.system("mkdir -p "+contourdir)


#######################################
# 1. read the image and detect the sepal
#######################################

# read the image
I = mahotas.imread(filename)


print '---------------------------------'
print 'Image ',filename,'read'
s = I.shape
# take the first channel
im = I[:,:,1]
pylab.imshow(im,cmap='gray')
imsave(imagedir+fileroot+'_output00.png',im)

s_im = im.shape

# smooth the image by a gaussian filter
fI = nd.filters.gaussian_filter(im,2*dist)
#pylab.imshow(fI,cmap='gray')
#pylab.show()

# compute threshold according to the Otsu method
T_otsu = mahotas.thresholding.otsu(im)
#print('T_otsu['+str(0)+']='+str(T_otsu))
# compute threshold according to the Riddler-Calvard method
T_rc = mahotas.thresholding.rc(im)
#print('T_rc['+str(0)+']='+str(T_rc))

# create binry image by taking the Otsu threshold
thresh = (fI> T_otsu)
pylab.imshow(thresh,cmap='gray')
imsave(imagedir+fileroot+'_output01.png',thresh)

# fill the holes
thresh = nd.morphology.binary_fill_holes(thresh)
pylab.imshow(thresh,cmap='gray')
imsave(imagedir+fileroot+'_output02.png',thresh)

# fill holes on the surface
struct = define_structure_element(2*dist)
thresh = nd.morphology.binary_closing(thresh, structure = struct)
pylab.imshow(thresh,cmap='gray')
imsave(imagedir+fileroot+'_output03.png',thresh)

# remove trychomes and background irregularities
thresh = nd.morphology.binary_opening(thresh, structure = struct)
pylab.imshow(thresh,cmap='gray')
imsave(imagedir+fileroot+'_output04.png',thresh)

print '... object detected'
#######################################
# 2. detect the contour
#######################################
# - as a list of ordered points
# - align along the Ox axis

# resample the sepal surface into vectors
y = surface2vector(thresh,int(dist/5)+1)
#y = surface2vector(thresh,dist)
print(len(y))
# change image (matrix) indices to XOy coordinates

y1 = np.zeros_like(y)
y1[:,0] = y[:,1]
y1[:,1] = -y[:,0]
y = y1
'''
y1 = np.zeros_like(y)
y1[:,0] = y[:,1]
y1[:,1] = y[:,0]
y = y1
'''

cm = np.array([np.mean(y[:,0]),np.mean(y[:,1])])
y[0,:]=cm

# detect its contour as a non-convex hull
contourname = contourdir+fileroot+'_contour.vx'
write_vector2vx(y,contourname)
os.system('../../SepalContour/ahull2py.R '+contourname+' '+str(2*dist))
y = read_contour(contourname)
print(len(y))
'''
################################
fig = plt.figure()
ax = fig.add_subplot(111)
ax.axis('equal')
ax.plot(y[:,0], y[:,1],'bo')
ax.plot(y[:,0], y[:,1],'b-')
ax.plot(y[0,0],y[0,1],'ro')
plt.grid(True)
fig.savefig(imagedir+fileroot+'_test2.png')	
##################################
'''


#fig_contour(y,'Original contour')
print '... contour detected'

from contours_home import *



# determine the shape principal components
coeff, score, latent = princomp(y)
m = np.mean(y.T,axis=1)


def rotate_on_Oy(z,n):
	"""
	rotates z so that unitary vector n becomes j vector on Oy
	"""
	R=np.array([-n[1],n[0],n[0],n[1]]).reshape((2,2))
	Rz=np.array([np.dot(R,y) for y in z])
	return Rz


# align the first principal component along the Ox axis
def place_contour(yy,m,latent,coeff):
	# Origo in the center of mass
	y = yy-m
	# the principal axis
	if latent[0]>=latent[1] : 
		n = coeff[:,0]
	else :
		n = coeff[:,1]
	# direction such that the basis of the sepal is in the positive Ox direction
	# -- we suppose that on the picture a) the sepal is vertical b) the base is on the bottom
	if n[1]<0 :
		n = -n
	# let us turn the principal axis on Oy
	y=rotate_on_Oy(y,n)
	return y

def get_main_direction(yy,m,latent,coeff):
	# Origo in the center of mass
	y = yy-m
	# the principal axis
	if latent[0]>=latent[1] : 
		n = coeff[:,0]
		lambd = latent[0]
	else :
		n = coeff[:,1]
		lambd = latent[1]
	print "PC1", lambd, "  ", n
	print latent, coeff
	# direction such that the basis of the sepal is in the positive Ox direction
	# -- we suppose that on the picture a) the sepal is vertical b) the base is on the bottom
	if n[1]<0 :
		n = -n
	return n


n = get_main_direction(y,m,latent,coeff)
n = n/np.sqrt(n[0]**2+n[1]**2)

#nnorm = np.array([-n[1],n[0]])


def fig_contour_axis(vect,titre,n,cm) :
	axis_vect = []
	axis_vect.append(cm)
	axis_vect.append(cm+500*n)
	axis_vect = np.array(axis_vect)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.axis('equal')
	ax.plot(vect[:,0], vect[:,1],'bo')
	ax.plot(vect[:,0], vect[:,1],'b-')
	ax.plot(vect[0,0],vect[0,1],'ro')
	ax.plot(axis_vect[:,0],axis_vect[:,1],'r-')
	ax.set_title(titre)
	plt.grid(True)
	fig.savefig(imagedir+fileroot+'_output05.png')
	#plt.show()
	return


#y = place_contour(y,m,latent,coeff, orientation)
#fig_contour(y,'Placed contour - using contour points')
fig_contour_axis(y,'Contour with main axis',n,cm)

y = place_contour(y,m,latent,coeff)
#######################################
# 3. reinit the contour so that the beginning is at the basis
#######################################

# detect the point with negative y coordinate, where the main axis sections the curve

min_dist_Oy = abs(y[0,1])+1
point_index = 0
for i in range(0,len(y)) :
	if (y[i,1]<0) and (abs(y[i,0])<min_dist_Oy) :
		min_dist_Oy = abs(y[i,0])
		point_index = i


def reinit_contour(y,index):
	"""
	reinits contour so that the first point will be that of the actual "ind" point
	"""
	z=np.zeros_like(y)
	z[:len(y)-index]=y[index:]
	print "reinit_contour : index=",index
	z[len(y)-index:]=y[:index]
	return z


y1 = reinit_contour(y,point_index)
#fig_contour(y1, 'toto')

# Place the origin in the base point of the sepal 
y1 = y1 - y1[0]


# resample the contour so that we have 2*Npoints equidistant points on it
y = resample(y1,2*Npoints)

def fig_contour_placed(vect,titre) :
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.axis('equal')
	ax.plot(vect[:,0], vect[:,1],'bo')
	ax.plot(vect[:,0], vect[:,1],'b-')
	ax.plot(vect[0,0],vect[0,1],'ro')
	ax.set_title(titre)
	plt.grid(True)
	fig.savefig(imagedir+fileroot+'_output06.png')
	#plt.show()
	return

y=y/scale

fig_contour_placed(y,'Sepal placed on Ox axis')


#######################################
# 4. measuring length, width, area...
#######################################
print '... measuring length, width, area...'

def compute_perimeter(y):
    p = 0
    for i in range(0,len(y)-1):
        p+=norm(y[i+1]-y[i])
    return p

base,pointe,length=compute_leaf_axis(y)
pointe_ind = closest(y,pointe)


length = max(y[:,1])-min(y[:,1])
width = max(y[:,0])-min(y[:,0])

#length = max(proj1)-min(proj1)
#width = max(proj2)-min(proj2)

area = sum(sum(thresh))/scale**2
perimeter = compute_perimeter(y)
circularity = 4*np.pi*area/perimeter**2
aspect_ratio = width/length

# CAUTION : measurements are in pixels....

FILE1 = open('measurements.txt',"a")
FILE1.write(fileroot+';')
FILE1.write(str(length)+';')
FILE1.write(str(width)+';')
FILE1.write(str(area)+';')
FILE1.write(str(perimeter)+';')
FILE1.write(str(aspect_ratio)+';')
FILE1.write(str(circularity))
FILE1.write('\n')
FILE1.close()


#######################################
# 5. write the contour to R
#######################################

print '... writing the contour of length ',len(y),'...'
write_vector2vx(y,contourname)

datadir = 'ContoursData_sample'+str(Npoints)

os.system("mkdir -p "+datadir)

lm_number = 0
ind_landmarks = [pointe_ind]

'''
print fileroot+' '
print str(length)+' '
print str(lm_number)+' '
print str(lm_number)+' '
print ind_landmarks
print y
'''

FILE = open(datadir+'/X_contours'+str(lm_number)+'.data',"a")
FILE.write(fileroot+' ') # 0
FILE.write(str(length)+' ') # 1
FILE.write(str(lm_number)+' ') # 2
for i in ind_landmarks :
	FILE.write(str(i)+' ')
for i in y :
	FILE.write(str(i[0])+' ') # x coordinates
FILE.write('\n')
FILE.close()
#
FILE = open(datadir+'/Y_contours'+str(lm_number)+'.data',"a")
FILE.write(fileroot+' ') # 0
FILE.write(str(length)+' ') # 1
FILE.write(str(lm_number)+' ') # 2
for i in ind_landmarks :
	FILE.write(str(i)+' ')
for i in y :
	FILE.write(str(i[1])+' ') # y coordinates
FILE.write('\n')
FILE.close()

