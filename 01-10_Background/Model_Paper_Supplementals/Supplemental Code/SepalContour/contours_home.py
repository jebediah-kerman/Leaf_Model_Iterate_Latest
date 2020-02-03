import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft
from scipy.linalg import norm
from scipy.signal import convolve
from scipy import interpolate

from contours import *
#import pywt

def extract_tooth(yd,ind_landmarks_d) :
	tooth_d = yd[(ind_landmarks_d[0]):(ind_landmarks_d[1])]
	# tooth_base
	tb_d = yd[ind_landmarks_d[1]]-yd[ind_landmarks_d[0]]
	toothlength_d = norm(tb_d)
	# put the origo in landmrk[1]
	tooth_d=tooth_d-tooth_d[0]
	# rotate so that the axis of the leaf is on Ox
	# ( the + direction is oriented towards lm0 )
	n=tooth_d[-1]/toothlength_d
	tooth_d=rotate(tooth_d,n)/toothlength_d
	return tooth_d, toothlength_d

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

def superpose_fig_contour(vect1,vect2,titre) :
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.axis('equal')
	ax.plot(vect1[:,0], vect1[:,1],'b')
	ax.plot(vect1[0,0],vect1[0,1],'bo')
	ax.plot(vect2[:,0], vect2[:,1],'r')
	ax.plot(vect2[0,0],vect2[0,1],'ro')
	ax.set_title(titre)
	plt.grid(True)
	plt.show()
	return

def compute_leaf_axis(z):
	"""
	argument : z=(N,2) array
	returns the coordiantes of the petiole-middle and of the leaftip
	"""
	start=z[0]
	end=z[-1]
	petiole_middle=(start+end)/2.0
	a=z-petiole_middle
	dist=[norm(a[i]) for i in xrange(0,len(a))]
	index=dist.index(max(dist))
	return petiole_middle,z[index],max(dist)

def rotate(z,n):
	"""
	rotates z so that unitary vector n becomes i vector on Ox
	"""
	R=np.array([n[0],n[1],-n[1],n[0]]).reshape((2,2))
	Rz=np.array([np.dot(R,y) for y in z])
	return Rz


def split_leaf_list(pointe_index,point_index_list):
	"""
	Splits the pointlist into two sublists delimited by its point 
	with index `pointe_index`. Both subcurves have their starting point at
	the delimiting point. The second subcurve (left) remains unchainged. The
	first one (right) is mirrored with respect to the axis defined by the delimiting point
	and the mean of the two extremity points of the original curve.

	**Parameters**

	pointe_index : positive integer
		index of point which defines the split of the curve

	point_index_list : list of indices (positive integers)
		indices on the y curve, to split

	**Returns**

	point_index_list_right : index_list
		defines the right list. Its indexation begins at the splitting point
		and it is reflected with respect to the axis defined by the delimiting point
		and the mean of the two extremity points of the original curve.

	point_index_list_left : index_list
		defines the left list of indices; 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates. Its indexation begins at the splitting point
		and its position is unchanged.

	**See also**

	glue_leaf()


	"""
	# the first half :
	list_right = [i for i in point_index_list if i<pointe_index]
	list_right = [pointe_index-i for i in list_right]
	list_right.reverse()
	# the second half - the indices are shifted
	list_left = [i for i in point_index_list if i>pointe_index]
	list_left = [i-pointe_index for i in list_left]
	return list_right,list_left


def find_tip_on_lobe_MaxMedian(lobeindexes,y):
	"""
	given a curve y on which lobeindexes define a tooth
	(!!! lobeindexes are defined on the curve y !!!)
	 - finds the coordinates of the point wich is the more distant
	 from the middle of the tooth basis
	"""
	tooth=y[lobeindexes[0]:lobeindexes[1]+1]
	# the middle of the tooth-basis :
	mid=(tooth[0]+tooth[-1])/2.
	# distances from the middle :
	return farest(tooth,mid)+lobeindexes[0]


def find_tip_on_lobe_farest(lobeindexes,y):
	"""
	given a curve y on which lobeindexes define a tooth
	(!!! lobeindexes are defined on the curve y !!!)
	 - finds the coordinates of the point wich is the more distant
	 from the tooth basis
	"""
	tooth=y[lobeindexes[0]:lobeindexes[1]+1]
	# the farest point from the basisline
	farest_dist,farest_index=find_most_distant_point(tooth)
	return farest_index+lobeindexes[0]


def find_tips(yy,ind_kmin) :
	# split the contour into lobes
	lobeindexes = find_lobes(ind_kmin)
	# on each lobe identify a tip as a point which is the farest from the middle of the base
	ind_kmax = [find_tip_on_lobe_farest(lobeindexes[i],yy) for i in range(0,len(lobeindexes)) ]
	return ind_kmax


def allDescendants(G,d1):
   d2 = []
   for d in d1:
       d2 += G.successors(d)
   return d2


def daughter_nodes(Gg,i):
	node_list_old = [i]
	descend = Gg.successors(i)
	node_list = node_list_old + descend
	while len(node_list)!=len(node_list_old) :
		node_list_old = node_list
		descend = allDescendants(Gg,descend)
		node_list = node_list + descend
	return node_list


def FirstPrimaryLobe(Gg,lg,zg,ind_kmin_g):
	"""
	finds the lobe corresponding to the whole 1st primary tooth
	"""
	# the primary teeth are all successors of [0] :
	Ptg=Gg.successors(0)
	Primary_lobes=[]
	for i in Ptg :
		# find all the daughter teeth
		# (together they construct the whole primary tooth)
		#daught_graphe=nx.dfs_tree(Gg,i)
		#daught_nodes=daught_graphe.nodes()
		daught_nodes = daughter_nodes(Gg,i)
		# list of teeth-extremities
		extr_list=[]
		for j in daught_nodes :
			extr_list=extr_list+lg[j].indexes
		Primary_lobes.append(leaf_lobe(1,[min(extr_list),max(extr_list)],zg,ind_kmin_g))
	# basewidths of the primary teeth :
	basewlist=[lobe.basewidth() for lobe in Primary_lobes]
	# find the index of the largest tooth :
	ind=basewlist.index(max(basewlist))
	# so the lobe corresponding to the largest primary tooth is
	P1lobe=Primary_lobes[ind]
	# and its label in the graphe
	P1label=Ptg[ind]
	return P1label, P1lobe


def FindPrimaryLobes(Gg,lg,zg,ind_kmin_g):
	"""
	finds the lobe corresponding to the whole primary teeth
	- omitting distal teeth
	"""
	# the primary teeth are all successors of [0] :
	Ptg=Gg.successors(0)
	Primary_lobes=[]
	for i in Ptg :
		# find all the daughter teeth
		# (together they construct the whole primary tooth)
		#daught_graphe=nx.dfs_tree(Gg,i)
		#daught_nodes=daught_graphe.nodes()
		daught_nodes = daughter_nodes(Gg,i)
		# list of teeth-extremities
		extr_list=[]
		for j in daught_nodes :
			extr_list=extr_list+lg[j].indexes
		Primary_lobes.append(leaf_lobe(1,[min(extr_list),max(extr_list)],zg,ind_kmin_g))
	# sort the primary teeth with respect to their position on the contour
	# first the one which is closest to the tip (might be a distal tooth)
	sort_key = lambda lobe: lobe.indexes[0]
	Primary_lobes.sort(key=sort_key)
	# compute basewidths of the primary teeth :
	basewlist=[lobe.basewidth() for lobe in Primary_lobes]
	# find the index of the largest tooth (this being the first primary tooth)
	ind=basewlist.index(max(basewlist))
	# collect primary teeth which come after the first primary tooth :
	PrimaryTeeth = Primary_lobes[ind:]
	# and their label in the graphe
	#PrimaryLabels=Ptg[ind:]
	return  PrimaryTeeth


def FindLobeLandmarks(PrimaryTeeth):
	"""
	finds the indexes of the points delimitating the teeth in 'PrimaryTeeth'
	"""
	landmarks = [lobe.indexes[0] for lobe in PrimaryTeeth]
	landmarks.append(PrimaryTeeth[-1].indexes[1])
	return landmarks

def check_singular_curve(y):
	"""
	if double points appear in a curve, then the second pops out
	"""
	y_list=list(y)
	pop_indices = []
	for i in range(0,len(y)-1):
		if y[i,0]==y[i+1,0] and y[i,1]==y[i+1,1]:
			pop_indices.append(i)
	pop_indices.reverse()
	for i in pop_indices:
		y_list.pop(i)
	return np.array(y_list)

def resample(y,N_sample):
	# k is the degree of the spline
	# s is a smoothing factor and we do not want to smooth at this stage
	tck,u=interpolate.splprep([y[:,0],y[:,1]],k=2, s=0)
	# evaluate the fitted curve on N_sample points
	u_sample = np.linspace(u.min(), u.max(), N_sample)
	# the resampled
	y_resampled=(np.array(interpolate.splev(u_sample,tck))).transpose()
	return y_resampled


def resample_smooth(y,N_sample,smooth):
	# k is the degree of the spline
	# s is a smoothing factor
	tck,u=interpolate.splprep([y[:,0],y[:,1]],k=2, s=smooth, nest=-1)
	# evaluate the fitted curve on N_sample points
	u_sample = np.linspace(u.min(), u.max(), N_sample)
	# the resampled
	y_resampled=(np.array(interpolate.splev(u_sample,tck))).transpose()
	return y_resampled
