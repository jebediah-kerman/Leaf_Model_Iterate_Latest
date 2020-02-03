import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import norm, solve, det
from scipy import argmin, argmax
from math import pi, acos
#import networkx as nx


def farest(y,p):
	"""
	Given a curve *y* as a pointlist, the function computes the index of the
	point on the curve which is the farest to a point *p*.

	**Parameters**

	y : 2D curve

	p : (1,2) dimensional array
		the cartesian coordinates of a point

	**Returns**

	index : positive integer
		the index of the farest point

	"""
	a=y-p
	dist=[norm(a[i]) for i in xrange(0,len(a))]
	index=dist.index(max(dist))
	return index


def closest(y,p):
	"""
	Given a curve *y* as a pointlist, the function computes the index of the
	point on the curve which is the closest to a point *p*.

	**Parameters**

	y : 2D curve

	p : (1,2) dimensional array
		the cartesian coordinates of a point

	**Returns**

	index : positive integer
		the index of the closest point

	"""
	a=y-p
	dist=[norm(a[i]) for i in xrange(0,len(a))]
	minimal_distance=min(dist)
	index=dist.index(minimal_distance)
	return index


def closest_list(y,p_list):
	"""
	Given a curve *y* as a pointlist, the function computes the 
	list of indexes of the
	points on the curve which are the closest to the points in  *p_list*.

	**Parameters**

	y : 2D curve

	p_list : (n,2) dimensional array
		the cartesian coordinates of *n* points

	**Returns**

	index_list : list of positive integers
		the list of indexes of the closest points on the curve

	"""
	index_list=[]
	for i in range(0,len(p_list)):
		index_list.append(closest(y,p_list[i]))
	return index_list


def point_on_the_curve(p,y):
	"""
	Given a curve *y* as a pointlist, checks if the point *p* is on the curve.

	**Parameters**

	p : (1,2) dimensional array
		the cartesian coordinates of a point

	y : 2D curve

	**Returns**

	on_the_curve : boolean
		*True* if the point is in the list of points defining the curve

	"""
	on_the_curve = 0
	i = 0
	while on_the_curve == 0 and i<len(y) :
		if y[i,0]==p[0] and y[i,1]==p[1] :
			on_the_curve = 1
		i+=1
	return on_the_curve


def compute_angle(m,o,n):
	"""
	Given three points m, o, n in 2D space, the function computes the angle
	between the segments om and on, in the direction from m towards n.

	**Parameters**

	m, o, n : 2 dimensional vectors
		2D cartesian coordinates of the three points

	**Returns**

	alpha : real number in the range [0;2pi[
		the angle defined by the points m, o and n

	**Examples**

	>>> import numpy as np
	>>> from contours import *
	>>> a = np.array([1,0])
	>>> b = np.array([0,0])
	>>> c = np.array([0,2])
	>>> compute_angle(a,b,c) == pi/2
	True
	"""
	v1 = m-o
	v2 = n-o
	cosalpha = np.dot(v1,v2)/(norm(v1)*norm(v2))
	# ... numerical approximations might cause that cosalpha is out of the domain [-1,1]
	# Therefore, we will force the right domain :
	if cosalpha<-1 :
		cosalpha = -1
	if cosalpha>1 :
		cosalpha = 1
	# acos gives angle values in the range [0;pi].
	# we need to compute angle values in the whole range [0;2pi[
	if np.cross(v1,v2) >= 0 :
		alpha = acos(cosalpha)
	else :
		alpha = 2*pi-acos(cosalpha)
	return alpha

def evaluate(y,indices):
	"""
	**Parameters**

	y : 2-column array
		defines a curve as an ordered list of points, 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates

	indices : list of positive integers
		list of indices

	**Returns**

	y_sub : 2-column array
		each row of the array contains the cartesian coordinates of the points
		with indices in the list `indices` on the curve `y`.

	"""
	values = []
	for i in indices :
		values.append(y[i])
	return np.array(values)

def split_leaf(y,pointe_index):
	"""
	Splits the 2D curve $Fz$ into two subcurves delimited by its point 
	with index `pointe_index`. Both subcurves have their starting point at
	the delimiting point. The second subcurve (left) remains unchainged. The
	first one (right) is mirrored with respect to the axis defined by the delimiting point
	and the mean of the two extremity points of the original curve.

	**Parameters**

	y : 2-column array
		defines a curve as an ordered list of points, 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates

	pointe_index : positive integer
		index of point which defines the split of the curve

	**Returns**

	y_right : 2-column array
		defines the right sub-curve of the split as an ordered list of points; 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates. Its indexation begins at the splitting point
		and it is reflected with respect to the axis defined by the delimiting point
		and the mean of the two extremity points of the original curve.

	y_left : 2-column array
		defines the left sub-curve of the split as an ordered list of points; 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates. Its indexation begins at the splitting point
		and its position is unchanged.

	**See also**

	glue_leaf()


	"""
	# define the mirroring axis of the curve
	start = y[0]
	end = y[-1]
	petiole_middle = (start+end)/2.0
	axis_p1 = y[pointe_index]
	axis_p2 = petiole_middle
	# the first half is reversed and mirrored on the axis
	y_right = list(y[:pointe_index+1])
	y_right.reverse()
	y_right = [symmetric_point(point,axis_p1,axis_p2) for point in y_right]
	y_right = np.array(y_right)
	# the second half is unchaged
	y_left = (y[pointe_index:]).copy()
	return y_right,y_left

def glue_leaf(index_d,index_g,pointe_index):
	"""
	Given the lists of point indexes `index_d` and `index_g` on the
	right and left half of an original curve splitted by the point
	with index `pointe_index`, this function
	gives the list of indeces on the entire contour.

	**Parameters**

	index_d : list of positive integers
		list of point-indexes on the right subcurve

	index_g : list of positive integers
		list of point-indexes on the left subcurve

	pointe_index : positive integer
		index of the point on the entire curve, which defined the split
		into the right and left subcurves.

	**Returns**

	list_index : list of positive integer
		list of point-indexes on the whole curve

	**See also**

	split_leaf()

	"""
	ii_d = list(index_d)
	ii_g = list(index_g)
	# the first half is reversed
	ii_d.reverse()
	# reconstitute the indexation on the whole curve
	ii_d = [pointe_index-ind for ind in ii_d]
	ii_g = [pointe_index+ind for ind in ii_g]
	return ii_d+ii_g

def symmetric_point(x,m,n):
	"""
	Finds the symmetric of the point `x` with 
	respect to the line defined by the points `m` and `n`

	**Parameters**

	x, m, n : (1,2) dimensional arrays
		cartesian coordinates of the three points

	**Returns**

	xs : (1,2) dimensional array
		cartesian coordinates of the reflection of `x` with respect to the
		axis defined by the points `m` and `n`

	"""
	e = (n-m)/norm(n-m) # vecteur directeur
	y = x-m
	yt = np.dot(y,e)*e
	xs = 2*yt-y+m
	return xs

def compute_leaf_point(z):
	"""
	Given a leafcontour `z` in 2D the function
	determines the leaftip as the most distant point on the leafcontour from the middle of the
	leaf basis. The input point-list defining the leafcontour is supposed to be
	constructed in such a way that its first and last points
	mark the basis of the leaf on the right respectively on the left side of the leaf.

	**Parameters**

	z : 2-column array
		defines a curve as an ordered list of points, 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates

	**Returns**

	index : positive integer
		the index of the point corresponding to the leaftip

	leaflength : positive real number
		the maximal distance between the middle of the leafbasis 
		and the most distant point from it

	"""
	start = z[0]
	end = z[-1]
	petiole_middle = (start+end)/2.0
	a = z-petiole_middle
	dist = [norm(a[i]) for i in xrange(0,len(a))]
	index = dist.index(max(dist))
	return index,max(dist)

def intersect_segments(u,v):
	"""
	Checks if two segments intersect each other.

	**Parameters**

	u, v : (2,2) dimensional arrays
		One (2,2) dimensional array represents a segment defined by its 
		two extremity points. The first row contains the cartesian coordinates of the
		first point, the second row contains the cartesian coordinates of the second point.

	**Returns**

	val : boolean
		`True` if the two segments intersect each-other, `False` otherwise.

	**Examples**

	>>> import numpy as np
	>>> from contours import intersect_segments
	>>> a = np.array([[0,0],[0,1]])
	>>> b = np.array([[-1,0.5],[1,0.5]])
	>>> intersect_segments(a,b)
	True
	>>> c = np.array([[0,0],[1,0]])
	>>> intersect_segments(c,b)
	False
	>>> d = np.array([[0,0],[0,0.5]])
	>>> intersect_segments(d,b)
	True

	"""
	val = False
	du = u[1]-u[0]
	dv = v[1]-v[0]
	A = (np.vstack((-du,dv))).transpose()
	if det(A) != 0 :
		b = u[0]-v[0]
		l = solve(A,b)
		if (l>= 0).sum()==2 and (l<=1).sum()==2 :
			val = True
	return val

def detect_loops(z) :
	"""
	Given a 2D curve defined by a list of points, the function
	gives a list of point-pairs on the curve which define a loop. 
	The point-pairs are designated by their indexes in the point-list.

	**Parameters**

	z : 2-column array
		defines a curve as an ordered list of points, 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates

	**Returns**

	list_loops : list of index-pairs, where the index-pairs are of the form of a list of two positive integers

	"""
	list_loops = []
	i = 0
	while i<len(z)-1 :
		u = z[i:i+2]
		k = i+2
		while k<len(z)-1 :
			v = z[k:k+2]
			if intersect_segments(u,v) :
				list_loops.append([i,k+1])
				i = k+1
			k += 1
		i += 1
	return list_loops

def find_lobes(ind_kmin):
	"""
	Given a list of indices returns a list of pairs of consecutive indices.

	**Parameters**

	ind_kmin : list of positive integers
		list of indexes

	**Returns**

	lobeindexes : list of pairs of consecutive indices

	**Examples**

	>>> import numpy as np
	>>> from contours import find_lobes
	>>> a = [2, 6, 45, 90, 100]
	>>> find_lobes(a)
	[[2, 6], [6, 45], [45, 90], [90, 100]]

	"""
	lobeindexes = []
	for i in xrange(0,len(ind_kmin)-1):
		i1 = ind_kmin[i]
		i2 = ind_kmin[i+1]
		lobeindexes.append([i1,i2])
	return lobeindexes


def find_tip_on_lobe(sinus_indices,k):
	"""
	Given a vector of real values `k`, the function identifies 
	the element of maximal value on the subset of values with
	indices between the indices `sinus_indices`.
	 
	**Parameters**

	sinus_indices : list of two positive integers

	k : vector of real values

	**Returns**

	tip_index : positive integer comprised between the indices `sinus_indices`
		the index of the highest value

	"""
	# determine the maximum of k on the lobe
	#ind_max = argmax(k[sinus_indices[0]+1:sinus_indices[1]])
	ind_max = argmax(k[sinus_indices[0]:sinus_indices[1]])
	return sinus_indices[0]+ind_max

def compute_remarkable_points(y_smoothed):
	"""
	Computes the remarkable points (sinuses and toothtips) on the 2D contour *y_smoothed*,
	The input contour *y_smoothed* is supposed to be a segment approximation of a curve.
	(for example by using the Douglas-Peucker algorithme). Sinuses and toothtips are 
	determined as local extremas of the curvature on the curve and curvature is approximated 
	by angles between two consecutive segments.

	*Definitions*

	- sinus = negative curvature point on the approximated curve
	- tooth = a portion of the curve defined by two consectutive sinuses
	- toothtip = maximal curvature point on a tooth

	*The algorithme*

	- compute the curvature in all points defining the curve
	- determine the list of sinuses as local minimum points of negative curvature
	- detect loops on the curve
	- for each detected loop a sinus is added to the list; this sinus is determined as the maximum curvature on the loop
	- if one has a positive curvature point between the last detected sinus and the end of the curve we consider that we omitted the last tooth and we add an additional sinus. The added sinus is the minimal curvature point of that portion - even if it is not necessarily negative.
	- define a tooth by a pair of successive sinuses and detect all teeth on the curve
	- for each tooth identify the toothtip as the global maximum point of the curvature

	**Parameters**

	y_smoothed : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	**Returns**

	ind_kmin : list of positive integers
		list of indices of sinus points on the contour

	ind_kmax : list of positive integers
		list of indices of toothpoints points on the contour


	**Examples**

	>>> from contours import *
	>>> import numpy as np
	>>> curve=np.zeros((200,2))
	>>> curve[:,0] = np.linspace(0,4*pi,200)
	>>> curve[:,1] = np.sin(curve[:,0])
	>>> sinuses,tips = compute_remarkable_points(curve)
	Sinus added at the end
	There are  3  sinuses identified.
	>>> sinuses, tips # the indices of sinuses and tips respectiveley
	([25, 124, 199], [75, 174]) 
	>>> evaluate(curve,sinuses) # the coordinates of the three sinuses
	array([[  1.57868978e+00,   9.99968847e-01],
		[  7.83030129e+00,   9.99719634e-01],
		[  1.25663706e+01,  -4.89842542e-16]])

	"""
	# compute the curvature
	k = compute_curvature_by_angles(y_smoothed)
	# consider the negative curvature values on the curve
	k_minus = (k <= 0.0)*k
	# identify local minimas of the signal "k_minus"
	ind_kmin=find_local_minima(k_minus)
	'''
	# on each side of the leaf at least one sinus has to be identified
	# -----------------------------------------------------------------
	if sum(kmin) == 0 :
		sin2 = argmin(k)
		kmin[sin2] = 1
		print 'existence of at least one sinus forced'
	'''
	# take into account loops :
	# -------------------------
	Floops = detect_loops(y_smoothed)
	# if loops are present, the maximum curvature on a loop is in fact a sinus
	# --- we have to identify one and only one sinus on the loop ---
	for loop in Floops :
		# if there was a curvature minimum identified on the loop, we remove it
		for sin_index in ind_kmin :
			if sin_index in range(loop[0],loop[1]):
				ind_kmin.remove(sin_index)
		# and append the maximal curvature as sinus on the loop
		ind_sinus = argmax(k[loop[0]:loop[1]])
		ind_kmin.append(loop[0]+ind_sinus)
	# put the indexes of sinuses in ascendent order
	ind_kmin.sort()
	#
	# last lobe omitted ?
	# --------------------
	# if one has a positive curvature point between the last detected sinus and the end of the curve
	# we consider that we omitted the last tooth and we add an additional sinus.
	# The added sinus is the minimal curvature point of that portion - even if it is not necessarily negative.
	if len(ind_kmin)>0 :
		if ind_kmin[-1]<len(y_smoothed)-1 :
			i1 = ind_kmin[-1]
			i2 = len(y_smoothed)-1
			if max(k[i1:i2])>0:
				indmax = argmin(-k[i1:i2])
				sin1 = argmin(k[i1+indmax:i2])
				ind_kmin.append(i1+indmax+sin1+1)
				print "Sinus added at the end"
	#
	nb_kmin = len(ind_kmin)
	#print 'There are ',nb_kmin,' sinuses identified.'
	# split the contour into an ensemble of lobes, separated by the minimas
	lobeindexes = find_lobes(ind_kmin)
	#
	# on each lobe identify tooth tip as global maximum of the curvature
	# ------------------------------------------------------------------
	ind_kmax = []
	for i in xrange(0,len(lobeindexes)):
		ind = find_tip_on_lobe(lobeindexes[i],k)
		ind_kmax.append(ind)
	return ind_kmin,ind_kmax


def compute_remarkable_points_on_curve(y_smoothed):
	"""
	Computes the remarkable points (sinuses and toothtips) on the 2D contour *y_smoothed*,
	The input contour *y_smoothed* is supposed to be a segment approximation of a curve.
	(for example by using the Douglas-Peucker algorithme). Sinuses and toothtips are 
	determined as local extremas of the curvature on the curve and curvature is approximated 
	by angles between two consecutive segments.

	*Definitions*

	- sinus = negative curvature point on the approximated curve
	- tooth = a portion of the curve defined by two consectutive sinuses
	- toothtip = maximal curvature point on a tooth

	*The algorithme*

	- compute the curvature in all points defining the curve
	- determine the list of sinuses as local minimum points of negative curvature
	- detect loops on the curve
	- for each detected loop a sinus is added to the list; this sinus is determined as the maximum curvature on the loop
	- if one has a positive curvature point between the last detected sinus and the end of the curve we consider that we omitted the last tooth and we add an additional sinus. The added sinus is the minimal curvature point of that portion - even if it is not necessarily negative.
	- define a tooth by a pair of successive sinuses and detect all teeth on the curve
	- for each tooth identify the toothtip as the global maximum point of the curvature

	**Parameters**

	y_smoothed : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	**Returns**

	ind_kmin : list of positive integers
		list of indices of sinus points on the contour

	ind_kmax : list of positive integers
		list of indices of toothpoints points on the contour


	**Examples**

	>>> from contours import *
	>>> import numpy as np
	>>> curve=np.zeros((200,2))
	>>> curve[:,0] = np.linspace(0,4*pi,200)
	>>> curve[:,1] = np.sin(curve[:,0])
	>>> sinuses,tips = compute_remarkable_points(curve)
	Sinus added at the end
	There are  3  sinuses identified.
	>>> sinuses, tips # the indices of sinuses and tips respectiveley
	([25, 124, 199], [75, 174]) 
	>>> evaluate(curve,sinuses) # the coordinates of the three sinuses
	array([[  1.57868978e+00,   9.99968847e-01],
		[  7.83030129e+00,   9.99719634e-01],
		[  1.25663706e+01,  -4.89842542e-16]])

	"""
	# compute the curvature
	k = compute_curvature_by_angles(y_smoothed)
	# consider the negative curvature values on the curve
	k_minus = (k <= 0.0)*k
	# identify local minimas of the signal "k_minus"
	ind_kmin=find_local_minima(k_minus)
	'''
	# on each side of the leaf at least one sinus has to be identified
	# -----------------------------------------------------------------
	if sum(kmin) == 0 :
		sin2 = argmin(k)
		kmin[sin2] = 1
		print 'existence of at least one sinus forced'
	'''
	# take into account loops :
	# -------------------------
	Floops = detect_loops(y_smoothed)
	# if loops are present, the maximum curvature on a loop is in fact a sinus
	# --- we have to identify one and only one sinus on the loop ---
	for loop in Floops :
		# if there was a curvature minimum identified on the loop, we remove it
		for sin_index in ind_kmin :
			if sin_index in range(loop[0],loop[1]):
				ind_kmin.remove(sin_index)
		# and append the maximal curvature as sinus on the loop
		ind_sinus = argmax(k[loop[0]:loop[1]])
		ind_kmin.append(loop[0]+ind_sinus)
	# put the indexes of sinuses in ascendent order
	ind_kmin.sort()
	#
	'''
	# last lobe omitted ?
	# --------------------
	# if one has a positive curvature point between the last detected sinus and the end of the curve
	# we consider that we omitted the last tooth and we add an additional sinus.
	# The added sinus is the minimal curvature point of that portion - even if it is not necessarily negative.
	if len(ind_kmin)>0 :
		if ind_kmin[-1]<len(y_smoothed)-1 :
			i1 = ind_kmin[-1]
			i2 = len(y_smoothed)-1
			if max(k[i1:i2])>0:
				indmax = argmin(-k[i1:i2])
				sin1 = argmin(k[i1+indmax:i2])
				ind_kmin.append(i1+indmax+sin1+1)
				print "Sinus added at the end"
	#
	'''
	nb_kmin = len(ind_kmin)
	print 'There are ',nb_kmin,' sinuses identified.'
	# split the contour into an ensemble of lobes, separated by the minimas
	lobeindexes = find_lobes(ind_kmin)
	#
	# on each lobe identify tooth tip as global maximum of the curvature
	# ------------------------------------------------------------------
	ind_kmax = []
	for i in xrange(0,len(lobeindexes)):
		ind = find_tip_on_lobe(lobeindexes[i],k)
		ind_kmax.append(ind)
	return ind_kmin,ind_kmax


def compute_remarkable_points_on_leaf(y,epsilon):
	"""
	Computes the remarkable points (sinuses and toothtips) on the 2D leaf-contour *y*
	using a Douglas-Peucker segment approximation of this curve with stop criterion *epsilon*.

	*Definitions*

	- sinus = negative curvature point on the approximated curve
	- tooth = a portion of the curve defined by two consectutive sinuses
	- toothtip = maximal curvature point on a tooth

	*The algorithme*

	- consider the Douglas-Peucker approximation of the curve *y* with stop criterion *epsilon*
	- determine the leaf's tip
	- split the approximated curve into the two halves defined by the leaftip
	- on each half determine the remarkable points
	- reconstitute the lists of sinuses and tips for the entire leafcontour

	**Parameters**

	y : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	epsilon : a real number
		the stop crierion distance in the Douglas-Peucker approximation

	**Returns**

	y_smoothed : 2-column array
		the Douglas-Peucker approximation of the input curve using *epsilon* as stop criterion

	ind_kmin : list of positive integers
		list of indices in the *y_smoothed* curve of sinus points on the contour. 

	ind_kmax : list of positive integers
		list of indices in the *y_smoothed* curve of toothpoints points on the contour. 

	"""
	# consider the Douglas-Peucker approximation of the curve y with stop criterion epsilon
	#y_smoothed = segment_approximation(y,epsilon)
	y_smoothed = DouglasPeucker(y,epsilon)
	# determine the leaf's tip
	pointe_index,leaflength = compute_leaf_point(y_smoothed)
	# split the leaf into the two halves defined by the leaftip
	zd,zg = split_leaf(y_smoothed,pointe_index)
	# on each half determine the remarkable points
	ind_kmin_d,ind_kmax_d = compute_remarkable_points(zd)
	kmin_d = evaluate(zd,ind_kmin_d)
	kmax_d = evaluate(zd,ind_kmax_d)
	ind_kmin_g,ind_kmax_g = compute_remarkable_points(zg)
	kmin_g = evaluate(zg,ind_kmin_g)
	kmax_g = evaluate(zg,ind_kmax_g)
	# reconstitute the lists of sinuses and tips for the entire leafcontour
	ind_kmin = glue_leaf(ind_kmin_d,ind_kmin_g,pointe_index)
	ind_kmax = glue_leaf(ind_kmax_d,ind_kmax_g,pointe_index)
	return y_smoothed,ind_kmin,ind_kmax

def find_local_minima(a):
	"""
	Finds local minima on a vector `a` and returns the list of indices of these minima.

	**Parameters**

	a : vector of real values

	**Returns**

	index_list : list of positive integers
		the list of indices of local minima of the vector `a`

	**Examples**

	>>> from contours import find_local_minima
	>>> import numpy as np
	>>> vect = np.array([1, 5, 7.8, 1e-7, 5,-8,-8,-8,6])
	>>> find_local_minima(vect)
	[0, 3, 6]

	"""
	index_list=[]
	if a[0]<a[1]:
		index_list.append(0)
	i = 1
	while i < len(a)-1:
		if a[i-1]>a[i] and a[i]<a[i+1] :
			index_list.append(i)
			i +=1
		else :
			if a[i-1]>a[i] and a[i]==a[i+1] :
				# if there is a plateau constituting a minimum
				# take as minimal value the middle of the plateau
				number = 1
				while a[i]==a[i+number] :
					number += 1
				if a[i]<a[i+number] :
					index_list.append(i+number/2)
					i += number+1
			else :
				i +=1
	if a[-1]<a[-2]:
		index_list.append(len(a)-1)
	return index_list

def compute_curvature_by_angles(z):
	"""
	Computes the curvature in each point of the curve `z`.
	The curvature in a point is approximated by the angle formed by the two
	consecutive segments of the curve which contain this point.

	The curvature values on the two extremities are taken to be zero.

	**Parameters**

	y : 2-column array
		defines a curve as an ordered list of points, 
		each row of the array represents a point on the curve and contains
		its 2D cartesian coordinates

	**Returns**

	k : vector of real values
		vector containing the angle values in each point of the curve. 
		If `y` is an (N,2) dimensional array, then `k` is an N dimensional vector.

	"""
	N = len(z)
	k = np.zeros(N)
	for i in range(1,N-1) :
		k[i] = compute_angle(z[i-1],z[i],z[i+1])-pi
	return k

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


def read_contourR(filename) :
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


def distance_point_droite(p,y):
	"""
	Gives the distance of the point `p` from
	the straight line defined by the first and last point in the point-list `y`.

	**Parameters**

	p : 2 dimensional vector
		cartesian coordinates of the point `p`

	y : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	**Returns**

	distance : positive real number
		the distance between the point and the straight line

	"""
	u = y[0]-y[-1]
	d = p-y[0]
	distance_up = abs(d[0]*u[1]-d[1]*u[0])
	distance = distance_up/norm(u)
	return distance

def find_most_distant_point(y):
	"""
	Given an ordered point-list y,
	the function determines the point of the list
	which is the most distant from the straight line
	defined by the first and last point in the list.

	**Parameters**

	y : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	**Returns**

	max_dist : positive real number
		the maximal distance

	max_index : positive integer
		the index of the most distant point


	"""
	dist = [distance_point_droite(y[i],y) for i in xrange(0,len(y))]
	max_index = argmax(dist)
	max_dist = dist[max_index]
	return max_dist, max_index

def DouglasPeucker(y,epsilon) :
	"""
	Given an ordered point-list y, which defines a curve,
	the function gives its simplification by the Douglas-Peucker (or `split-and-merge`) algorithme.

	The purpose of this recursive simplification algorithm is, to find a similar curve,
	which consists of a subset of the points that defined the original curve. 
	Similarity is controled by a maximal distance `epsilon` between the original
	and the simplified curves.

	**Parameters**

	y : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	**Returns**

	y_smoothed : 2-column array
		defines the output curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	"""
	# find the point on the curve with the maximum distance
	# from the segment defined by the curve's extremity points
	max_distance,max_index=find_most_distant_point(y)
	# If max_distance is greater than epsilon, recursively simplify
	if max_distance >= epsilon :
		arc1 = np.copy(DouglasPeucker(y[:max_index+1],epsilon))
		arc2 = np.copy(DouglasPeucker(y[max_index:],epsilon))
		# build the result list
		y_smoothed = np.array(list(arc1)+list(arc2[1:]))
	else :
		y_smoothed = np.array([y[0],y[-1]])
	return y_smoothed

def reinit_contour(y,thresh):
	"""
	There is a `gap` in a contour if there is a segment in its
	point-list definition, which is longer than a threshold `thresh`.
	In case of a gap, the contour is reinitialised so that the gap falls 
	at the end. This is useful to reinitialise the point-list in such a 
	way that the first and last point of the list marks the base of the leafblade.
	"""
	# identifies the discontinuity point as two points
	# being at distance which is superior to a threshold
	distances=[norm(y[i+1,:]-y[i,:]) for i in range(0,len(y)-1)]
	distances = [norm(y[0,:]-y[-1,:])]+distances
	a=[(i,distances[i]) for i in xrange(0,len(y)) if distances[i]>thresh]
	# we only treat the case when there is only ONE discontinuity under the given threshold !!!
	if len(a)==0 :
		print 'no discontinuity'
		z=y
	else :
		if len(a)>1 :
			print 'more than one discontinuity'
			index=-1
		else :
			print 'there is one discontinuity and the curve is reinitialised'
			index=a[0][0]
			# the discontinuity is between index-1 and index
		z=np.zeros_like(y)
		z[:len(y)-index]=y[index:]
		print "reinit_contour : index=",index
		print len(z[len(y)-index:]), len(y[:index])
		z[len(y)-index:]=y[:index]
	return z


def compute_optimal_epsilon(y,epsilon_step):
	"""
	the idea is to make a scan with a list of different DP parameters
	and choose the one which corresponds to the most robust detected 
	remarkable pointnumber

	**Parameters**

	y : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	epsilon_step : positive real number
		distance between two consecutive epsilon parameters in the parameter scan

	**Returns**

	epsilon : positive real number
		the optimal parameter value

	"""
	pointe_index,leaflength=compute_leaf_point(y)
	# the parameter values in the parameter scan
	epsilons=np.arange(2,int(leaflength/20),epsilon_step)
	# for all epsilon value we count the detected sinuses
	number_of_points=[]
	for epsilon in epsilons :
		print '-------'
		print 'epsilon=',epsilon
		y_DP = DouglasPeucker(y,epsilon)
		ind_kmin,ind_kmax = compute_remarkable_points(y_DP)
		number_of_points.append(len(ind_kmin))
	number_of_points=np.array(number_of_points)
	epsilons=np.array(epsilons)
	# find the relative length of plateaux
	plateau_length=[]
	plateau_epsilon=[]
	current_epsilon=epsilons[0]
	current_plateau_value=number_of_points[0]
	for i in xrange(1,len(epsilons)) :
		if number_of_points[i]<current_plateau_value or i==len(epsilons)-1:
			length=epsilons[i]-current_epsilon
			plateau_length.append(length/current_epsilon)
			plateau_epsilon.append(current_epsilon)
			current_epsilon=epsilons[i]
			current_plateau_value=number_of_points[i]
	plateau_length=np.array(plateau_length)
	plateau_epsilon=np.array(plateau_epsilon)
	# look for local maxima on plateau-length
	local_maxima_indices=find_local_minima(-plateau_length)
	# take the first one
	if len(local_maxima_indices)>1 :
		first_maxi_index=local_maxima_indices[1]
	else :
		first_maxi_index=local_maxima_indices[0]
	# the robust parametervalue
	epsilon=plateau_epsilon[first_maxi_index]
	return epsilon


def compute_optimal_epsilon_distfactscan(y, dist_factor_step):
	"""
	the idea is to make a scan with a list of different DP parameters
	and choose the one which corresponds to the most robust detected 
	remarkable pointnumber. However, the scan is done with respect to
	the "dist_factor", which is the ratio epsilon/leaflength.

	**Parameters**

	y : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	dist_factor_step : positive real number
		step in the dist_factor parameter scan

	**Returns**

	epsilon : positive real number
		the optimal parameter value

	"""
	pointe_index,leaflength=compute_leaf_point(y)
	# the parameter values in the parameter scan
	epsilon_step = leaflength*dist_factor_step
	epsilons=np.arange(epsilon_step*2,int(leaflength/20),epsilon_step)
	# for all epsilon value we count the detected sinuses
	number_of_points=[]
	for epsilon in epsilons :
		print '-------'
		print 'epsilon=',epsilon
		y_DP = DouglasPeucker(y,epsilon)
		ind_kmin,ind_kmax = compute_remarkable_points(y_DP)
		number_of_points.append(len(ind_kmin))
	number_of_points=np.array(number_of_points)
	epsilons=np.array(epsilons)
	# find the relative length of plateaux
	plateau_length=[]
	plateau_epsilon=[]
	current_epsilon=epsilons[0]
	current_plateau_value=number_of_points[0]
	for i in xrange(1,len(epsilons)) :
		if number_of_points[i]<current_plateau_value or i==len(epsilons)-1:
			length=epsilons[i]-current_epsilon
			plateau_length.append(length/current_epsilon)
			plateau_epsilon.append(current_epsilon)
			current_epsilon=epsilons[i]
			current_plateau_value=number_of_points[i]
	plateau_length=np.array(plateau_length)
	plateau_epsilon=np.array(plateau_epsilon)
	# look for local maxima on plateau-length
	local_maxima_indices=find_local_minima(-plateau_length)
	# take the first one
	if len(local_maxima_indices)>1 :
		first_maxi_index=local_maxima_indices[1]
	else :
		first_maxi_index=local_maxima_indices[0]
	# the robust parametervalue
	epsilon=plateau_epsilon[first_maxi_index]
	return epsilon



################ hierarchy ##################

class leaf_lobe:
	"""
	A leaf_lobe is a portion of a leaf-contour, defined by two
	(not necessarily consecutive) sinuses.

	**Attributes**

	deg : positive integer
		degree of the lobe

	indexes : list of two indexes
		the indices of the first and last sinus of the lobe

	z : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	ind_sin : list of positive integers
		the indices of all sinuses on the curve

	"""
	def __init__(self, degree,lobeindexes,z,ind_sin):
		self.deg=degree
		self.indexes=lobeindexes
		self.z=z
		self.ind_sin=ind_sin
	#
	def cuttable(self) :
		# the lobe is cuttable, if there is at least one
		# intermediate sinus
		lobe=self.indexes
		sinuses_in_lobe=[i for i in self.ind_sin if (i>lobe[0] and i<lobe[1])]
		return len(sinuses_in_lobe)>0
	#
	def basewidth(self) :
		# the distance between the two extremities of the lobe
		lobe=self.indexes
		z=self.z
		return norm(z[lobe[0]]-z[lobe[1]])
	#
	def cut_lobe_old(self) :
		# if the lobe is cuttable, it creates two lobes.
		# The created lobes might be sisters, and in this case both have the same degree 
		# as the original lobe. The other situation is 
		#alpha_threshold=1.2*pi
		alpha_threshold=1.2*pi
		if self.cuttable()==0 :
			toto = [self]
		else :
			deg=self.deg
			lobe=self.indexes
			z=self.z
			sinuses_in_lobe=[i for i in self.ind_sin if (i>lobe[0] and i<lobe[1])]
			#
			alpha_min=2*pi
			for i in sinuses_in_lobe :
				alpha=compute_angle(z[lobe[0]],z[i],z[lobe[1]])
				if alpha<alpha_min :
					alpha_min=alpha
					i_cut=i
			lobe1=[lobe[0],i_cut]
			lobe2=[i_cut,lobe[1]]
			deg1=deg
			deg2=deg
			# if alpha_min>alpha_threshold one of the teeth is secondary
			# the secondary tooth is the one which has the lower median
			if alpha_min>alpha_threshold :
				lobe_mid=(lobe[0]+lobe[1])/2.0
				tip1=farest(z[lobe1[0]:lobe1[1]],lobe_mid)+lobe1[0]
				l1=norm(z[tip1]-lobe_mid)
				tip2=farest(z[lobe2[0]:lobe2[1]],lobe_mid)+lobe2[0]
				l2=norm(z[tip2]-lobe_mid)
				if l1>l2 :
					deg2=deg+1
				else :
					deg1=deg+1
			lobe_object1=leaf_lobe(deg1,lobe1,z,self.ind_sin)
			lobe_object2=leaf_lobe(deg2,lobe2,z,self.ind_sin)
			toto=[lobe_object1, lobe_object2]
		return toto


	def cut_lobe(self) :
		# if the lobe is cuttable, it creates two lobes.
		# The created lobes might be sisters, and in this case both have the same degree 
		# as the original lobe. The other situation is 
		#alpha_threshold=1.2*pi
		alpha_threshold=1.25*pi
		if self.cuttable()==0 :
			toto = [self]
		else :
			deg=self.deg
			lobe=self.indexes
			z=self.z
			sinuses_in_lobe=[i for i in self.ind_sin if (i>lobe[0] and i<lobe[1])]
			#
			alpha_min=2*pi
			for i in sinuses_in_lobe :
				alpha=compute_angle(z[lobe[0]],z[i],z[lobe[1]])
				if alpha<alpha_min :
					alpha_min=alpha
					i_cut=i
			lobe1=[lobe[0],i_cut]
			lobe2=[i_cut,lobe[1]]
			deg1=deg
			deg2=deg
			# if alpha_min>alpha_threshold one of the teeth is secondary
			# the secondary tooth is the one which has the lower median
			if alpha_min>alpha_threshold :
				#deg1, deg2 = compute_lobe_degree_using_median(lobe1, lobe2, deg, z)
				deg1, deg2 = compute_lobe_degree_using_areas(lobe1, lobe2, deg, z)
			lobe_object1=leaf_lobe(deg1,lobe1,z,self.ind_sin)
			lobe_object2=leaf_lobe(deg2,lobe2,z,self.ind_sin)
			toto=[lobe_object1, lobe_object2]
		return toto

def compute_lobe_degree_using_median(lobe1, lobe2, deg, z) :
	deg1 = deg
	deg2 = deg
	lobe_mid=(z[lobe1[0]]+z[lobe2[1]])/2.0
	tip1=farest(z[lobe1[0]:lobe1[1]],lobe_mid)+lobe1[0]
	l1=norm(z[tip1]-lobe_mid)
	tip2=farest(z[lobe2[0]:lobe2[1]],lobe_mid)+lobe2[0]
	l2=norm(z[tip2]-lobe_mid)
	if l1>l2 :
		deg2=deg+1
	else :
		deg1=deg+1
	return deg1, deg2


def compute_lobe_degree_using_areas(lobe1, lobe2, deg, z) :
	"""
	Given a lobe of degree "deg" on a curve "z", 
	which is cut into two consecutive lobes having their list of starting 
	and ending indices "lobe1" and "lobe2" respectively. We already know
	that one of the emerging lobe is the daughter of the other.
	This function gives the degrees of the two emerging lobes based
	on a criterion of the daughter/mother area ratio.
	"""
	deg1 = deg
	deg2 = deg
	# T is the area of the triangle defined by the points lobe1[0], lobe1[1] and lobe2[1]
	x1 = z[lobe1[0]]-z[lobe1[1]]
	x2 = z[lobe2[1]]-z[lobe1[1]]
	T = abs(x1[0]*x2[1]-x2[0]+x1[1])/2.0
	# A1 is the area under lobe1
	A1 = compute_area_of_lobe(z[lobe1[0]:lobe1[1]])
	# A2 is the area under lobe2
	A2 = compute_area_of_lobe(z[lobe2[0]:lobe2[1]])
	# Ri is the daughter/mother area ratio if i is the mother
	R1 = A2/(A1+T)
	R2 = A1/(A2+T)
	if R1<R2 :
		deg2=deg+1
	else :
		deg1=deg+1
	return deg1, deg2



def treat_node(G,node_nb,lobe_dict):
	"""
	This is the main function of the hierarchy algorithm.
	Its structure is analogous to the Douglas-Peucker algorithm in the 
	sense that we iteratively cut and recut lobes.

	**Parameters**

	G : an oriented graphe
		contains lobes as nodes, while the edges are such that they point 
		from mother to daughter lobes.

	node_nb : positive integer
		the label of the node to treat

	lobe_dict : dictionnary (integer, leaf_lobe)
		for all node-label of the graphe (key) defines the leaf_lobe which corresponds to it.

	**Returns**

	"""
	lobe=lobe_dict[node_nb]
	lobe_list=lobe.cut_lobe()
	if len(lobe_list)==2 :
		l1=lobe_list[0]
		l2=lobe_list[1]
		if l1.deg==l2.deg :
			# l1 replaces the original
			lobe_dict[node_nb]=l1
			# l2 is added
			i2=len(lobe_dict)
			lobe_dict[i2]=l2
			pred=G.predecessors(node_nb)
			G.add_edge(pred[0],i2)
			treat_node(G,node_nb,lobe_dict)
			treat_node(G,i2,lobe_dict)
		else :
			if l2.deg>l1.deg :
				# l1 remains in the graphe (with modified indices)
				# and we append l2 as a daughter of l1
				lobe_dict[node_nb]=l1
				i2=len(lobe_dict)
				lobe_dict[i2]=l2
				G.add_edge(node_nb,i2)
				treat_node(G,node_nb,lobe_dict)
				treat_node(G,i2,lobe_dict)
			else :
				# l2 remains and we append l1 as a daughter of l2
				lobe_dict[node_nb]=l2
				i1=len(lobe_dict)
				lobe_dict[i1]=l1
				G.add_edge(node_nb,i1)
				treat_node(G,node_nb,lobe_dict)
				treat_node(G,i1,lobe_dict)
	return


def create_lobe_graphe(zd,ind_sin) :
	'''
	The function constructs the hierarchised structure of the half-leaf 
	- represented by a tree (oriented graphe)
	- the root of the graphe is the half-leaf - degree 0
	- primary teeth are daughters of the root - degree 1
	- secondary teeth are daughters of the tooth on which they live - degree 2

	**Parameters**

	zd : 2-column array
		defines the input curve as an ordered list of points, 
		each row of the array representing a point on the curve and containing
		its 2D cartesian coordinates

	ind_sin : list of positive integers
		the list of indices of sinuses on the input contour

	**Returns**

	G : an oriented graphe
		represents the hierarchy of the teeths.

	lobe_dictionnary : dictionnary (integer, leaf_lobe)
		the dictionnary which for all node-label of the graphe (key) 
		gives the corresponding leaf_lobe.

	'''
	# initialisation of the graphe
	# and of the dictionnary for the lobe <-> node-label correspondence
	G = nx.DiGraph()
	lobe_dictionnary={}
	# define the half of the leaf as a degree 0 lobe
	lobe0=leaf_lobe(0,[0,len(zd)-1],zd,ind_sin)
	# add it to the graphe and to the dictionnary (with label 0)
	G.add_node(0)
	lobe_dictionnary[0]=lobe0
	# if there is at least one tooth on the half-leaf
	if len(ind_sin)>1 :
		# define the lobe defined by the first and last sinus as a degree 1 lobe
		lobe0=leaf_lobe(1,[ind_sin[0],ind_sin[-1]],zd,ind_sin)
		# add it to the graphe and to the dictionnary (with label 1)
		G.add_edge(0,1)
		lobe_dictionnary[1]=lobe0 
		treat_node(G,1,lobe_dictionnary)
	return G, lobe_dictionnary


def find_graphe_labels_for_tips(G,lobedict,ind_kmax):
	"""
	Supposing that there is only one tip on each lobe, the function
	finds for each tip its label in the graphe `G`
	"""
	tip_dict={}
	teeth_nodes=G.nodes()
	teeth_nodes.remove(0)
	for i in teeth_nodes:
		lobe=lobedict[i].indexes
		tip=[j for j in ind_kmax if (lobe[0]<j and j<lobe[1])]
		tip_dict[tip[0]]=i
	return tip_dict


def compute_area_of_lobe(z):
	"""
	Supposing that there is only one tip on each lobe, the function
	finds for each tip its label in the graphe `G`
	"""
	x1 = z[0]
	x2 = z[-1]
	x = z - x2
	A = -(x1[0]-x2[0])*(x1[1]-x2[1])/2.0
	for i in xrange(1,len(z)) :
		A -= (x[i,1]+x[i-1,1])*(x[i,0]-x[i-1,0])/2.0 
	return A


################################

def locmin(a):
	m=np.r_[True,a[1:]<a[:-1]]&np.r_[a[:-1]<=a[1:],True]
	return np.int0(m)

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
