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


import numpy as np
from numpy import linalg


class point:
	def __init__(self):
		self.Index_Vertex = None
		self.nbTriangle = None
		self.nbInTr = None
		self.coords = [None, None]
		self.displ = [None,None] #ux, uy
		self.Elast = None
		self.RotatedMatrix = [None, None, None, None, None, None] #AR, BR, IR, GR, ER, FR
		self.HookeMatrix = [None, None, None, None] #A1, A2, B12, C3
		self.Deformation = [None, None, None, None] # dxux, dxuy, dyux, dyuy

	def initiate(self, row):
		self.setIndex(int(row[0]))
		self.setnbTr(int(row[1]))
		self.setnbInTr(int(row[2]))
		self.setCoords(tofloat(row[3:5]))
		self.setDispl(tofloat(row[5:7]))
		self.setElast(float(row[7]))
		self.setRotMat(tofloat(row[8:14]))
		self.setHookeMat(tofloat(row[14:18]))
		self.setDeformation(tofloat(row[19:23]))

	def getIndex(self):
		return self.Index_Vertex

	def setIndex(self, i):
		self.Index_Vertex = i

	def getnbTr(self):
		return self.nbTriangle

	def setnbTr(self, i):
		self.nbTriangle = i

	def getnbInTr(self):
		return self.nbInTr

	def setnbInTr(self, i):
		self.nbInTr = i

	def getCoords(self): # x, y
		return self.coords

	def setCoords (self, c):
		self.coords = c

	def getDispl(self): # ux, uy
		return self.displ

	def setDispl(self, u):
		self.displ = u

	def getElast(self): # Young's modulus
		return self.Elast

	def setElast(self, Y):
		self.Elast = Y

	def getRotMat(self):
		return self.RotatedMatrix

	def setRotMat(self, mat):
		self.RotatedMatrix = mat

	def getHookeMat(self):
		return self.HookeMatrix

	def setHookeMat(self, mat):
		self.HookeMatrix = mat

	def getDeformation(self):
		return self.Deformation

	def setDeformation(self, a):
		self.Deformation = a

	def getGrArea(self):
		'''
		returns a scalar (for a scalar plot)
		1+dux/dx+duy/dy = growth rate vertice based
		'''
		Def = self.getDeformation()
		Def = Sym(Def)
		VP = Diag(Def)
		return 1+VP[0]+VP[1]

	def getOne(self):
		return 1



def Diag(Matrix):
	'''
	Returns [eigenVal1, egienVal2] eigenVal1 > egenVal2.
	@ Argument can be:
	- 4 values array [a,b,c,d] for the matrix [a,b;c,d]
	- a 2-dimensional array [[a,b],[c,d]]
	- 3 values array, [a,b,c] for the matrix [a,c;c,b] (symetric matrix)
	'''
	if len(Matrix) == 3:
		Matrix = [Matrix[0], Matrix[2], Matrix[2], Matrix[1]]
	if len(Matrix)==4:
		tr = Matrix[0]+Matrix[3]
		det = Matrix[0]*Matrix[3]-Matrix[1]*Matrix[2]
	elif Matrix.shape == (2,2):
		tr = Matrix[0,0]+Matrix[1,1]
		det = Matrix[0,0]*Matrix[1,1]-Matrix[0,1]*Matrix[1,0]
	# 0=lambda^2 - tr*lambda + det
	delta = np.sqrt(tr**2 - 4.*det)
	VP = [0,0]
	VP[0] = (tr+delta)/2.
	VP[1] = (tr-delta)/2.
	return VP

def Sym(Matrix):
	Matrix[1] = (Matrix[1]+Matrix[2])/2.
	Matrix[2] = Matrix[1]
	return Matrix

def tofloat(liste):
	result = []
	for l in liste:
		result.append(float(l))
	return result
