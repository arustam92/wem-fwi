import pyCMap
import SepVector
import Hypercube
import numpy as np

class ConformalMap:
	def __init__(self,hyper):
		self.cppMode = None

	def forward(self,inSpace):
		return SepVector.floatVector(
		fromCpp=self.cppMode.forward(inSpace.getCpp()))

	def inverse(self,outSpace):
		return SepVector.floatVector(
		fromCpp=self.cppMode.inverse(outSpace.getCpp()))

	def getInHyper(self):
		return Hypercube.hypercube(
		hypercube=self.cppMode.getInHyper())
	def getOutHyper(self):
		return Hypercube.hypercube(
		hypercube=self.cppMode.getOutHyper())



class ComplexLog (ConformalMap):
	def __init__(self, hyper, eps=0.01):
		self.cppMode = pyCMap.ComplexLog(hyper.getCpp(),eps)

class ShiftAndScale (ConformalMap):
	def __init__(self, hyper, a=1, b=0):
		self.cppMode = pyCMap.ShiftAndScale(hyper.getCpp(), a, b)

class MapChain (ConformalMap):
	def __init__(self, map1, map2):
		self.cppMode = pyCMap.MapChain(map1.cppMode,map2.cppMode)

class Rotation (ConformalMap):
	def __init__(self, hyper, rad):
		self.cppMode = pyCMap.Rotation(hyper.cppMode,rad)
