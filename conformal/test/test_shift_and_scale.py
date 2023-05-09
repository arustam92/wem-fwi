import Conformal
import genericIO
import SepVector
import Hypercube
import numpy as np
import time

"""
	Input has to be a unit circle
	Output maps to [0,2pi] x [-eps,0]
"""
n1 = 500
n2 = 500
o1 = -0.8
omax1 = 0.2
o2 = -0.2
omax2 = 0.5
d1= (omax1-o1)/n1
d2 = (omax2-o2)/n2

input = SepVector.getSepVector(Hypercube.hypercube(ns=[n1,n2],ds=[d1,d2],os=[o1,o2]), storage='dataComplex')
a = 2
b = 1j
inputNd = input.getNdArray()
inputNd[:] = 1

start = time.time()
sas = Conformal.ShiftAndScale(input.getHyper(),a=a, b=b)
print('Constructor: %f' % (time.time() - start))

start = time.time()
output = sas.forward(input)
print('Forward: %f' % (time.time() - start))
output.writeVec('output.H')


"""
	Check inverse
"""
start = time.time()
input = sas.inverse(output)
print('Inverse: %f' %(time.time() - start))
input.writeVec('input.H')
