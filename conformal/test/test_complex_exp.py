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
o1 = 0.
omax1 = 0.2
o2 = 0.
omax2 = 0.5
d1= (omax1-o1)/n1
d2 = (omax2-o2)/n2

input = SepVector.getSepVector(Hypercube.hypercube(ns=[n1,n2],ds=[d1,d2],os=[o1,o2]), storage='dataComplex')

start = time.time()
cexp = Conformal.ComplexLog(input.getHyper(),eps=0.02)
print('Constructor: %f' % (time.time() - start))
"""
	Create lines
"""
inputNd = input.getNdArray()
# inputNd[:,50] = -1
for i in range(0,inputNd.shape[0],2):
	inputNd[i,:] = 1
	inputNd[:,i] = -1

inputNd[:] = 1

start = time.time()
output = cexp.forward(input)
print('Forward: %f' % (time.time() - start))

output.writeVec('output.H')


"""
	Check inverse
"""
outputNd = output.getNdArray()
for i in range(outputNd.shape[0]):
	outputNd[:,i] = i

start = time.time()
input = cexp.inverse(output)
print('Inverse: %f' %(time.time() - start))
input.writeVec('input.H')
