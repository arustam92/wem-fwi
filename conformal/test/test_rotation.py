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
omax1 = 0.5
o2 = 0.
omax2 = 0.5
d1= (omax1-o1)/n1
d2 = (omax2-o2)/n2

input = SepVector.getSepVector(Hypercube.hypercube(ns=[n1,n2],ds=[d1,d2],os=[o1,o2]), storage='dataComplex')

inputNd = input.getNdArray()
inputNd[0:250,:] = 1
inputNd[250:,:] = -1

angles = np.linspace(np.pi/2, -np.pi/2, 256)

out = SepVector.getSepVector(Hypercube.hypercube(ns=[n1,n2,angles.size],ds=[d1,d2,(angles[1]-angles[0])],os=[o1,o2,angles[0]]), storage='dataComplex')
outNd = out.getNdArray()

for i in range(angles.size):
	print(i)
	rot = Conformal.Rotation(input.getHyper(),angles[i])
	output = rot.forward(input)
	outputNd = output.getNdArray()
	outNd[i,:,:] = outputNd[:]


#
out.writeVec('output.H')
input.writeVec('input.H')
#
#
# """
# 	Check inverse
# """
# outputNd = output.getNdArray()
# for i in range(outputNd.shape[0]):
# 	outputNd[:,i] = i
#
# start = time.time()
# input = cexp.inverse(output)
# print('Inverse: %f' %(time.time() - start))
# input.writeVec('input.H')
