import Conformal
import genericIO
import SepVector
import Hypercube
import numpy as np
import time

n1 = 100
n2 = 100
o1 = -1
omax1 = 100
o2 = -1
omax2 = 100
d1= float(omax1-o1)/n1
d2 = float(omax2-o2)/n2
max1 = o1 + (n1-1)*d1
max2 = o2 + (n2-1)*d2
diag = max2 + 1j*max1

input = SepVector.getSepVector(Hypercube.hypercube(ns=[n1,n2],ds=[d1,d2],os=[o1,o2]), storage='dataComplex')
"""
	Create lines
"""
inputNd = input.getNdArray()
# inputNd[:,50] = -1
for i in range(0,inputNd.shape[0],2):
	inputNd[i,:] = 1
	inputNd[:,i] = -1

start = time.time()
a = 1/np.abs(diag)
b = 80+50j
sas = Conformal.ShiftAndScale(input.getHyper(),a=a, b=b)
print('Constructor 1: %f' % (time.time() - start))

start = time.time()
cexp = Conformal.ComplexLog(sas.getOutHyper(),eps=0.02)
print('Constructor 2: %f' % (time.time() - start))

start = time.time()
chain = Conformal.MapChain(sas,cexp)
print('Constructor 3: %f' % (time.time() - start))

start = time.time()
output = chain.forward(input)
print('Forward: %f' % (time.time() - start))
output.writeVec('output.H')


"""
	Check inverse
"""
outputNd = output.getNdArray()
for i in range(outputNd.shape[0]):
	outputNd[:,i] = i

start = time.time()
input = chain.inverse(output)
print('Inverse: %f' %(time.time() - start))
input.writeVec('input.H')
