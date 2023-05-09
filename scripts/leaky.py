import SepVector
import Hypercube
import genericIO
import Operator
import pyOperator as Op
import json
import sys

double = 0
for par in sys.argv[1:]:
    if ('in=' in par):
        input = par.split('in=')[1]
    if ('out=' in par):
        out = par.split('out=')[1]
    if ('alpha=' in par):
        alpha = float(par.split('alpha=')[1])
        alpha_final=alpha
    if ('alpha_final=' in par):
        alpha_final = float(par.split('alpha_final=')[1])
    if ('double=' in par):
        double = float(par.split('double=')[1])

slow = genericIO.defaultIO.getVector(input)
slowNd = slow.getNdArray()
output = slow.clone()
outputNd = output.getNdArray()

ax = slow.getHyper().axes
temp = SepVector.getSepVector(Hypercube.hypercube(ns=[ax[0].n,ax[1].n,ax[2].n],
													os=[ax[0].o, ax[1].o, ax[2].o],
													ds=[ax[0].d,ax[1].d,ax[2].d]),storage='dataComplex')
tempNd = temp.getNdArray()

if double:
    leak = Operator.DoubleLeakyIntegration(temp,temp,alpha=alpha,alpha_final=alpha_final)
else:
    leak = Operator.LeakyIntegration(temp,temp,alpha=alpha,alpha_final=alpha_final)

if len(slowNd.shape) < 4:
    leak.forward(False,slow,output)
else:
    temp2 = temp.clone()
    temp2nd = temp2.getNdArray()
    for iter in range(slowNd.shape[0]):
        print(iter)
        temp2nd[:] = slowNd[iter,:]
        leak.forward(False,temp2,temp)
        outputNd[iter,:] = tempNd[:]

output.writeVec(out)
