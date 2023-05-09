import SepVector
import Hypercube
import genericIO
import Operator
import pyOperator as Op
import json
import sys

for par in sys.argv[1:]:
    if ('in=' in par):
        input = par.split('in=')[1]
    if ('space=' in par):
        space = par.split('space=')[1]
out = list(input)
out.insert(-2,'_interpol')
out = ''.join(out)

slow = genericIO.defaultIO.getVector(input)
space = genericIO.defaultIO.getVector(space)

slowNd = slow.getNdArray()
output = space.clone()
outputNd = output.getNdArray()

ax = slow.getHyper().axes
ax2 = space.getHyper().axes
ns = []
os = []
ds = []
for i in ax:
    ns.append(i.n)
    os.append(i.o)
    ds.append(i.d)

temp = SepVector.getSepVector(Hypercube.hypercube(ns=ns, os=os, ds=ds),storage='dataComplex')
tempNd = temp.getNdArray()
ns[-1] = ax2[3].n
print(ns)
#
leak = Operator.CosineInterpolation(temp,temp,points=points)
if len(slowNd.shape) < 5:
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
