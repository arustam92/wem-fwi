import SepVector
import Hypercube
import genericIO
import numpy as np
import Operator
import sys

inp = str(sys.argv[1])
bg = str(sys.argv[2])
output = str(sys.argv[3])

inp = genericIO.defaultIO.getVector(inp)
bg = genericIO.defaultIO.getVector(bg)

inpNd = inp.getNdArray()
bgNd = bg.getNdArray()
n = []
d = []
o = []
for i in range(inp.getHyper().getNdim()):
    n.append(inp.getHyper().getAxis(i+1).n)
    d.append(inp.getHyper().getAxis(i+1).d)
    o.append(inp.getHyper().getAxis(i+1).o)

ndim = len(n)
if ndim > 2:
    del n[2], o[2], d[2]

out = SepVector.getSepVector(Hypercube.hypercube(ns=n,os=o,ds=d))
outNd = out.getNdArray()

outNd[:] = np.sum(inpNd[:] - bgNd[:],axis=ndim-3) + bgNd[:]

genericIO.defaultIO.writeVector(output,out)
