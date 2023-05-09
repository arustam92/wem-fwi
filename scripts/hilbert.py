import SepVector
import Hypercube
import genericIO
import Operator
import pyOperator as Op
import json
import sys
import numpy as np

for par in sys.argv[1:]:
    if ('in=' in par):
        input = par.split('in=')[1]
    if ('out=' in par):
        output = par.split('out=')[1]

model = genericIO.defaultIO.getVector(input)
modelNd = model.getNdArray()
ax = model.getHyper().axes

if len(modelNd.shape) > 3:
    full = model.clone()
    fullNd = full.getNdArray()
    data = SepVector.getSepVector(Hypercube.hypercube(ns=[ax[0].n,ax[1].n,ax[2].n],
                                            os=[ax[0].o,ax[1].o,ax[2].o], ds=[ax[0].d,ax[1].d,ax[2].d]),storage='dataComplex')
    dataNd = data.getNdArray()
    mod = data.clone()
    modNd = mod.getNdArray()
    hilbert = Operator.Hilbert(mod,data)

    for i in range(model.shape[0]):
        print(i)
        modNd[:] = modelNd[i,:]
        hilbert.forward(False,mod,data)
        fullNd[i,:] = dataNd[:]
else:
    full = model.clone()
    hilbert = Operator.Hilbert(model,full)
    hilbert.forward(False,model,full)

full.writeVec(output)
