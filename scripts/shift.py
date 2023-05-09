import SepVector
import Hypercube
import genericIO
import WEM
import newExtWEM
import json
import sys
import numpy as np
import Operator
import pyOperator as Op

for par in sys.argv[1:]:
    if ('shift' in par):
        shift = int(par.split('shift=')[1])
    if ('in' in par):
        input = par.split('in=')[1]
        input = genericIO.defaultIO.getVector(input)
    if ('out' in par):
        out = par.split('out=')[1]

output = input.clone()
waveNd = input.getNdArray()
outNd = output.getNdArray()
outNd[:] = np.roll(waveNd,shift)

output.writeVec(out)
