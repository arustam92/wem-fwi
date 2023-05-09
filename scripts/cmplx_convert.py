import SepVector
import Hypercube
import genericIO
import WEM
import Operator
import sys

for par in sys.argv[1:]:
    if ('in' in par):
        input = par.split('in=')[1]
    if ('out' in par):
        out = par.split('out=')[1]

data = genericIO.defaultIO.getVector(input)

n = []
d = []
o = []
for i in range(data.getHyper().getNdim()):
    n.append(data.getHyper().getAxis(i+1).n)
    d.append(data.getHyper().getAxis(i+1).d)
    o.append(data.getHyper().getAxis(i+1).o)

d[0] = float(1/(n[0]*d[0]))

dataF = SepVector.getSepVector(Hypercube.hypercube(ns=n,os=o,ds=d),storage="dataComplex")


fftOp = Operator.FFT(data,dataF,1,1)
fftOp.forward(False,data,dataF)

genericIO.defaultIO.writeVector(out,dataF)
