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
    if ('wave=' in par):
        wave = par.split('wave=')[1]
        wave = genericIO.defaultIO.getVector(wave)
        nt = wave.getHyper().getAxis(1).n
    if ('nt=' in par):
        nt = int(par.split('nt=')[1])
    if ('f0=' in par):
        f0 = float(par.split('f0=')[1])

input = genericIO.defaultIO.getVector(input)
inputNd = input.getNdArray()
waveNd = wave.getNdArray()

wax = wave.getHyper().axes
iax = input.getHyper().axes

outFreq = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,iax[0].n,iax[1].n],
                                        os=[wax[0].o,iax[0].o,iax[1].o], ds=[wax[0].d,iax[0].d,iax[1].d]),storage='dataComplex')
outTau = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,iax[0].n,iax[1].n],
                                        os=[wax[0].o,iax[0].o,iax[1].o], ds=[wax[0].d,iax[0].d,iax[1].d]))
outFreqNd = outFreq.getNdArray()
outTauNd = outTau.getNdArray()


iter = 1
if len(input.shape) > 3:
    iter = input.shape[0]
    full = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,iax[0].n,iax[1].n,iter],
                                            os=[wax[0].o,iax[0].o,iax[1].o,0], ds=[wax[0].d,iax[0].d,iax[1].d,1]),storage='dataComplex')
    fullNd = full.getNdArray()
    for i in range(iter):
        print(i)
        outTauNd[:] = inputNd[i,:]
        fft = Operator.FFT(outTau,outFreq,1,1)
        fft.forward(False,outTau,outFreq)
        fullNd[i,:] = outFreq[:]
else:
    full = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,iax[1].n,iax[2].n],
                                            os=[wax[0].o,iax[1].o,iax[2].o], ds=[wax[0].d,iax[1].d,iax[2].d]),storage='dataComplex')
    fft = Operator.FFT(input,full,1,1)
    fft.forward(False,input,full)

full.writeVec(output)
