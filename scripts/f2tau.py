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


if0 = int(f0 * (wax[0].n-1)*wax[0].d)
iter = 1
if len(input.shape) > 3:
    iter = input.shape[0]
    full = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,iax[0].n,iax[1].n,iter],
                                            os=[wax[0].o,iax[0].o,iax[1].o,0], ds=[wax[0].d,iax[0].d,iax[1].d,1]))
    fullNd = full.getNdArray()
    for i in range(iter):
        print(i)
        # fill with average value first
        tmp = np.transpose(inputNd[i,:],(1,2,0))
        # for j in range(outFreqNd.shape[2]):
        #     outFreqNd[:,:,j] = np.sum(tmp,axis=2) / tmp.shape[2] * np.sqrt(wax[0].n)
        # outFreqNd[:,:,:if0] = tmp[:,:,0]
        # outFreqNd[:,:,inputNd.shape[1]+if0:] = tmp[:,:,-1]
        # get the actual values
        outFreqNd[:,:,if0:inputNd.shape[1]+if0] = tmp[:]

        fft = Operator.FFT(outTau,outFreq,1,1)
        fft.adjoint(False,outTau,outFreq)
        fullNd[i,:] = np.fft.fftshift(outTauNd,axes=2)
else:
    full = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,iax[0].n,iax[1].n],
                                            os=[wax[0].o,iax[0].o,iax[1].o], ds=[wax[0].d,iax[0].d,iax[1].d]))
    fullNd = full.getNdArray()
    outFreqNd[:,:,if0:inputNd.shape[0]+if0] = np.transpose(inputNd[:],(1,2,0))
    fft = Operator.FFT(outTau,outFreq,1,1)
    fft.adjoint(False,outTau,outFreq)
    fullNd[:] = np.fft.fftshift(outTauNd,axes=2)

full.writeVec(output)
