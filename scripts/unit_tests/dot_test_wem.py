import SepVector
import Hypercube
import genericIO
import pyOperator as Op
import newExtWEM
import WEM
import Operator
import json
import os
import sys
import numpy as np

ns = [300,300]
ds = [10,10]
os = [0,0]

nt = 1500
dt = 0.002
ot = 0
f0 = 5

vel = 2000

par = {
    "ns" : 1,
    "osx" : os[0],
    "dsx" : 20,
    "osz" : os[1],

    "nr" : ns[0],
    "orx" : os[0],
    "drx" : ds[0],
    "orz" : os[1],

    "prop" : "ssf",
    "nref" : 1,
    "tap" : 50,
    "fmin" : 1,
    "fmax" : 15,
    "ng" : 0,
    "ntaylor" : 1,

    "onepass" : 0,
}

t = np.linspace(ot,(nt-1)*dt,nt)
ricker = (1-2*(np.pi*f0*t)**2)*np.exp(-(np.pi*f0*t)**2)

slow = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
data = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,par["nr"],par["ns"]],
                                                ds=[dt,par["drx"],par["dsx"]],
                                                os=[ot,par["orx"],par["osx"]]))
wave = SepVector.getSepVector(Hypercube.hypercube(ns=[nt],ds=[dt],os=[ot]))
parObj = genericIO.pythonParams(par).getCpp()

dataRand = data.clone()
waveRand = wave.clone()

slowNd = slow.getNdArray()
slowNd[:] = 1 / vel
dataNd = data.getNdArray()
dataRandNd = dataRand.getNdArray()
waveNd = wave.getNdArray()
waveNd[:] = ricker[:]

waveRandNd = waveRand.getNdArray()

waveRandNd[:] = np.random.rand(waveRandNd.shape[0])
dataRandNd[:] = np.random.rand(dataRandNd.shape[0],dataRandNd.shape[1],dataRandNd.shape[2])

wemOp = WEM.WEM(slow,data,wave,parObj)


wemOp.forward(False,waveRand,data)
# wemOp.adjoint(False,image,dataRand)

print("Model: %f" % imageRand.dot(image))
print("Data: %f" % dataRand.dot(data))
print("Abs. Error: %f" % (imageRand.dot(image)-dataRand.dot(data)))
print("Rel. Error: %f" % (1-imageRand.dot(image)/dataRand.dot(data)))


data.writeVec("test.H")



#
