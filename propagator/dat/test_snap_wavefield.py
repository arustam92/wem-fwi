import genericIO
import SepVector
import Hypercube
import numpy as np
import WEM
import json

"""Create velocity model"""
n1 = 100
n2 = 100
o1 = -100
omax1 = 1000
o2 = -100
omax2 = 1000
d1= float(omax1-o1)/n1
d2 = float(omax2-o2)/n2

slow = SepVector.getSepVector(Hypercube.hypercube(ns=[n1,n2],ds=[d1,d2],os=[o1,o2]),storage='dataComplex')
slowNd = slow.getNdArray()
slowNd[:] = 1 / 1000

"""Read in parameters"""
with open('geometry.json') as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()

"""Read in wavelet"""
wave = genericIO.defaultIO.getVector('wave.H')

"""Create wavefield"""
nt = wave.getHyper().getAxis(1).n
dt = wave.getHyper().getAxis(1).d
ot = 0.
oshot = parObj.getFloat('osx')
dshot = parObj.getFloat('dsx')
nshot = parObj.getInt('ns')

data = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,n1,n2,nshot],
													os=[ot,o1,o2,oshot],
													ds=[dt,d1,d2,dshot]))

wem = WEM.WEM(slow,data,wave,parObj)
wem.wavefield(slow,data)

genericIO.defaultIO.writeVector('const_wavefield_two.H',data)
