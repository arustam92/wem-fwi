import Conformal
import genericIO
import SepVector
import Hypercube
import numpy as np
import WEM
import json

"""Create velocity model"""
path = '/home/user/project/projects/marmousiWEM'
vel = genericIO.defaultIO.getVector(path + '/Vel/velTrue.H')
velNd = vel.getNdArray()

n1 = vel.getHyper().getAxis(1).n
d1= vel.getHyper().getAxis(1).d
n2 = vel.getHyper().getAxis(2).n
d2 = vel.getHyper().getAxis(2).d
o1 = vel.getHyper().getAxis(1).o
o2 = vel.getHyper().getAxis(2).o
slow = SepVector.getSepVector(Hypercube.hypercube(hypercube=vel.getHyper()),storage='dataComplex')
slowNd = slow.getNdArray()
slowNd[:] = 1e-3/velNd[:]
#
# map = Conformal.ComplexExp(slow.getHyper(),eps=0.01)
# mappedSlow = map.forward(slow)
# mappedSlow.writeVec('marm_mapped.H')

"""Read in parameters"""
with open('geometry.json') as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()

"""Read in wavelet"""
wave = genericIO.defaultIO.getVector(path + '/Wav/MarmWave_yinbin.H')

"""Create wavefield"""
nt = wave.getHyper().getAxis(1).n
dt = wave.getHyper().getAxis(1).d
ot = 0.
oshot = 0
dshot = 0
nshot = parObj.getInt('ns')

data = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,n1,n2,nshot],
													os=[ot,o1,o2,oshot],
													ds=[dt,d1,d2,dshot]))

wem = WEM.WEM(slow,data,parObj,wave)
wem.wavefield(slow,data)

genericIO.defaultIO.writeVector('marm_wavefield_two.H',data)
