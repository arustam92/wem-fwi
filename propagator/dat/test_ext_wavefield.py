import Conformal
import genericIO
import SepVector
import Hypercube
import numpy as np
import newExtWEM
import WEM
import json

"""Create velocity model"""
path = '/home/user/project/projects/constWrongBg'
slow = genericIO.defaultIO.getVector(path + '/temp.H')

"""Read in wavelet"""
wave = genericIO.defaultIO.getVector(path+'/Wav/constWav.H')

"""Read in parameters"""
with open(path + '/Par/geom.json') as f:
    par = json.load(f)
par['osx'] = 5000
par['ns'] = 1
parObj = genericIO.pythonParams(par).getCpp()

"""Create data"""
nt = wave.getHyper().getAxis(1).n
dt = wave.getHyper().getAxis(1).d
ot = 0.
oshot = parObj.getFloat('osx')
dshot = parObj.getFloat('dsx')
nshot = parObj.getInt('ns')
n1 = slow.getHyper().getAxis(1).n
o1 = slow.getHyper().getAxis(1).o
d1 = slow.getHyper().getAxis(1).d
n2 = slow.getHyper().getAxis(2).n
o2 = slow.getHyper().getAxis(2).o
d2 = slow.getHyper().getAxis(2).d

data = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,n1,n2,nshot],
													os=[ot,o1,o2,oshot],
													ds=[dt,d1,d2,dshot]))

wem = newExtWEM.extWEM(slow,data,parObj,wave)
wem.wavefield(slow,data)

genericIO.defaultIO.writeVector('ext_wavefield.H',data)
