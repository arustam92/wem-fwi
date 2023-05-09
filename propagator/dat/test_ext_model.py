import Conformal
import genericIO
import SepVector
import Hypercube
import numpy as np
import newExtWEM
import WEM
import json

"""Read in parameters"""
with open('geometry.json') as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()

"""Create velocity model"""
n1 = 500
n2 = 500
o1 = 0
omax1 = n1
o2 = 0
omax2 = n2
d1= float(omax1-o1)/n1
d2 = float(omax2-o2)/n2

slow = SepVector.getSepVector(Hypercube.hypercube(ns=[n1,n2],ds=[d1,d2],os=[o1,o2]),storage='dataComplex')
slowNd = slow.getNdArray()
slowNd[:] = 1 / 1000
slowNd[int(n2/2)::,:] = 1 / 2000

tmax = 0.2
fmax = parObj.getFloat('fmax')
model = newExtWEM.createExtModel(slow,tmax,fmax)
print(model.getNdArray().shape)

"""Read in wavelet"""
wave = genericIO.defaultIO.getVector('wave.H')

"""Create data"""
nt = wave.getHyper().getAxis(1).n
dt = wave.getHyper().getAxis(1).d
ot = 0.
oshot = parObj.getFloat('osx')
dshot = parObj.getFloat('dsx')
nshot = parObj.getInt('ns')
orec = parObj.getFloat('orx')
drec = parObj.getFloat('drx')
nrec = parObj.getInt('nr')

data = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,nrec,nshot],
													os=[ot,orec,oshot],
													ds=[dt,drec,dshot]))

wem = newExtWEM.extWEM(model,data,parObj,wave)
wem.forward(False,model,data)

image = model.clone()
born = newExtWEM.extBorn(model,data,wave,parObj)
born.adjoint(False,image,data)

genericIO.defaultIO.writeVector('ext_data.H',data)
image.writeVec('ext_image.H')
