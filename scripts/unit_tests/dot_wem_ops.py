import SepVector
import Hypercube
import genericIO
import pyOperator as Op
import wem_ops
import WEM
import Operator
import json
import os
import sys
import numpy as np

ns = [100,100]
ds = [10,10]
os = [0,0]
vel = 1000
freq = 10
loc = np.random.randint(0,ns[0]-1,size=10).tolist()
s0 = 1/vel * .1 + 1j*1
depth = 50

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
    "tap" : 0,
    "fmin" : 1,
    "fmax" : 15,
    "ng" : 0,
    "ntaylor" : 1,

    "onepass" : 0,
}
parObj = genericIO.pythonParams(par).getCpp()

slow = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
slowNd = slow.getNdArray()
slowNd[:] = 1/vel
# slowNd.real[:] += 1/100 * .01* np.random.rand(ns[0],ns[1])
# slowNd.imag[:] += 1/100 * .01* np.random.rand(ns[0],ns[1])
slowNd[int(ns[1]/2),:] *= .5
#
# slow.rand()

print(np.amax(1/slowNd),np.amin(1/slowNd))

ref = wem_ops.RefSampler(slow,par["nref"])

model = SepVector.getSepVector(Hypercube.hypercube(ns=[ns[0]],ds=[ds[0]],os=[os[0]]),storage='dataComplex')
data = model.clone()

# print("SplitStep: ")
# ss = wem_ops.SplitStep(model,data,slow,freq,loc,s0,depth)
# ss.dotTest(verbose=True)

print("Phshift: ")
ph = wem_ops.Phshift(model,data,slow.getHyper(),ds[1],freq,s0,0)
ph.dotTest(verbose=True)

print("Scatter: ")
scat = wem_ops.Scatter(model,data,slow,par["ntaylor"],freq,depth,parObj)
scat.dotTest(verbose=True)
#
print("SSF: ")
scat = wem_ops.SSF(model,data,slow,freq,depth,parObj,ref)
scat.dotTest(verbose=True)

print("IC: ")
ic = wem_ops.IC(model,data,slow,depth)
ic.dotTest(verbose=True)

print("Reflect: ")
model = slow.clone()
data = slow.clone()
refl = wem_ops.Reflect(model,data,slow)
refl.dotTest(verbose=True)

print("Down: ")
model = slow.clone()
data = slow.clone()
down = wem_ops.Down(model,data,slow,parObj,ref,freq)
down.dotTest(verbose=True)

print("Up: ")
model = slow.clone()
data = slow.clone()
up = wem_ops.Up(model,data,slow,parObj,ref,freq)
up.dotTest(verbose=True)

print("LinDown: ")
model = slow.clone()
data = slow.clone()
bg = slow.clone()
bgNd = bg.getNdArray()
bgNd[:] = 1
lindown = wem_ops.LinDown(model,data,slow,parObj,bg,down,freq)
lindown.dotTest(verbose=True)

print("dReflect: ")
model = slow.clone()
data = slow.clone()
drefl = wem_ops.dReflect(model,data,slow)
drefl.dotTest(verbose=True)
