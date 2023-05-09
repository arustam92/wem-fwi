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

snap = 0
for par in sys.argv[1:]:
    if ('vel' in par):
        vel = par.split('vel=')[1]
        slow = genericIO.defaultIO.getVector(vel)
    if ('wave' in par):
        wave = par.split('wave=')[1]
        wave = genericIO.defaultIO.getVector(wave)
    if ('geom' in par):
        geom = par.split('geom=')[1]
        with open(geom) as f:
            p = json.load(f)
        parObj = genericIO.pythonParams(p)
    if ('out' in par):
        out = par.split('out=')[1]
    if ('refl' in par):
        refl = par.split('refl=')[1]
        refl = genericIO.defaultIO.getVector(refl)
    if ('snap' in par):
        snap = par.split('snap=')[1]

# modeling
oshot = float(parObj.pars["osx"])
dshot = float(parObj.pars["dsx"])
nshot = int(parObj.pars["ns"])
orec = float(parObj.pars["orx"])
drec = float(parObj.pars["drx"])
nrec = int(parObj.pars["nr"])
nt = wave.getHyper().getAxis(1).n
dt = wave.getHyper().getAxis(1).d
ot = 0.

# create a WEM propagator
if snap:
    ax = slow.getHyper().axes
    data = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,ax[0].n,ax[1].n,nshot],
													os=[ot,ax[0].o,ax[1].o,oshot],
													ds=[dt,ax[0].d,ax[1].d,dshot]))
else:
    data = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,nrec,nshot],
                                                        os=[ot,orec,oshot],
                                                        ds=[dt,drec,dshot]))
if len(slow.getNdArray().shape) < 3:
    propWEM = WEM.WEM(slow,data,wave,parObj)
else:
    propWEM = newExtWEM.extWEM(slow,data,wave,parObj)

######## test #############
hyp = Operator.HyperbolicPenalty(slow,slow,l=1,tau=p["tau"])
propWEM = Op.ChainOperator(hyp,propWEM)
######## test #############

if snap:
    propWEM.wavefield(slow,data)
else:
    propWEM.forward(False,slow,data)

if ('refl' in par):
    bornRefl = WEM.BornRefl(refl_pad,slow,data,wave,parObj)
    bornRefl.forward(True,refl_pad,data)

# dataWin = genericIO.regFile.writeWindow(data,n1=int(nt/2))
data.writeVec(out)
