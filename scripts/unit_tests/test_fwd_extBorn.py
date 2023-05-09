import SepVector
import Hypercube
import genericIO
import newExtWEM
import json

with open('geometry.json') as f:
    geom = json.load(f)
parObj = genericIO.pythonParams(geom).getCpp()

slow = genericIO.defaultIO.getVector("ext.H")
wave = genericIO.defaultIO.getVector("wave.H")

# modeling
oshot = parObj.getFloat("osx")
dshot = parObj.getFloat("dsx")
nshot = parObj.getInt("ns")
orec = parObj.getFloat("orx")
drec = parObj.getFloat("drx")
nrec = parObj.getInt("nr")
nt = wave.getHyper().getAxis(1).n
dt = wave.getHyper().getAxis(1).d
ot = 0.

# create a WEM propagator
data = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,nrec,nshot],
													os=[ot,orec,oshot],
													ds=[dt,drec,dshot]))
spike = slow.clone()
spike.zero()
spikeNd = spike.getNdArray()
spikeNd[:, int(spikeNd.shape[1]/2), int(spikeNd.shape[2]/2)] = 1

propWEM = newExtWEM.extBorn(slow,data,wave,parObj)

propWEM.forward(False,spike,data)
genericIO.defaultIO.writeVector("test_fwd_extBorn.H", data)
propWEM.adjoint(False,spike,data)
genericIO.defaultIO.writeVector("test_adj_extBorn.H", spike)
