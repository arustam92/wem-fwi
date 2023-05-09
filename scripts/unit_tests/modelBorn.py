import SepVector
import Hypercube
import genericIO
import ExtWEM
import json

par = {
    "ns" : 20,
    "osx" : 0,
    "osz" : 0,
    "dsx" : 50,

    "nr" : 100,
    "orx" : 0,
    "orz" : 0,
    "drx" : 10,

    "prop" : 'ssf',
    "nref" : 1,
    "tap" : 10,
    "fmin" : 1,
    "fmax" : 15,
    "nfreq" : 1,
    "ng" : 1,

    "illum" : 1
}
parObj = genericIO.pythonParams(par).getCpp()

slow = genericIO.defaultIO.getVector("extBg_wrong.H")
dslow = genericIO.defaultIO.getVector("mig_wrong.H")
wave = genericIO.defaultIO.getVector("c_wave_1-15.H")

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
													ds=[dt,drec,dshot]),storage='dataComplex')

propWEM = ExtWEM.ExtBorn(slow,data,wave,parObj)

# write the output
propWEM.forward(False,dslow,data)

# dataWin = genericIO.regFile.writeWindow(data,n1=int(nt/2))
genericIO.defaultIO.writeVector("unit_test_spike_data_Born.H", data)
