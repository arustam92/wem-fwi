import SepVector
import Hypercube
import genericIO
import WEM
import Operator

par = {
    "ns" : 1,
    "osx" : 650,
    "osz" : 708,
    "dsx" : 50,

    "nr" : 445,
    "orx" : 650,
    "orz" : 708,
    "drx" : 12,

    "prop" : 'ssf',
    "nref" : 10,
    "tap" : 50,
    "fmin" : 1,
    "fmax" : 15,
    "nfreq" : 1,
    "ng" : 1,

    "illum" : 1
}
parObj = genericIO.pythonParams(par).getCpp()


slow = genericIO.defaultIO.getVector("Vel/MarmSlow.H")
wave = genericIO.defaultIO.getVector("Wav/c_wave_3-40.H")

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
dataF = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,nrec,nshot],
													os=[ot,orec,oshot],
													ds=[dt,drec,dshot]), storage="dataComplex")

fftOp = Operator.FFT(data,dataF,1,1)
propWEM = WEM.WEM(slow,dataF,parObj,wave)

# write the output

propWEM.forward(False,slow,dataF)

genericIO.defaultIO.writeVector("test.H",dataF)
