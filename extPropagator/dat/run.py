import SepVector
import Hypercube
import genericIO
import WEM

p = "/net/server2/homes/sep/arust'am/projects/testWEM"
par = ["","wave=wave.H","par=geom.par"]
io=genericIO.pyGenericIO.ioModes(par)
parObj=io.getDefaultIO().getParamObj()

slow = genericIO.defaultIO.getVector("gauss.H")
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

propWEM = WEM.WEM(slow,data,parObj,wave)

# write the output
D = "./Dat"

propWEM.forward(False,slow,data)

# dataWin = genericIO.regFile.writeWindow(data,n1=int(nt/2))
genericIO.defaultIO.writeVector(D+"/data.H",data)