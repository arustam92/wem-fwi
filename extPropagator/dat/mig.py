import SepVector
import Hypercube
import genericIO
import WEM

D = "./Dat"
p = "/net/server2/homes/sep/arustam/projects/testWEM"
par = ["","wave=wave.H","par=geom.par"]

io=genericIO.pyGenericIO.ioModes(par)
parObj=io.getDefaultIO().getParamObj()

# sigs_path = "gauss.H"
slow = genericIO.defaultIO.getVector("bg.H")
wave = genericIO.defaultIO.getVector("wave.H")

# dimensions 
oshot = parObj.getFloat("osx")
dshot = parObj.getFloat("dsx")
nshot = parObj.getInt("ns") 
orec = parObj.getFloat("orx")
drec = parObj.getFloat("drx")
nrec = parObj.getInt("nr")
nt = wave.getHyper().getAxis(1).n
dt = wave.getHyper().getAxis(1).d
ot = 0.

image = slow.clone()
image.zero()
# write the output

data = genericIO.defaultIO.getVector(D+"/data.H")
# create a WEM propagator
bornOp = WEM.Born(slow,data,wave,parObj)
bornOp.adjoint(False,image,data)

# for i in range (nshots):
# 	data = genericIO.defaultIO.getVector(D+"/syn%s.H"%i)
# 	print("Migrating shot %d" % i)
# 	bornOp.adjoint(True,image,data)
# 	bornOp.nextSrc()

genericIO.defaultIO.writeVector(D+"/image.H",image)