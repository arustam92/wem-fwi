import SepVector
import Hypercube
import genericIO
import WEM

import numpy as np
import matplotlib.pyplot as plt

D = "./Dat"
p = "/net/server2/homes/sep/arustam/projects/testWEM"
par = ["","wave=wave.H","par=geom.par"]

io=genericIO.pyGenericIO.ioModes(par)
parObj=io.getDefaultIO().getParamObj()

# sigs_path = "gauss.H"
slow = genericIO.defaultIO.getVector("slow.H")
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


dataWEM = SepVector.getSepVector(Hypercube.hypercube(ns=[nt,nrec,nshot],
													os=[ot,orec,oshot],
													ds=[dt,drec,dshot]))

# create a WEM propagator
wemOp = WEM.WEM(slow,dataWEM,parObj,wave)
bornOp = WEM.Born(slow,dataWEM,wave,parObj)
dataBorn = dataWEM.clone()
data0 = dataWEM.clone()
wemOp.forward(False,slow,data0)

############## perturbation #################
loc2=int(slow.getHyper().getAxis(2).n/2)
loc1=int(slow.getHyper().getAxis(1).n/2)
dvel0 = slow.clone()
dvel0.zero()
dvel0N = dvel0.getNdArray()
s0 = slow.getNdArray()[loc2][loc1]
perc = 0.01
# dvel0.scale(smax*perc)
dvel0N[loc2][loc1] = perc*s0
print (dvel0.dot(dvel0))

dvel = dvel0.clone()

res = np.zeros(5)

for i in range(1,5):

	bornOp.forward(False,dvel,dataBorn)	#Bdx

	dataBorn.scaleAdd(data0,1.,1.)

	slow.scaleAdd(dvel,1.,1.)	#s + ds
	wemOp.forward(False,slow,dataWEM)	#f(s+ds)

	dataBorn.scaleAdd(dataWEM,1.,-1.)

	res[i] = dataBorn.dot(dataBorn)
	print(res[i])

	dvel.scaleAdd(dvel0,1.,1.);


res = res/res[1]
x = np.linspace(0,4,5)
x2 = x**2
fig = plt.figure()
plt.plot(x,x2,'r',label='Theoretical')
plt.plot(x,res,'g',label='Computed')
plt.xlabel('Percent of perturbation')
plt.ylabel('Absolute error')
plt.legend()
plt.show()

fig.savefig("linTest.pdf",bbox_inches='tight')
