import sys
sys.path.insert(0,"/sep/arustam/package/python-solver/lib/python")

import SepVector
import Hypercube
import genericIO
import WEM

import pyVector as Vec
import pyOperator as Op
import pyLCGsolver as LCG
import pySymLCGsolver as SymLCGsolver
import pyProblem as Prblm
import pyStopperBase as Stopper
import sep_util as sep

D = "./Dat"
p = "/net/server2/homes/sep/arust'am/projects/testWEM"
par = ["","wave=wave.H","par=geom.par"]
io=genericIO.pyGenericIO.ioModes(par)
parObj=io.getDefaultIO().getParamObj()

# sigs_path = "gauss.H"
slow = genericIO.defaultIO.getVector("bg.H")
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
data = genericIO.defaultIO.getVector(D+"/data.H")
bornOp = WEM.Born(slow,data,wave,parObj)

model = slow.clone()
model.zero()
#Create L2-norm linear problem
L2Prob = Prblm.ProblemL2Linear(model,data,bornOp)
#Create stopper
niter = 10
Stop  = Stopper.BasicStopper(niter=niter)
#Create solver
LCGsolver = LCG.LCGsolver(Stop)
#Running the solver
LCGsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,
					iter_buffer_size=100,save_model=True,prefix="Inv/") 
LCGsolver.run(L2Prob,verbose=True,restart=False)

