import sys
sys.path.insert(0,"/sep/arustam/package/python-solver/lib/python")

import SepVector
import Hypercube
import genericIO
import WEM

import pyVector as Vec
import pyOperator as Op
import pyNLCGsolver as NLCG
import pySymLCGsolver as SymLCGsolver
import pyProblem as Prblm
import pyStopperBase as Stopper
from sys_util import logger


D = "./Dat"
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

data = genericIO.defaultIO.getVector(D+"/dataGauss1.H")

# create a WEM propagator
bornOp = WEM.Born(slow,data,wave,parObj)
wemOp = WEM.WEM(slow,data,parObj,wave)
nlOp = Op.NonLinearOperator(wemOp,bornOp,set_background_func=bornOp.setBgSlow)

# model = slow.clone()
# model.zero()
#Create L2-norm linear problem
L2Prob = Prblm.ProblemL2NonLinear(slow,data,nlOp)
#Create stopper
niter = 2
Stop  = Stopper.BasicStopper(niter=niter)
#Create solver
NLCGsolver = NLCG.NLCGsolver(Stop,logger=logger("gauss1_NLCG_log.txt"))
#Running the solver
NLCGsolver.setDefaults(iter_sampling=1,save_obj=True,iter_buffer_size=5,save_model=True,save_res=True,prefix="NL_Inv/") 
NLCGsolver.run(L2Prob,verbose=True,restart=False)

