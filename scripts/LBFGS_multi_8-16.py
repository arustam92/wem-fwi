import numpy as np
#sys.path.insert(0,"/sep/arustam/package/python-solver/lib/python")
# sys.path.insert(0,"/sep/arustam/python-solver/GenericSolver/python")

import SepVector
import Hypercube
import genericIO
import WEM
import Operator

import pyVector as Vec
import pyOperator as Op
import pyLBFGSsolver as LBFGS
import pyStepperCvSrch as CVstepper
#import pySymLCGsolver as SymLCGsolver
import pyProblem as Prblm
import pyStopperBase as Stopper
from sys_util import logger


D = "./Dat"
par = ["","par=geom.par"]
io=genericIO.pyGenericIO.ioModes(par)
parObj=io.getDefaultIO().getParamObj()

slow = genericIO.defaultIO.getVector("Multi/ricker15_1-8Hz_inv_mod.H")
wave = genericIO.defaultIO.getVector("Wav/c_ricker15_8-16Hz.H",storage="dataComplex")
data = genericIO.defaultIO.getVector("Dat/c_dataGauss_ricker15_8-16Hz.H",storage="dataComplex")

# create a WEM propagator
rect=[1,7]
repeat=3
smoothOp = WEM.Smooth(slow,slow,rect,repeat)

bornOp = WEM.Born(slow,data,wave,parObj)
wemOp = WEM.WEM(slow,data,parObj,wave)
bornOpSmooth = Op.ChainOperator(smoothOp,bornOp)
nlOp = Op.NonLinearOperator(wemOp,bornOpSmooth,set_background_func=bornOp.setBgSlow)

#Create L2-norm linear problem
niter = 100
L2Prob = Prblm.ProblemL2NonLinear(slow,data,nlOp)
#Create stopper
Stop  = Stopper.BasicStopper(niter=niter)
#Create solver
cvStep = CVstepper.CvSrchStep()
LBFGSsolver = LBFGS.LBFGSsolver(Stop,stepper=cvStep,logger=logger("LBFGS_multi_freq.txt"))
#Running the solver
LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,
			flush_memory=True,iter_buffer_size=5,save_model=True,prefix="Multi/ricker15_mix_8-16Hz")
LBFGSsolver.run(L2Prob,verbose=True,restart=False)
