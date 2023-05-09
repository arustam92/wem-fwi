import numpy as np
import SepVector
import Hypercube
import genericIO
import WEM
import Operator

import pyVector as Vec
import pyOperator as Op
import pyNonLinearSolver as Solver
#import pySymLCGsolver as SymLCGsolver
import pyProblem as Prblm
import pyStopper as Stopper
from sys_util import logger
import json

D = "./Dat"
with open('Par/geom_refl.json') as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()

with open('Par/geom.json') as f:
    par = json.load(f)
parObjTomo = genericIO.pythonParams(par).getCpp()

slow = genericIO.defaultIO.getVector("Vel/bodyBgSlow.H")
wave = genericIO.defaultIO.getVector("Wav/medWave.H")
dataTomo = genericIO.defaultIO.getVector("Dat/body_datWEM.H")
data = genericIO.defaultIO.getVector("Dat/body_datWEM_refl.H")

# create a WEM propagator
rect=([1,1,7],[1,1,7])
repeat=1
# smoothOp = Operator.Smooth(slow,slow,rect,repeat,mod='top')
smoothOp = Operator.Taper(slow,slow,tap=[20],axes=[1],zero=[20])

bornOp = WEM.Born(slow,data,wave,parObj)
wemOp = WEM.WEM(slow,data,parObj,wave)
bornOpSmooth = Op.ChainOperator(smoothOp,bornOp)
nlOp = Op.NonLinearOperator(wemOp,bornOpSmooth,set_background_func=bornOp.setBgSlow)

bornOpTomo = WEM.Born(slow,data,wave,parObjTomo)
wOpTomo = WEM.WEM(slow,data,parObjTomo,wave)
bornOpSmoothTomo = Op.ChainOperator(smoothOp,bornOpTomo)
wemOpTomo = Op.NonLinearOperator(wOpTomo,bornOpSmoothTomo,set_background_func=bornOpTomo.setBgSlow)

#Create L2-norm linear problem
niter = 50
L2Prob = Prblm.ProblemL2NonLinearReg(slow,data,nlOp,0.1,reg_op=wemOpTomo,prior_model=dataTomo)
#Create stopper
Stop  = Stopper.BasicStopper(niter=niter)
#Create solver
LBFGSsolver = Solver.LBFGSsolver(Stop,logger=logger("fwi.txt"))
#Running the solver
LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
			flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/body_TR_0.1")
LBFGSsolver.run(L2Prob,verbose=True,restart=False)
