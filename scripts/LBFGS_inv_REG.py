import SepVector
import Hypercube
import genericIO
import WEM
import ExtWEM
import Operator
import json

import pyOperator as Op
import pyNonLinearSolver as Solver
import pyProblem as Prblm
import pyStopper as Stopper
from sys_util import logger


D = "./"
with open('Par/geom.json') as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()

slow = genericIO.defaultIO.getVector("Vel/constSlowBg2000_ext.H")

wave = genericIO.defaultIO.getVector(D+"Wav/c_constWav_3-25.H")
data = genericIO.defaultIO.getVector(D+"Dat/c_constData_3-25.H")

# create a WEM propagator
rect=[5,5]
repeat=3
smoothOp = Operator.Smooth(slow,slow,rect,repeat)
bornOp = ExtWEM.ExtBorn(slow,data,wave,parObj)
# bornOp.setTomoMod()
bornOpSm = Op.ChainOperator(smoothOp,bornOp)
wemOp = ExtWEM.ExtWEM(slow,data,parObj,wave)
nlOp = Op.NonLinearOperator(wemOp,bornOpSm,set_background_func=bornOp.setBgSlow)

der = Operator.Derivative(slow,slow)
derOp = Op.NonLinearOperator(der,der)
#
#Create L2-norm linear problem
eps = 0.5
#Create stopper
niter = 10
Stop  = Stopper.BasicStopper(niter=niter)
#Create solver
LBFGSsolver = Solver.LBFGSsolver(Stop,logger=logger("LBFGS_egVzExt"))

L2ProbReg = Prblm.ProblemL2NonLinearReg(slow,data,nlOp,eps,reg_op=derOp)
#Running the solver
LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
			flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/const2_reg_eps%s" % eps)
LBFGSsolver.run(L2ProbReg,verbose=True,restart=False)
