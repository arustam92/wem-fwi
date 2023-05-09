import SepVector
import Hypercube
import genericIO
import WEM
import newExtWEM
import Operator
import json
import numpy as np

import pyOperator as Op
import pyNonLinearSolver as Solver
import pyLinearSolver as LinearSolver
import pyProblem as Prblm
import pyStopper as Stopper
import pyStepper as Stepper
from sys_util import logger

eps = 0.1
slow2 = genericIO.defaultIO.getVector("Vel/constSlowBg1800_ext.H")
slow = slow2.clone()
niter = 10

######## prepare the input #############
D = "./"
with open('Par/geom.json') as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()
wave = genericIO.defaultIO.getVector(D+"Wav/constWav.H")
data = genericIO.defaultIO.getVector(D+"Dat/constData_tap.H")
######## prepare the input #############

######## prepare the regularization term #############
der = Operator.Derivative(slow,slow)
regOp = Op.NonLinearOperator(der,der)
######## prepare the regularization term #############

######## prepare the data fitting term #############
smoothOp = Operator.Taper(slow,slow,tap=[10],axes=[2])
bornOp = newExtWEM.extBorn(slow,data,wave,parObj)
bornOpSm = Op.ChainOperator(smoothOp,bornOp)
wemOp = newExtWEM.extWEM(slow,data,parObj,wave)
nlOp = Op.NonLinearOperator(wemOp,bornOpSm,set_background_func=bornOp.setBgSlow)

######## prepare the data fitting term #############
######## prepare the FWI problem #############
L2ProbReg = Prblm.ProblemL2NonLinearReg(slow,data,nlOp,eps,reg_op=regOp)
Stop  = Stopper.BasicStopper(niter=niter)
stepper = Stepper.CvSrchStep(alpha_min=1e-7)
LBFGSsolver = Solver.LBFGSsolver(Stop,stepper=stepper,m_steps=50,logger=logger("Log/test" ))
LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
			flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/test")
LBFGSsolver.run(L2ProbReg,verbose=True,restart=False)
