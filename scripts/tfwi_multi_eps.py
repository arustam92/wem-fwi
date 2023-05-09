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

eps = 1
# slow = genericIO.defaultIO.getVector("Vel/constSlowBg1800.H")
alpha = 0.
niter = 5000

######## prepare the input #############
D = "./"
with open('Par/geom.json') as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()
wave = genericIO.defaultIO.getVector(D+"Wav/constWav.H")
data = genericIO.defaultIO.getVector(D+"Dat/constData_tap.H")
######## prepare the input #############
model = genericIO.defaultIO.getVector("Vel/constSlowBg1800_ext.H")
# model = newExtWEM.createExtModel(slow,0.5,25)


######## prepare the precondioned variable #############
# slow2 = SepVector.getSepVector(Hypercube.hypercube(ns=[axes[0].n,axes[1].n,axes[2].n,len(points[i])],
# 													os=[axes[0].o,axes[1].o,axes[2].o,0],
# 													ds=[axes[0].d,axes[1].d,axes[2].d,1]),storage='dataComplex')

# LinStop  = Stopper.BasicStopper(niter=subiter)
# CGsolver = LinearSolver.LCGsolver(LinStop)
# L2Prob = Prblm.ProblemL2Linear(slow2,slow,int)
# CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
# 			flush_memory=True,iter_buffer_size=1,save_model=False )
# CGsolver.run(L2Prob,verbose=True)
######## prepare the precondioned variable #############

######## prepare the regularization term #############
der = Operator.Derivative(model,model)
# derOp = Op.NonLinearOperator(der,der)
regOp = Op.NonLinearOperator(der,der)
######## prepare the regularization term #############

######## prepare the data fitting term #############
smoothOp = Operator.Taper(model,model,tap=[5],axes=[1],zero=[5])
# rect=([3,1,3],[3,1,3])
# repeat=1
# smoothOp = Operator.Smooth(model,model,rect,repeat,mod='all')

bornOp = newExtWEM.extBorn(model,data,wave,parObj)
bornOpSm = Op.ChainOperator(smoothOp,bornOp)
# bornOpSm = Op.ChainOperator(realSmooth,bornOpSm)
wemOp = newExtWEM.extWEM(model,data,parObj,wave)
nlOp = Op.NonLinearOperator(wemOp,bornOpSm,set_background_func=bornOp.setBgSlow)
# intOp = Op.NonLinearOperator(int,int)
# nlOp = Op.CombNonlinearOp(intOp,nlOp)
######## prepare the data fitting term #############

######## prepare the FWI problem #############
L2ProbReg = Prblm.ProblemL2NonLinearReg(model,data,nlOp,eps,reg_op=regOp)
Stop  = Stopper.BasicStopper(niter=niter)
stepper = Stepper.CvSrchStep(alpha_min=1e-7)
LBFGSsolver = Solver.LBFGSsolver(Stop,stepper=stepper,m_steps=50,logger=logger("Log/1800_eps_%s.txt" % (eps)))
LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
            flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/1800_full_eps_%s" % (eps) )
LBFGSsolver.run(L2ProbReg,verbose=True,restart=False)

# slow2 = genericIO.defaultIO.getVector("Inv/1800_lin_%s_%s_inv_mod.H" % (eps,len(points[i])) )
# int.forward(False,slow2,slow)
# del slow2
