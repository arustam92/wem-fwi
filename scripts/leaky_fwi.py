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

slow2 = genericIO.defaultIO.getVector("Vel/constSlowBg1800_ext.H")
slow = slow2.clone()
alpha = [0.5, 0.8, 0.999, 1]
# alpha_final = [0.1,0.5,0.7,1]
eps = [1e-5,0.01,0.01,0.01]
niter = 20
subiter = 50

for i in range(len(alpha)):
    ######## prepare the input #############
    D = "./"
    with open('Par/geom.json') as f:
        par = json.load(f)
    parObj = genericIO.pythonParams(par).getCpp()
    wave = genericIO.defaultIO.getVector(D+"Wav/constWav.H")
    data = genericIO.defaultIO.getVector(D+"Dat/constData_tap.H")
    ######## prepare the input #############

    ######## prepare the precondioned variable #############
    slow.zero()
    int = Operator.LeakyIntegration(slow,slow2,alpha=alpha[i])
    LinStop  = Stopper.BasicStopper(niter=subiter)
    CGsolver = LinearSolver.LCGsolver(LinStop)
    L2Prob = Prblm.ProblemL2Linear(slow,slow2,int)
    CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
    			flush_memory=True,iter_buffer_size=1,save_model=False )
    CGsolver.run(L2Prob,verbose=True)
    ######## prepare the precondioned variable #############

    ######## prepare the regularization term #############
    der = Operator.Derivative(slow,slow)
    reg = Op.ChainOperator(int,der)
    # derOp = Op.NonLinearOperator(der,der)
    regOp = Op.NonLinearOperator(reg,reg)
    ######## prepare the regularization term #############

    ######## prepare the data fitting term #############
    smoothOp = Operator.Taper(slow,slow,tap=[10],axes=[2],zero=[10])
    bornOp = newExtWEM.extBorn(slow,data,wave,parObj)
    bornOpSm = Op.ChainOperator(smoothOp,bornOp)
    # bornOpSm = Op.ChainOperator(realSmooth,bornOpSm)
    wemOp = newExtWEM.extWEM(slow,data,parObj,wave)
    nlOp = Op.NonLinearOperator(wemOp,bornOpSm,set_background_func=bornOp.setBgSlow)
    intOp = Op.NonLinearOperator(int,int)
    nlOp = Op.CombNonlinearOp(intOp,nlOp)
    ######## prepare the data fitting term ############

    ######## prepare the FWI problem #############
    L2ProbReg = Prblm.ProblemL2NonLinearReg(slow,data,nlOp,eps[i],reg_op=regOp)
    ######## find eps #############
    L2ProbReg.epsilon = eps[i] * L2ProbReg.estimate_epsilon(verbose=True)
    print(L2ProbReg.epsilon)
    ######## find eps #############
    Stop  = Stopper.BasicStopper(niter=niter)
    stepper = Stepper.CvSrchStep(alpha_min=1e-7)
    LBFGSsolver = Solver.LBFGSsolver(Stop,stepper=stepper,m_steps=50,logger=logger("Log/1800_leaky_%s_alpha_%s.txt" % (eps[i],alpha[i])))
    LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
    			flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/1800_leaky_%s_%s" % (eps[i],alpha[i]) )
    LBFGSsolver.run(L2ProbReg,verbose=True,restart=False)

    slow = genericIO.defaultIO.getVector("Inv/1800_leaky_%s_%s_inv_mod.H" % (eps[i],alpha[i]) )
    int.forward(False,slow,slow2)
