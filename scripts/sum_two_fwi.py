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
slow = genericIO.defaultIO.getVector("Vel/constSlowBg1800_ext.H")
slow2 = genericIO.defaultIO.getVector("Vel/constSlowBg1800.H")
alpha = 0.
niter = 50
subiter = 30
axes = slow.getHyper().axes
# smooth = [11,1,7,7]
points = []
cases = [20]
for i in cases:
    points.append(np.linspace(1,25,i).tolist())

for i in range(len(cases)):

    ######## prepare the precondioned variable #############
    # prep = Operator.LeakyIntegrationInverse(slow,slow,alpha=alpha)
    # prep.forward(slow2,slow)
    ######## prepare the precondioned variable #############

    D = "./"
    with open('Par/geom.json') as f:
        par = json.load(f)
    parObj = genericIO.pythonParams(par).getCpp()
    wave = genericIO.defaultIO.getVector(D+"Wav/constWav.H")
    data = genericIO.defaultIO.getVector(D+"Dat/constData_tap.H")

    # REGULARIZATION
    # rect=(smooth,[1,1,1,1])
    # repeat=1
    # realSmooth = Operator.Smooth(slow,slow,rect,repeat,mod='all')
    smoothOp = Operator.Taper(slow,slow,tap=[10],axes=[2])
    # int = Operator.LinearInterpolation(slow,slow,points=points[i])
    der = Operator.Derivative(slow,slow)
    # reg = Op.ChainOperator(int,der)
    # derOp = Op.NonLinearOperator(der,der)
    # intOp = Op.NonLinearOperator(int,int)
    derOp = Op.NonLinearOperator(der,der)

    bornOp = newExtWEM.extBorn(slow,data,wave,parObj)
    bornOpSm = Op.ChainOperator(smoothOp,bornOp)
    # bornOpSm = Op.ChainOperator(realSmooth,bornOpSm)
    wemOp = newExtWEM.extWEM(slow,data,parObj,wave)
    nlOp1 = Op.NonLinearOperator(wemOp,bornOpSm,set_background_func=bornOp.setBgSlow)
    # nlOp = Op.CombNonlinearOp(intOp,nlOp)

    stack = Operator.Stack(slow,slow2)
    stackOp = Op.NonLinearOperator(stack,stack)
    wemOp = WEM.WEM(slow2,data,parObj,wave)
    bornOp = WEM.Born(slow2,data,wave,parObj)
    # bornOpSm = Op.ChainOperator(smoothOp,bornOp)
    nlOp2 = Op.NonLinearOperator(wemOp,bornOp,set_background_func=bornOp.setBgSlow)
    nlOp2 = Op.CombNonlinearOp(stackOp,nlOp2)

    nlOp = Op.sumNlOperator(nlOp2,nlOp1)

    #
    #Create L2-norm linear problem
    #Create stopper
    Stop  = Stopper.BasicStopper(niter=niter)
    stepper = Stepper.CvSrchStep(alpha_min=1e-7)
    # stepper.alpha =
    LBFGSsolver = Solver.LBFGSsolver(Stop,stepper=stepper,m_steps=50,logger=logger("sum_ext1800_%s_alpha_%s.txt" % (eps,alpha)))

    L2ProbReg = Prblm.ProblemL2NonLinearReg(slow,data,nlOp,eps,reg_op=derOp)
    #Running the solver
    # eps[i] = eps[i] * L2ProbReg.estimate_epsilon(verbose=True)
    # L2ProbReg.epsilon = eps[i]
    LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
    			flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/sum_%s" % eps )
    LBFGSsolver.run(L2ProbReg,verbose=True,restart=False)
