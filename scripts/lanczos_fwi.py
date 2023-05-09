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

slow = genericIO.defaultIO.getVector("Vel/constSlowBg1800_ext.H")
ax = slow.getHyper().axes
ns = [[51,51,51]]   # have to be odd numbers
# alpha_final = [0.1,0.5,0.7,1]
eps = [0.1,0.1,0.1]
niter = 200
subiter = 15
a = [[11,11,11],[11,11,11],[3,3,3]]

for i in range(len(ns)):
    ######## prepare the input #############
    D = "./"
    with open('Par/geom.json') as f:
        par = json.load(f)
    parObj = genericIO.pythonParams(par).getCpp()
    wave = genericIO.defaultIO.getVector(D+"Wav/constWav.H")
    data = genericIO.defaultIO.getVector(D+"Dat/constData.H")

    ds = []
    for j in range(3):
        ds.append((ax[j].n-1)*ax[j].d / (ns[i][j] - 1))

    sub_slow = SepVector.getSepVector(Hypercube.hypercube(ns=[ns[i][0],ns[i][1],ns[i][2]],
                                            ds=[ds[0],ds[1],ds[2]],
                                            os=[ax[0].o,ax[1].o,ax[2].o]),storage='dataComplex')
    ######## prepare the input #############

    ######## prepare the precondioned variable #############
    sub_slow.zero()
    int = Operator.LanczosInterpolation3D(sub_slow,slow,a=a[i])
    LinStop  = Stopper.BasicStopper(niter=subiter)
    CGsolver = LinearSolver.LCGsolver(LinStop)
    L2Prob = Prblm.ProblemL2Linear(sub_slow,slow,int)
    CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
    			flush_memory=True,iter_buffer_size=1,save_model=False,prefix="temp")
    CGsolver.run(L2Prob,verbose=True)

    ######## prepare the precondioned variable #############

    ####### prepare the regularization term #############
    rect = [
        [1,5,5],[1,5,5]
    ]
    diff = Operator.SmoothDiff(slow,slow,rect)
    der = Operator.Derivative(slow,slow)
    reg = Op.ChainOperator(int,diff)
    reg = Op.ChainOperator(reg,der)
    regOp = Op.NonLinearOperator(reg,reg)
    ######## prepare the regularization term #############

    ######## prepare the data fitting term #############
    smoothOp = Operator.Taper(slow,slow,tap=[5],axes=[1],zero=[5])
    bornOp = newExtWEM.extBorn(slow,data,wave,parObj)
    bornOpSm = Op.ChainOperator(smoothOp,bornOp)
    wemOp = newExtWEM.extWEM(slow,data,parObj,wave)
    nlOp = Op.NonLinearOperator(wemOp,bornOpSm,set_background_func=bornOp.setBgSlow)
    intOp = Op.NonLinearOperator(int,int)
    nlOp = Op.CombNonlinearOp(intOp,nlOp)
    ######## prepare the data fitting term ############

    ######## prepare the FWI problem #############
    L2ProbReg = Prblm.ProblemL2NonLinearReg(sub_slow,data,nlOp,eps[i],reg_op=regOp)
    ######## find eps #############
    # L2ProbReg.epsilon = eps[i] * L2ProbReg.estimate_epsilon(verbose=True)
    # print(L2ProbReg.epsilon)
    ######## find eps #############
    Stop  = Stopper.BasicStopper(niter=niter)
    stepper = Stepper.CvSrchStep(alpha_min=1e-7)
    LBFGSsolver = Solver.LBFGSsolver(Stop,stepper=stepper,m_steps=50,logger=logger("Log/lanc5_%s_%s-%s-%s_%s.txt" % (a,ns[i][0],ns[i][1],ns[i][2],eps[i])))
    LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
    			flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/lanc5_%s-%s-%s_%s-%s-%s_%s" % (ns[i][0],ns[i][1],ns[i][2],a[i][0],a[i][1],a[i][2],eps[i]) )
    LBFGSsolver.run(L2ProbReg,verbose=True,restart=False)

    sub_slow = genericIO.defaultIO.getVector("Inv/lanc5_%s-%s-%s_%s-%s-%s_%s_inv_mod.H" % (ns[i][0],ns[i][1],ns[i][2],a[i][0],a[i][1],a[i][2],eps[i]) )
    int.forward(False,sub_slow,slow)
    slow.writeVec("Inv/lanc5_%s-%s-%s_%s-%s-%s_%s_inv_mod_int.H" % (ns[i][0],ns[i][1],ns[i][2],a[i][0],a[i][1],a[i][2],eps[i]))
