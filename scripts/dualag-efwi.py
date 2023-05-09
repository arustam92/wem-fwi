import numpy as np
import SepVector
import Hypercube
import genericIO
import newExtWEM
import Operator

import pyVector as Vec
import pyOperator as Op
import pyNonLinearSolver as Solver
import pySolverConstrained as ConSolver
import pyLinearSolver as LinearSolver
import pyProblem as Prblm
import pyProblemConstrained as CPrblm
import pyStopper as Stopper
import pyStepper as Stepper
from sys_util import logger
import json
import sys

for p in sys.argv[1:]:
    if ('parfile' in p):
        parfile = p.split('parfile=')[1]
        with open(parfile) as f:
            par = json.load(f)
    else:
      raise Exception("Must provide parfile!")

with open(par['geom']) as f:
    geom = json.load(f)
parObj = genericIO.pythonParams(geom).getCpp()

wave = genericIO.defaultIO.getVector(par['wavelet'])
data = genericIO.defaultIO.getVector(par['data'])
model = genericIO.defaultIO.getVector(par['start_model'])

######## prepare the precondioned variable #############
if 'decon' in par and par['decon'] == 1:
    mod = model.clone()
    decon = Operator.Decon(model, model, wave, eps=1, f0=geom["fmin"])
    deconOp = Op.NonLinearOperator(decon,decon)
    LinStop  = Stopper.BasicStopper(niter=100)
    CGsolver = LinearSolver.LCGsolver(LinStop)
    L2Prob = Prblm.ProblemL2Linear(model,mod,decon)
    CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
                flush_memory=True,iter_buffer_size=1,save_model=True,prefix='temp')
    CGsolver.run(L2Prob,verbose=True)

######## prepare the precondioned variable #############
model_pre = model
if ('pre' in par):
    if ('leaky' in par['pre']['type']):
        model_pre = model.clone()
        model_pre.zero()
        int = Operator.LeakyIntegration(model_pre,model,alpha=par['pre']['alpha'])
    if ('double-' in par['pre']['type']):
        model_pre = model.clone()
        model_pre.zero()
        int = Operator.DoubleLeakyIntegration(model_pre,model,alpha=par['pre']['alpha'])
    else:
        ax = model.getHyper().axes
        ns = par['pre']['ns']
        ds = []
        for j in range(len(ns)):
            ds.append((ax[j].n-1)*ax[j].d / (ns[j] - 1))
        model_pre = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=[ax[0].o,ax[1].o,ax[2].o]),storage='dataComplex')
        model_pre.zero()
        if ("lanczos" in par['pre']['type'] ):
            int = Operator.LanczosInterpolation3D(model_pre,model,a=par['pre']['a'],taper=par['pre']['taper'])
        if ('spline' in par['pre']['type']):
            int = Operator.Spline3D(model_pre,model,type=par['pre']['type'])
    # else:
    #     raise Exception('Unknown intepolation!')
    intOp = Op.NonLinearOperator(int,int)
    LinStop  = Stopper.BasicStopper(niter=par['pre']['niter'])
    CGsolver = LinearSolver.LCGsolver(LinStop)
    L2Prob = Prblm.ProblemL2Linear(model_pre,model,int)
    CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
                flush_memory=True,iter_buffer_size=1,save_model=True,prefix='temp')
    CGsolver.run(L2Prob,verbose=True)

# rect = [1,1,5]
# repeat = 1
# smoothOp = Operator.Smooth(model,model,rect,repeat,mod='all')
smoothMod = Operator.Taper(model,model,tap=[10],axes=[1],zero=[5])


bornOp1 = newExtWEM.extBorn(model,data,wave,parObj)
wemOp = newExtWEM.extWEM(model,data,wave,parObj)
######## prepare the data tapering #############
data_tap = data
if ('data_taper' in par):
    data_tap = data.clone()
    tap = par['data_taper']['tap']
    zero = par['data_taper']['zero']
    sm1 = Operator.Taper(data,data,tap=[tap[0]],axes=[0],zero=[zero[0]])
    sm2 = Operator.Taper(data,data,tap=[tap[1]],axes=[1],zero=[zero[1]])
    dataSm = Op.ChainOperator(sm1,sm2)
    if ('scale_offset' in par['data_taper']):
        scaleOp = Operator.ScaleOffset(data,data,geom)
        dataSm = Op.ChainOperator(dataSm,scaleOp)
    wemOp = Op.ChainOperator(wemOp,dataSm)
    bornOp = Op.ChainOperator(bornOp1,dataSm)
    dataSm.forward(False,data,data_tap)

if ('scale_offset' in par):
    scaleOp = Operator.ScaleOffset(data,data,par)
    dataScale = Op.ChainOperator(sm1,sm2)
    wemOp = Op.ChainOperator(wemOp,dataSm)
    bornOp = Op.ChainOperator(bornOp1,dataSm)
    scaleOp.forward(False,data,data_tap)

if ('grad_pre' in par):
    if ('leaky' in par['grad_pre']['type']):
        int = Operator.LeakyIntegration(model,model,alpha=par['grad_pre']['alpha'])
    if ('double-' in par['grad_pre']['type']):
        int = Operator.DoubleLeakyIntegration(model,model,alpha=par['grad_pre']['alpha'])
    bornOp = Op.ChainOperator(int,bornOp)

bornOpSmooth = Op.ChainOperator(smoothMod,bornOp)
nlOp = Op.NonLinearOperator(wemOp,bornOpSmooth,set_background_func=bornOp1.setBgSlow)

######## test #############
# hyp = Operator.HyperbolicPenalty(model,model,l=1,tau=1e-6)
# sofclip = Operator.Softclip(model,model,l=1,tau=1e-6)
# hypOp = Op.NonLinearOperator(hyp,sofclip,set_background_func=sofclip.setBg)
# nlOp = Op.CombNonlinearOp(hypOp,nlOp)
######## test #############

######## test #############
# alpha = 0.005
# expNl = Operator.ComplexExp(model,model)
# expLin = Operator.LinComplexExp(model,model)
# # Regularization operators
# expOp = Op.NonLinearOperator(expNl,expLin,set_background_func=expLin.setBg)
# nlOp = Op.CombNonlinearOp(expOp, nlOp)
######## test #############
if 'decon' in par and par['decon'] == 1:
    nlOp = Op.CombNonlinearOp(deconOp,nlOp)

if 'hilbert' in par and par['hilbert'] == 1:
    hilbert = Operator.Hilbert(model,model)
    hilbertOp = Op.NonLinearOperator(hilbert,hilbert)
    nlOp = Op.CombNonlinearOp(hilbertOp,nlOp)

if ('pre' in par):
    nlOp = Op.CombNonlinearOp(intOp,nlOp)

if ('av_space' in par and par['av_space'] == 1):
    stack = Operator.StackAndSpread(model_pre,model_pre)
    stackOp = Op.NonLinearOperator(stack,stack)
    nlOp = Op.CombNonlinearOp(stackOp, nlOp)
######## prepare the regularization term #############
der = Operator.Derivative(model_pre,model_pre,which=1,order=4)
# Regularization operators
dsoNlOp = Op.NonLinearOperator(der,der)
# if ('pre' in par): dsoNlOp = Op.CombNonlinearOp(intOp,dsoNlOp)

dual = None
if ("start_dual" in par):
    dual = genericIO.defaultIO.getVector(par['start_dual'])
m_zero = model_pre.clone()
m_zero.zero()
L2Prob = CPrblm.ProblemAugLagrangian(model_pre,m_zero,dsoNlOp,0,eq_op=nlOp,eq_rhs=data_tap,dual_prior=dual)

#Create stopper
Stop  = Stopper.BasicStopper(niter=par['niter'])
#Create solver
rho = [50,50,50]
log = par['non_linear_solver']['prefix'].split('/')
stepper = Stepper.CvSrchStep(maxfev=5)
LBFGSsolver = Solver.NLCGsolver(Stop,stepper=stepper,logger=logger("Log/%s.txt" % log[-1]))
LBFGSsolver.setDefaults(**par['non_linear_solver'])
AugLagSolver = ConSolver.AugLagrangianSolver(LBFGSsolver,rho=par['epsilon'],p_rho=par["p_rho"],constraint_tol=par['c_tol'],m_rho=par["m_rho"], outer=par["outer"])
AugLagSolver.run(L2Prob,verbose=True, restart=False)
# Running the solver
