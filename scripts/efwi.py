import numpy as np
import SepVector
import Hypercube
import genericIO
import newExtWEM
import Operator

import pyVector as Vec
import pyOperator as Op
import pyLinearSolver as LinearSolver
import pyNonLinearSolver as Solver
import pyProblem as Prblm
import pyStopper as Stopper
import pyStepper as Stepper
from sys_util import logger
import json
import sys

import dask_util
from __pyDaskVector import DaskVector, DaskSuperVector
from __pyDaskOperator import DaskOperator


for p in sys.argv[1:]:
    if ('parfile' in p):
        parfile = p.split('parfile=')[1]
        with open(parfile) as f:
            par = json.load(f)
    else:
      raise Exception("Must provide parfile!")

client = dask_util.load(par['dask_client'])
chunks = par['chunks']

with open(par['geom']) as f:
    geom = json.load(f)
parObj = genericIO.pythonParams(geom)

wave = genericIO.defaultIO.getVector(par['wavelet'])

data = genericIO.defaultIO.getVector(par['data'])
data = DaskVector(client, from_vector=data, chunks=chunks)

loc_model = genericIO.defaultIO.getVector(par['start_model'])
loc_den = genericIO.defaultIO.getVector(par['start_density'])
model = DaskSuperVector(loc_model, loc_den)

######## prepare the precondioned variable #############
model_pre = model
if ('pre' in par):
    vecs_pre = []
    ops_pre = []
    for i, mod in enumerate([loc_model, loc_den]):
        ax = mod.getHyper().axes
        ns = par['pre']['ns'][i]
        ds = []
        for j in range(len(ns)):
            ds.append((ax[j].n-1)*ax[j].d / (ns[j] - 1))
        mod_pre = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=[ax[0].o,ax[1].o,ax[2].o]),storage='dataComplex')
        mod_pre.zero()
        vecs_pre.append(mod_pre)

        if (par['pre']['type'][i] == 'lanczos'):
            int = Operator.LanczosInterpolation3D(mod_pre, mod, a=par['pre']['a'][i], taper=par['pre']['taper'][i])
        elif ('spline' in par['pre']['type'][i]):
            int = Operator.Spline3D(mod_pre, mod, type=par['pre']['type'][i])
        else:
            raise Exception('Unknown intepolation!')
        ops_pre.append(int)
    
    int = Op.Dstack(*ops_pre)
    model_pre = DaskSuperVector(*vecs_pre)
    int.setDomainRange(model_pre, model)

    intOp = Op.NonLinearOperator(int,int)
    LinStop  = Stopper.BasicStopper(niter=par['pre']['niter'])
    CGsolver = LinearSolver.LCGsolver(LinStop)
    L2Prob = Prblm.ProblemL2Linear(model_pre,model,int)
    CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
                flush_memory=True,iter_buffer_size=1,save_model=True,prefix='temp')
    CGsolver.run(L2Prob,verbose=True)

# # # symmetric imaging condition
model_pad = model
# if 'symmetric' in par and par['symmetric'] == 1:
#     padOp = Operator.Pad.from_params(loc_model, beg=[0,0,0], end=[model.getHyper().getAxis(3).n, 0, 0], mode='edge')
#     # prepare for interpolation
#     model_pad = padOp.get_padded()
#     ax = model_pad.getHyper().axes.copy()
#     ax[-1].n *= 2
#     ax[-1].d /= 2

#     model_int = SepVector.getSepVector(axes=ax,storage='dataComplex')
#     model_int = DaskVector(client, from_vector=model_int, chunks=(1,1,1))
#     model_pad = DaskVector(client, from_vector=model_pad, chunks=(1,1,1))
    
#     padOp = DaskOperator(client, Operator.Pad, model, model_pad, beg=[0,0,0], end=[model.getHyper().getAxis(3).n, 0, 0], mode='constant')
#     splineOp = DaskOperator(client, Operator.Spline3D, model_pad, model_int, type="CR-spline")
#     # masking operator
#     maskOp = DaskOperator(client, Operator.Pad, model_pad, model_int, beg=[0,0,0], end=[model_pad.getHyper().getAxis(3).n, 0, 0], mode='constant')
#     maskOp = maskOp.H

#     symmOp = Op.ChainOperator(padOp, Op.ChainOperator(splineOp, maskOp))
#     symmOp.forward(False, model, model_pad)

bornOp = DaskOperator(client, newExtWEM.extBorn, model_pad,data,wave,parObj)
wemOp = DaskOperator(client, newExtWEM.extWEM, model_pad,data,wave,parObj)
# if ('grad_pre' in par):
#     if ('leaky' in par['grad_pre']['type']):
#         gradOp = DaskOperator(client, Operator.LeakyIntegration, model_pad,model_pad,alpha=par['grad_pre']['alpha'])
#     elif ('double-' in par['grad_pre']['type']):
#         gradOp = DaskOperator(client, Operator.DoubleLeakyIntegration, model_pad,model_pad,alpha=par['grad_pre']['alpha'])
#     else:
#         ax = model_pad.getHyper().axes
#         ns = par['grad_pre']['ns']
#         ds = []
#         for j in range(len(ns)):
#             ds.append((ax[j].n-1)*ax[j].d / (ns[j] - 1))
#         grad_pre = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=[ax[0].o,ax[1].o,ax[2].o]),storage='dataComplex')
#         grad_pre = DaskVector(client, from_vector=grad_pre, chunks=(1,1,1))
#         grad_pre.zero()
#         if ("lanczos" in par['grad_pre']['type'] ):
#             grad = DaskOperator(client, Operator.LanczosInterpolation3D, grad_pre,model_pad,a=par['grad_pre']['a'],taper=par['grad_pre']['taper'])
#         if ('spline' in par['grad_pre']['type']):
#             grad = DaskOperator(client, Operator.Spline3D, grad_pre,model_pad,type=par['grad_pre']['type'])
#         gradOp = Op.ChainOperator(grad.H, grad)
#     bornOp = Op.ChainOperator(gradOp,bornOp)

######## prepare the data tapering #############
data_tap = data
# if ('data_taper' in par):
#     data_tap = data.clone()
#     tap = par['data_taper']['tap']
#     zero = par['data_taper']['zero']
#     sm1 = Operator.Taper(data,data,tap=[tap[0]],axes=[0],zero=[zero[0]])
#     sm2 = Operator.Taper(data,data,tap=[tap[1]],axes=[1],zero=[zero[1]])
#     dataSm = Op.ChainOperator(sm1,sm2)
#     if ('scale_offset' in par['data_taper']):
#         scaleOp = Operator.ScaleOffset(data,data,geom)
#         dataSm = Op.ChainOperator(dataSm,scaleOp)
#     wemOp = Op.ChainOperator(wemOp,dataSm)
#     bornOp = Op.ChainOperator(bornOp,dataSm)
#     dataSm.forward(False,data,data_tap)

if ('mask' in par):
    mask_vec = genericIO.defaultIO.getVector(par['mask'])
    mask = DaskSuperVector(mask_vec, mask_vec)
    maskOp = Operator.Mask(model_pad, model_pad, mask)
    bornOp = Op.ChainOperator(maskOp, bornOp)

# if ('scale_offset' in par):
#     data = data_tap
#     scaleOp = Operator.ScaleOffset(data,data, geom, scale=par["scale_offset"])
#     wemOp = Op.ChainOperator(wemOp,scaleOp)
#     bornOp = Op.ChainOperator(bornOp,scaleOp)
#     data_tap = data_tap.clone()
#     scaleOp.forward(False,data,data_tap)

nlOp = Op.NonLinearOperator(wemOp,bornOp)

####### test #############
hyp = DaskOperator(client, Operator.HyperbolicPenalty, model_pad,model_pad,l=1,tau=geom["tau"])
sofclip = DaskOperator(client, Operator.Softclip, model_pad,model_pad,l=1,tau=geom["tau"])
hypOp = Op.NonLinearOperator(hyp,sofclip)
nlOp = Op.CombNonlinearOp(hypOp,nlOp)
####### test #############

####### test #############
if "vmin" in par and "vmax" in par:
    clipOp = Operator.SoftMinMax(model_pad, model_pad, 1/par["vmax"]**2, 1/par["vmin"]**2, tau=[0.01,0.001])
    nlOp = Op.CombNonlinearOp(clipOp,nlOp)
####### test #############

# if 'decon' in par and par['decon'] == 1:
#     temp = data_tap.clone()
#     decon = DaskOperator(client, Operator.Decon, data_tap, data_tap, wave, eps=1e-2)
#     deconOp = Op.NonLinearOperator(decon,decon)
#     decon.forward(False,temp,data_tap)
#     nlOp = Op.CombNonlinearOp(nlOp,deconOp)

# if 'symmetric' in par and par['symmetric'] == 1:
#     nlOp = Op.CombNonlinearOp(Op.NonLinearOperator(symmOp), nlOp)

if ('pre' in par):
    nlOp = Op.CombNonlinearOp(intOp,nlOp)

# if 'av_space' in par and par["av_space"] == 1:
#     stackOp = Operator.StackAndSpread(model,model)
#     nlOp = Op.CombNonlinearOp(Op.NonLinearOperator(stackOp),nlOp)

######## prepare the regularization term #############
dt = wave.getHyper().getAxis(1).d
tmax = (wave.getHyper().getAxis(1).n - 1)*dt + wave.getHyper().getAxis(1).o
der = DaskOperator(client, Operator.Derivative, model,model,which=1,order='exact',f0=geom["fmin"],tmax=tmax/2,dt=dt,alpha=par["alpha"],beta=par["beta"],mode=par['der_mode'])

# Regularization operators
dsoNlOp = Op.NonLinearOperator(der,der)
if ('pre' in par): dsoNlOp = Op.CombNonlinearOp(intOp,dsoNlOp)

#Create stopper
Stop  = Stopper.BasicStopper(niter=par['niter'])
#Create solver
rho = [50,50,50]
log = par['non_linear_solver']['prefix'].split('/')
stepper = Stepper.CvSrchStep(maxfev=5)

# if "vmin" in par and "vmax" in par:
#     minBound = model_pre.clone()
#     minBound[0].set(1/par.get('vmax')**2 - 1j)
#     minBound[1].set(-1e6 - 1e6j)
#     maxBound = model_pre.clone()
#     maxBound[0].set(1/par.get('vmin')**2 + 1j)
#     maxBound[1].set(1e6 + 1e6j)
# else:
#     minBound = None
#     maxBound = None

if par['solver'] == 'lbfgs':
    L2Prob = Prblm.ProblemL2NonLinearReg(model_pre,data_tap,nlOp,par["epsilon"],reg_op=dsoNlOp)
    if minBound and maxBound:
        solver = Solver.LBFGSBsolver(Stop,m_steps=10,stepper=stepper,logger=logger("Log/%s.txt" % log[-1]))
    else:
        solver = Solver.LBFGSsolver(Stop,m_steps=10,stepper=stepper,logger=logger("Log/%s.txt" % log[-1]))
    solver.setDefaults(**par['non_linear_solver'])
elif par['solver'] == 'nlcg':
    L2Prob = Prblm.ProblemL2NonLinearReg(model_pre,data_tap,nlOp,par["epsilon"],reg_op=dsoNlOp)
    solver = Solver.NLCGsolver(Stop,beta_type='HZ',logger=logger("Log/%s.txt" % log[-1]))
    solver.setDefaults(**par['non_linear_solver'])

elif par['solver'] == 'auglag':
    import pySolverConstrained as CSolver
    import pyProblemConstrained as CPrblm

    dual = None
    if ("start_dual" in par):
        dual = genericIO.defaultIO.getVector(par['start_dual'])
    L2Prob = CPrblm.ProblemAugLagrangian(model_pre,data_tap,nlOp,0,eq_op=dsoNlOp,dual_prior=dual)
    nlcg = Solver.NLCGsolver(Stop,beta_type='HZ',logger=logger("Log/%s.txt" % log[-1]))
    nlcg.setDefaults(**par['non_linear_solver'])
    solver = CSolver.AugLagrangianSolver(nlcg,rho=par['epsilon'],p_rho=par["p_rho"],constraint_tol=par['c_tol'],m_rho=par["m_rho"], outer=par["outer"], save_dual=par['save_dual'])

elif par['solver'] == 'auglag-dual':
    import pySolverConstrained as CSolver
    import pyProblemConstrained as CPrblm

    if ("start_dual" in par):
        dual = genericIO.defaultIO.getVector(par['start_dual'])
    else:
        dual = model.clone()
        dual.zero()
    L2Prob = CPrblm.ProblemAugLagrangian(model_pre,dual,dsoNlOp,0,eq_op=nlOp,eq_rhs=data_tap)
    nlcg = Solver.NLCGsolver(Stop,beta_type='HZ',logger=logger("Log/%s.txt" % log[-1]))
    nlcg.setDefaults(**par['non_linear_solver'])
    solver = CSolver.AugLagrangianSolver(nlcg,rho=par['epsilon'],p_rho=par["p_rho"],constraint_tol=par['c_tol'],m_rho=par["m_rho"], outer=par["outer"], save_dual=par['save_dual'])

solver.run(L2Prob,verbose=True, restart=False)
