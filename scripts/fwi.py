import numpy as np
import SepVector
import Hypercube
import genericIO
import WEM
import Operator

import pyVector as Vec
import pyOperator as Op
import pyNonLinearSolver as Solver
import pyLinearSolver as LinearSolver
import pyProblem as Prblm
import pyStopper as Stopper
from sys_util import logger
import json
import sys

import dask_util
from __pyDaskVector import DaskVector
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

model = genericIO.defaultIO.getVector(par['start_model'])
model = DaskVector(client, from_vector=model, chunks=(1,1))

######## prepare the precondioned variable #############
model_pre = model
if ('pre' in par):
    ax = model.getHyper().axes
    ns = par['pre']['ns']
    ds = []
    for j in range(2):
        ds.append((ax[j].n-1)*ax[j].d / (ns[j] - 1))
    model_pre = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=[ax[0].o,ax[1].o]),storage='dataComplex')
    model_pre = DaskVector(client, from_vector=model_pre, chunks=(1,1))
    model_pre.zero()
    if (par['pre']['type'] == 'lanczos'):
        int = DaskOperator(client, Operator.LanczosInterpolation2D, model_pre, model, a=par['pre']['a'])
    if ('spline' in par['pre']['type']):
        int = DaskOperator(client, Operator.Spline2D, model_pre, model, type=par['pre']['type'])
    else:
        raise Exception('Unknown intepolation!')
    intOp = Op.NonLinearOperator(int,int)
    LinStop  = Stopper.BasicStopper(niter=par['pre']['niter'])
    CGsolver = LinearSolver.LCGsolver(LinStop)
    L2Prob = Prblm.ProblemL2Linear(model_pre,model,int)
    CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
                flush_memory=True,iter_buffer_size=1,save_model=True,prefix='temp')
    CGsolver.run(L2Prob,verbose=True)

bornOp = DaskOperator(client, WEM.Born, model, data, wave, parObj, mode='real')
wemOp = DaskOperator(client, WEM.WEM, model, data, wave, parObj)
if ('mask' in par):
    mask = genericIO.defaultIO.getVector(par['mask'])
    maskOp = DaskOperator(client, Operator.Mask, model, model, mask)
    bornOp = Op.ChainOperator(maskOp, bornOp)
    
######## prepare the data tapering #############
data_tap = data
# if ('data_taper' in par):
#     data_tap = data.clone()
#     tap = par['data_taper']['tap']
#     zero = par['data_taper']['zero']
#     sm1 = Operator.Taper(data,data,tap=[tap[0]],axes=[0],zero=[zero[0]])
#     sm2 = Operator.Taper(data,data,tap=[tap[1]],axes=[1],zero=[zero[1]])
#     dataSm = Op.ChainOperator(sm1,sm2)
#     wemOp = Op.ChainOperator(wemOp,dataSm)
#     bornOp = Op.ChainOperator(bornOp,dataSm)
#     dataSm.forward(False,data,data_tap)

nlOp = Op.NonLinearOperator(wemOp,bornOp)

if 'decon' in par and par['decon'] == 1:
    temp = data_tap.clone()
    decon = DaskOperator(client, Operator.Decon, data_tap, data_tap, wave, eps=1e-2)
    deconOp = Op.NonLinearOperator(decon,decon)
    decon.forward(False,temp,data_tap)
    nlOp = Op.CombNonlinearOp(nlOp,deconOp)

if ('pre' in par):
    nlOp = Op.CombNonlinearOp(intOp,nlOp)

#Create L2-norm linear problem
niter = par['niter']
if "vmin" in par and "vmax" in par:
    minBound = model_pre.clone()
    minBound.set(1/par.get('vmax')**2)
    maxBound = model_pre.clone()
    maxBound.set(1/par.get('vmin')**2)
else:
    minBound = None
    maxBound = None

L2Prob = Prblm.ProblemL2NonLinear(model_pre,data_tap,nlOp, minBound=minBound, maxBound=maxBound)
#Create stopper
Stop  = Stopper.BasicStopper(niter=niter)
#Create solver
log = par['non_linear_solver']['prefix'].split('/')
if "solver" in par["non_linear_solver"]:
    if par["non_linear_solver"]["solver"] == 'nlcg':
        if "beta" in par["non_linear_solver"]:
            beta_type = par['non_linear_solver']['beta']
            del par['non_linear_solver']['beta']
        else:
            beta_type = 'FR'
        solver = Solver.NLCGsolver(Stop,beta_type=beta_type, logger=logger("Log/%s.txt" % log[-1]))
    del par['non_linear_solver']['solver']
else:
    solver = Solver.LBFGSsolver(Stop,logger=logger("Log/%s.txt" % log[-1]))
#Running the solver
solver.setDefaults(**par['non_linear_solver'])
solver.run(L2Prob,verbose=True,restart=False)
