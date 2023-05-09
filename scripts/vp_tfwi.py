import SepVector
import Hypercube
import genericIO
import WEM
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
import sys
import copy

for p in sys.argv[1:]:
    if ('parfile' in p):
        parfile = p.split('parfile=')[1]
        with open(parfile) as f:
            par = json.load(f)
    else:
      raise Exception("Must provide parfile!")

eps = par['epsilon']
with open(par['geom']) as f:
    geom = json.load(f)
######## prepare the input #############

parObj = genericIO.pythonParams(geom).getCpp()
wave = genericIO.defaultIO.getVector(par['wavelet'])
data = genericIO.defaultIO.getVector(par['data'])
######## prepare the input #############
model = genericIO.defaultIO.getVector(par['start_model'])
refl = genericIO.defaultIO.getVector(par['start_ext_refl'])
refl.zero()
model_pre = model
######## prepare the precondioned variable #############
if ('pre' in par):
    ax = model.getHyper().axes
    ns = par['pre']['ns']
    ds = []
    for j in range(2):
        ds.append((ax[j].n-1)*ax[j].d / (ns[j] - 1))
    model_pre = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=[ax[0].o,ax[1].o]),storage='dataComplex')
    model_pre.zero()
    if (par['pre']['type'] == 'lanczos'):
        int = Operator.LanczosInterpolation2D(model_pre,model,a=par['pre']['a'])
    if ('spline' in par['pre']['type']):
        int = Operator.Spline2D(model_pre,model,type=par['pre']['type'])
    else:
        raise Exception('Unknown intepolation!')
    intOp = Op.NonLinearOperator(int,int)
    LinStop  = Stopper.BasicStopper(niter=par['pre']['niter'])
    CGsolver = LinearSolver.LCGsolver(LinStop)
    L2Prob = Prblm.ProblemL2Linear(model_pre,model,int)
    CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
                flush_memory=True,iter_buffer_size=1,save_model=True,prefix='temp')
    CGsolver.run(L2Prob,verbose=True)

######## pad the model #############
refl_pad = refl
model_pad = model
model_pad = [[0,0],[0,0]]
reflPad = [[0,0],[0,0]]
if ('tap' in geom):
    if ('model_pad' in par): model_pad = par['model_pad']
    ds = np.array([model.getHyper().getAxis(1).d, model.getHyper().getAxis(2).d])
    os = np.array([model.getHyper().getAxis(1).o, model.getHyper().getAxis(2).o])
    ns = np.array([model.getHyper().getAxis(1).n, model.getHyper().getAxis(2).n])
    beg = np.array([geom['tap'],0]) + np.array(model_pad[0])
    end = np.array([geom['tap'],0]) + np.array(model_pad[1])
    os = (os - beg*ds).tolist()
    ns = (ns + beg + end).tolist()
    ds = ds.tolist()
    model_pad = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    ds.append(refl.getHyper().getAxis(3).d)
    os.append(refl.getHyper().getAxis(3).o)
    ns.append(refl.getHyper().getAxis(3).n)
    refl_pad = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    modPadOp = Operator.Pad2d(model,model_pad,beg=beg[::-1],end=end[::-1])
    modPadOp.forward(False,model,model_pad)
    modPadOpNL = Op.NonLinearOperator(modPadOp,modPadOp)

    if ('refl_pad' in par): reflPad = par['refl_pad']
    beg = np.array([geom['tap'],0]) + np.array(reflPad[0])
    end = np.array([geom['tap'],0]) + np.array(reflPad[1])
    reflPadOp = Operator.Pad2d(refl,refl_pad,beg=beg[::-1],end=end[::-1])

######## prepare the regularization term #############
der = Operator.Derivative(refl_pad,refl_pad)
if ('tap' in geom): der = Op.ChainOperator(reflPadOp,der)
# Regularization operators
dsoNlJac = Op.ZeroOp(model_pad,refl_pad)
dsoNlDummy = Op.ZeroOp(model_pad,refl_pad)
dsoNlOp = Op.NonLinearOperator(dsoNlDummy,dsoNlJac)
if ('tap' in geom): dsoNlOp = Op.CombNonlinearOp(modPadOpNL,dsoNlOp)
if ('pre' in par): dsoNlOp = Op.CombNonlinearOp(intOp,dsoNlOp)
# Variable projection operator for the regularization term
vpRegOp = Op.VpOperator(dsoNlOp,der,Op.dummy_set_background,Op.dummy_set_background)
######## prepare the regularization term #############

######## prepare the data fitting term #############
smoothOp = Operator.Taper(refl_pad,refl_pad,tap=[5],axes=[1],zero=[5])
smoothOp2 = Operator.Taper(model_pad,model_pad,tap=[10],axes=[0],zero=[10])

bornTomo = WEM.BornTomo(model_pad,data,wave,parObj)
bornTomoSm = Op.ChainOperator(smoothOp2,bornTomo)

bornRefl = WEM.BornRefl(refl_pad,model_pad,data,wave,parObj)
if ('pre' in par): bornRefl.add_prec(int)
if ('tap' in geom):
    bornRefl.add_pad(modPadOp)
    bornTomo.add_pad(reflPadOp)
bornReflSm = Op.ChainOperator(smoothOp,bornRefl)
if ('tap' in geom): bornReflSm = Op.ChainOperator(reflPadOp,bornReflSm)

bornFull = WEM.Born(model_pad,data,wave,parObj)
bornFullSm = Op.ChainOperator(smoothOp2,bornFull)
wemOp = WEM.WEM(model_pad,data,wave,parObj)

wemOp = Op.NonLinearOperator(wemOp,bornFull,set_background_func=bornFull.setBgSlow)
nlOp = Op.NonLinearOperator(wemOp,bornTomo,set_background_func=bornTomo.setBgSlow)
if ('tap' in geom):
    wemOp = Op.CombNonlinearOp(modPadOpNL,wemOp)
    nlOp = Op.CombNonlinearOp(modPadOpNL,nlOp)

######## prepare the data tapering #############
data_tap = data
if ('data_taper' in par):
    data_tap = data.clone()
    tap = par['data_taper']['tap']
    zero = par['data_taper']['zero']
    sm1 = Operator.Taper(data,data,tap=[tap[0]],axes=[0],zero=[zero[0]])
    sm2 = Operator.Taper(data,data,tap=[tap[1]],axes=[1],zero=[zero[1]])
    dataSm = Op.ChainOperator(sm1,sm2)
    dataSmNl = Op.NonLinearOperator(dataSm,dataSm)
    bornReflSm = Op.ChainOperator(bornReflSm,dataSm)
    bornTomoSm = Op.ChainOperator(bornTomoSm,dataSm)
    wemOp = Op.CombNonlinearOp(wemOp,dataSmNl)
    dataSm.forward(False,data,data_tap)

if ('pre' in par):
    wemOp = Op.CombNonlinearOp(intOp,wemOp)
    nlOp = Op.CombNonlinearOp(intOp,nlOp)
vpOp = Op.VpOperator(nlOp,bornReflSm,bornRefl.setBgSlow,bornTomo.setBgRefl)
#
######## prepare the linear solver #############
LinStop  = Stopper.BasicStopper(niter=par['subiter'])
log = par['linear_solver']['prefix'].split('/')
CGsolver = LinearSolver.LCGsolver(LinStop,logger=logger("Log/%s.txt"%par['save_as_prefix']))
L2Prob = Prblm.ProblemL2Linear(refl,data,bornReflSm)
CGsolver.setDefaults(**par['linear_solver'])
######## prepare the precondioned variable #############

######## prepare the FWI problem #############
g_op = None
if (par["g_op"]): g_op = wemOp
VpProbReg = Prblm.ProblemL2VpReg(model_pre,refl,vpOp, data_tap,CGsolver,g_op=g_op,h_op_reg=vpRegOp,epsilon=eps)
Stop  = Stopper.BasicStopper(niter=par['niter'])
stepper = Stepper.CvSrchStep(alpha_min=1e-7)
log = par['non_linear_solver']['prefix'].split('/')
LBFGSsolver = Solver.LBFGSsolver(Stop,stepper=stepper,m_steps=50,logger=logger("Log/%s.txt" % log[-1]))
LBFGSsolver.setDefaults(**par['non_linear_solver'])
LBFGSsolver.run(VpProbReg,verbose=True,restart=False)
