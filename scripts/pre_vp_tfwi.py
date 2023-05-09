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

eps = 10
# slow = genericIO.defaultIO.getVector("Vel/constSlowBg1800.H")
alpha = 0.
niter = 5
subiter = 25
ns = [51,51]   # have to be odd numbers
a = [11,11]

######## prepare the input #############
D = "./"
with open('Par/geom.json') as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()
wave = genericIO.defaultIO.getVector(D+"Wav/constWav.H")
data = genericIO.defaultIO.getVector(D+"Dat/constData_tap.H")
######## prepare the input #############
model = genericIO.defaultIO.getVector("Vel/constSlowBg1800.H")
refl = genericIO.defaultIO.getVector("Vel/constSlowBg1800_ext.H")
refl.zero()
ax = model.getHyper().axes
# model = newExtWEM.createExtModel(slow,0.5,25)

ds = []
for j in range(2):
    ds.append((ax[j].n-1)*ax[j].d / (ns[j] - 1))

sub_slow = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=[ax[0].o,ax[1].o]),storage='dataComplex')
######## prepare the precondioned variable #############
sub_slow.zero()
int = Operator.LanczosInterpolation2D(sub_slow,model,a=a)
intOp = Op.NonLinearOperator(int,int)
LinStop  = Stopper.BasicStopper(niter=25)
CGsolver = LinearSolver.LCGsolver(LinStop)
L2Prob = Prblm.ProblemL2Linear(sub_slow,model,int)
CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
            flush_memory=True,iter_buffer_size=1,save_model=True,prefix="temp")
CGsolver.run(L2Prob,verbose=True)

######## prepare the regularization term #############
der = Operator.Derivative(refl,refl)
# Regularization operators
dsoNlJac = Op.ZeroOp(model,refl)
dsoNlDummy = Op.ZeroOp(model,refl)
dsoNlOp = Op.NonLinearOperator(dsoNlDummy,dsoNlJac)
dsoNlOp = Op.CombNonlinearOp(intOp,dsoNlOp)
# Variable projection operator for the regularization term
vpRegOp = Op.VpOperator(dsoNlOp,der,Op.dummy_set_background,Op.dummy_set_background)
######## prepare the regularization term #############

######## prepare the data fitting term #############
smoothOp = Operator.Taper(refl,refl,tap=[5],axes=[1],zero=[5])
smoothOp2 = Operator.Taper(model,model,tap=[5],axes=[0],zero=[5])

bornTomo = WEM.BornTomo(model,data,wave,parObj)
bornTomoSm = Op.ChainOperator(smoothOp2,bornTomo)
bornRefl = WEM.BornRefl(refl,model,data,wave,parObj)
bornRefl.add_prec(int)
bornReflSm = Op.ChainOperator(smoothOp,bornRefl)
# bornOpSm = Op.ChainOperator(realSmooth,bornOpSm)
wemOp = WEM.WEM(model,data,wave,parObj)
nlOp = Op.NonLinearOperator(wemOp,bornTomoSm,set_background_func=bornTomo.setBgSlow)
nlOp = Op.CombNonlinearOp(intOp,nlOp)
vpOp = Op.VpOperator(nlOp,bornReflSm,bornRefl.setBgSlow,bornTomo.setBgRefl)
# intOp = Op.NonLinearOperator(int,int)
# nlOp = Op.CombNonlinearOp(intOp,nlOp)
######## prepare the data fitting term #############

######## prepare the linear solver #############
LinStop  = Stopper.BasicStopper(niter=subiter)
CGsolver = LinearSolver.LCGsolver(LinStop)
L2Prob = Prblm.ProblemL2Linear(refl,data,bornReflSm)
CGsolver.setDefaults(iter_sampling=1,save_obj=False,save_res=False,save_grad=False,
			flush_memory=True,iter_buffer_size=1,save_model=False )
######## prepare the precondioned variable #############

######## prepare the FWI problem #############
VpProbReg = Prblm.ProblemL2VpReg(sub_slow,refl,vpOp,data,CGsolver,h_op_reg=vpRegOp,epsilon=eps)
Stop  = Stopper.BasicStopper(niter=niter)
stepper = Stepper.CvSrchStep(alpha_min=1e-7)
LBFGSsolver = Solver.LBFGSsolver(Stop,stepper=stepper,m_steps=50,logger=logger("Log/1800_pre_vp_eps_%s.txt" % (eps)))
LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
            flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/1800_pre_vp_eps_%s" % (eps) )
LBFGSsolver.run(VpProbReg,verbose=True,restart=False)
# bornRefl.adjoint(False,refl,data)
# bornTomo.setBgRefl(refl)
# bornTomo.setBgSlow(model)
# bornTomo.adjoint(False,model,data)
# model.writeVec('temp.H')

# slow2 = genericIO.defaultIO.getVector("Inv/1800_lin_%s_%s_inv_mod.H" % (eps,len(points[i])) )
# int.forward(False,slow2,slow)
# del slow2
