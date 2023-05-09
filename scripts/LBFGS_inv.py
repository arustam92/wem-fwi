import json

import SepVector
import Hypercube
import genericIO
import WEM
import Operator

import pyVector as Vec
import pyOperator as Op
import pyNonLinearSolver as Solver
import pyStepper as Stepper
import pyProblem as Prblm
import pyStopper as Stopper
from sys_util import logger


D = "./Dat"

with open('Par/geom.json') as f:
    par = json.load(f)

bands = ['0','3-5','3-10','3-15']

max = [0,5,10,15,20]
par["fmin"] = 3
for i in range(1,len(bands)):

    par["fmax"] = max[i] + 10
    parObj = genericIO.pythonParams(par).getCpp()

    slow = genericIO.defaultIO.getVector("Inv/wrong_marm_%s_inv_mod.H" % bands[i-1])
    wave = genericIO.defaultIO.getVector("Wav/MarmWave_%s.H" % bands[i])
    data = genericIO.defaultIO.getVector("Dat/test_%s.H" % bands[i])

    # create a WEM propagator
    rect=([1,2,10],[1,2,10])
    repeat=1
    smoothOp = Operator.Smooth(slow,slow,rect,repeat,mod='top')
    bornOp = WEM.Born(slow,data,wave,parObj)
    bornOpSmooth = Op.ChainOperator(smoothOp,bornOp)
    wemOp = WEM.WEM(slow,data,parObj,wave)
    nlOp = Op.NonLinearOperator(wemOp,bornOpSmooth,set_background_func=bornOp.setBgSlow)

    L2Prob = Prblm.ProblemL2NonLinear(slow,data,nlOp)
    #Create stopper
    niter = 50
    Stop  = Stopper.BasicStopper(niter=niter)
    #Create solver
    stepper = Stepper.CvSrchStep(alpha_min=1e-7)
    LBFGSsolver = Solver.LBFGSsolver(Stop,stepper=stepper,logger=logger("LBFGS_%sHz.txt" % bands[i]))
    #Running the solver
    LBFGSsolver.setDefaults(iter_sampling=1,save_obj=True,save_res=True,save_grad=True,
    			flush_memory=True,iter_buffer_size=1,save_model=True,prefix="Inv/wrong_marm_%s" % bands[i])
    LBFGSsolver.run(L2Prob,verbose=True,restart=False)
