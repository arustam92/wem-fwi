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


model = genericIO.defaultIO.getVector(par['start_model'])

hypOp = Operator.HyperbolicPenalty(model,model,l=1,tau=1e-8)

out = model.clone()
hypOp.forward(False,model,out)

out.writeVec("test.H")
