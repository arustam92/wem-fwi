import SepVector
import Hypercube
import genericIO
import pyOperator as Op
import WEM
import Operator

import numpy as np


pre = genericIO.defaultIO.getVector("Inv/1800_lin_3_100_inv_mod.H")
full = genericIO.defaultIO.getVector("Vel/constSlowBg1800_ext.H")
fmax = full.getHyper().getAxis(4).o + full.getHyper().getAxis(4).n*full.getHyper().getAxis(4).d
points = np.linspace(1,fmax,pre.getHyper().getAxis(4).n).tolist()


lin = Operator.LinearInterpolation(pre,full,points=points)
lin.forward(False,pre,full)

full.writeVec("Inv/1800_lin_3_100_inv_mod_interpol.H")
