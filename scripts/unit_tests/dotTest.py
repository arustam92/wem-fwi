import SepVector
import Hypercube
import genericIO
import pyOperator as Op
import newExtWEM
import WEM
import Operator
import json
import os
import sys

ext = 0
for par in sys.argv[1:]:
    if ('model' in par):
        model = par.split('model=')[1]
    if ('data' in par):
        data = par.split('data=')[1]
    # if ('wave' in par):
    #     wave = par.split('wave=')[1]
    # if ('geom' in par):
    #     geom = par.split('geom=')[1]


# with open(geom) as f:
#     par = json.load(f)
# parObj = genericIO.pythonParams(par).getCpp()
model = genericIO.defaultIO.getVector(model)
# wave = genericIO.defaultIO.getVector(wave)
data = genericIO.defaultIO.getVector(data)

#
# bornOp = WEM.BornRefl(model,data,wave,parObj)
# bornOp.dotTest()

int = Operator.LanczosInterpolation2D(model,data,a=[3,3])
int.dotTest(verbose=True)
# realizations = 100
# for i in range(realizations):
#     print("Realization %s" % i)
#     slowness = genericIO.defaultIO.getVector("Vel2/rand_%s.H" % i)
#     bornOp.setBgSlow(slowness)
#     bornOp.adjoint(False,image,data)
#     genericIO.defaultIO.writeVector("Image2/mig_rand_%s.H" % i,image)
#     os.system('python scripts/f2tau.py "Image2/mig_rand_%s.H" "Image2/mig_rand_%s_tau.H"' % (i,i))

#
