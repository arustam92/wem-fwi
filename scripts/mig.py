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
spline = None
image = 0
for par in sys.argv[1:]:
    if ('vel' in par):
        vel = par.split('vel=')[1]
    if ('wave' in par):
        wave = par.split('wave=')[1]
    if ('geom' in par):
        geom = par.split('geom=')[1]
    if ('data' in par):
        data = par.split('data=')[1]
    if ('out' in par):
        out = par.split('out=')[1]
    if ('image' in par):
        image = par.split('image=')[1]

with open(geom) as f:
    par = json.load(f)
parObj = genericIO.pythonParams(par).getCpp()
slow = genericIO.defaultIO.getVector(vel)
wave = genericIO.defaultIO.getVector(wave)
data = genericIO.defaultIO.getVector(data)


# out = genericIO.defaultIO.getVector("Dat/mig.H")

if not image:
    image = slow.clone()
else:
    image = genericIO.defaultIO.getVector(image)
#
if (len(slow.getHyper().axes) > 2):
    bornOp = newExtWEM.extBorn(slow,data,wave,parObj)
    if spline:
        ax = slow.getHyper().axes
        ns = list(spline)
        ds = []
        for j in range(len(ns)):
            ds.append((ax[j].n-1)*ax[j].d / (ns[j] - 1))
        model_pre = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=[ax[0].o,ax[1].o,ax[2].o]),storage='dataComplex')
        int = Operator.Spline3D(model_pre,slow,type='CR-spline')
        bornOp = Op.ChainOperator(int,bornOp)
        image = model_pre.clone()
else:
    bornOp = WEM.Born(slow,data,wave,parObj)

# taper = Operator.Taper(slow,slow,tap=[50],axes=[1])

# taper.forward(False,out,slow)
bornOp.setBgSlow(slow)
bornOp.adjoint(False,image,data)
image.writeVec(out)

# realizations = 100
# for i in range(realizations):
#     print("Realization %s" % i)
#     slowness = genericIO.defaultIO.getVector("Vel2/rand_%s.H" % i)
#     bornOp.setBgSlow(slowness)
#     bornOp.adjoint(False,image,data)
#     genericIO.defaultIO.writeVector("Image2/mig_rand_%s.H" % i,image)
#     os.system('python scripts/f2tau.py "Image2/mig_rand_%s.H" "Image2/mig_rand_%s_tau.H"' % (i,i))

#
