import SepVector
import Hypercube
import genericIO
import Operator
import pyOperator as Op
import json
import sys

a = 0
type = 0
taper = 0
for par in sys.argv[1:]:
    if ('in=' in par):
        input = par.split('in=')[1]
    if ('space=' in par):
        space = par.split('space=')[1]
    if ('a=' in par):
        a = par.split('a=')[1]
        a = list(map(int,a.split(',')))
    if ('type=' in par):
        type = par.split('type=')[1]
    if ('taper=' in par):
        taper = par.split('taper=')[1]
        taper = list(map(float,taper.split(',')))


out = list(input)
out.insert(-2,'_int')
out = ''.join(out)

slow = genericIO.defaultIO.getVector(input)
space = genericIO.defaultIO.getVector(space)
slowNd = slow.getNdArray()
spaceNd = space.getNdArray()

sax = space.getHyper().axes
ax = slow.getHyper().axes
ns = []
os = []
ds = []
sns = []
sos = []
sds = []
for i in range(len(sax)):
    ns.append(ax[i].n)
    os.append(ax[i].o)
    ds.append(ax[i].d)
    sns.append(sax[i].n)
    sos.append(sax[i].o)
    sds.append(sax[i].d)


smol = SepVector.getSepVector(Hypercube.hypercube(ns=ns, os=os, ds=ds),storage='dataComplex')
smolNd = smol.getNdArray()

if (len(ns) == 3):
    if a:
        lanc = Operator.LanczosInterpolation3D(smol,space,a=a,taper=taper)
    elif type:
        lanc = Operator.Spline3D(smol,space,type=type,taper=taper)
else:
    if a:
        lanc = Operator.LanczosInterpolation2D(smol,space,a=a)
    elif type:
        lanc = Operator.Spline2D(smol,space,type=type)
if len(slowNd.shape) == len(sax):
    output = space.clone()
    lanc.forward(False,slow,output)
else:
    sns.append(ax[-1].n)
    sos.append(ax[-1].o)
    sds.append(ax[-1].d)
    output = SepVector.getSepVector(Hypercube.hypercube(ns=sns, os=sos, ds=sds),storage='dataComplex')
    outputNd = output.getNdArray()
    for iter in range(slowNd.shape[0]):
        print(iter)
        smolNd[:] = slowNd[iter,:]
        lanc.forward(False,smol,space)
        outputNd[iter,:] = spaceNd[:]

output.writeVec(out)
