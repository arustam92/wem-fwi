import SepVector
import Hypercube
import genericIO
import Operator
import pyOperator as Op
import json
import sys
import numpy as np

for par in sys.argv[1:]:
    if ('in=' in par):
        input = par.split('in=')[1]
    if ('out=' in par):
        out = par.split('out=')[1]

data = genericIO.defaultIO.getVector(input)
output = data.clone()
output2 = data.clone()

mask = np.zeros((data.getHyper().getAxis(2).n,data.getHyper().getAxis(1).n))
x = np.linspace(data.getHyper().getAxis(1).o, (data.getHyper().getAxis(1).n-1)*data.getHyper().getAxis(1).d, data.getHyper().getAxis(1).n)
z = np.linspace(data.getHyper().getAxis(2).o, (data.getHyper().getAxis(2).n-1)*data.getHyper().getAxis(2).d, data.getHyper().getAxis(2).n)
X,Z = np.meshgrid(x,z)
center = [35000,65000]
radius = np.sqrt((X-center[1])**2 + (Z-center[0])**2)
rmax = 27000
mask = np.where(radius > rmax, 0, 1)
taper = Operator.Mask(output,data,mask,7, repeat=5)

# taper = Operator.Mask(output2,data,mask=,sigma=(3,3))
# # taper = Operator.Taper(data,data,tap=[100],axes=[2],zero=[320])
taper.forward(False,data,output2)
# taper = Operator.Taper(data,data,tap=[50],axes=[1],zero=[50])
# taper.forward(False,output,output2)

genericIO.defaultIO.writeVector(out,output2)
