import SepVector
import Hypercube
import genericIO
import Operator
import ExtWEM
import numpy as np
from scipy.ndimage import gaussian_filter
import os

PATH = "/home/user/project/projects/extWEM/"
par = ["","par=geom.par"]

io=genericIO.pyGenericIO.ioModes(par)
parObj=io.getDefaultIO().getParamObj()

# sigs_path = "gauss.H"
slow = genericIO.defaultIO.getVector(PATH + "Vel/extBg.H")
blurred = slow.clone()
n1=slow.getHyper().getAxis(1).n
n2=slow.getHyper().getAxis(2).n
n3=slow.getHyper().getAxis(3).n
blurNd = blurred.getNdArray()
out = blurred.clone()
outNd = out.getNdArray()

# l = [15,15]
# repeat=3
# sm = Operator.Smooth(slow,blurred,l,repeat)
# sm.forward(False,slow,blurred)



num = [5,10]
sigmas = [0.03*n1,0.07*n1]
amp = [0.1,0.3]

sign = 1
realizations = 100
spike = np.zeros((n2,n1))
for i in range(realizations):
    outNd[:] = blurNd[:]
    n = np.random.randint(num[0],num[1])
    for j in range(n):
        spike[:] = 0
        sig = np.random.uniform(sigmas[0],sigmas[1])
        perc = np.random.uniform(amp[0],amp[1])
        ind1 = np.random.randint(0,n1)
        ind2 = np.random.randint(0,n2)
        spike[ind2,ind1] = sign
        spike = gaussian_filter(spike,sigma=(sig,sig))
        a = perc * outNd.real[0,ind2,ind1]
        spike = a * spike/np.amax(np.abs(spike))
        outNd += spike
        sign *= -1
    print("Realization %s" %i )
    genericIO.defaultIO.writeVector("Vel2/rand_%s.H" % i,out)
