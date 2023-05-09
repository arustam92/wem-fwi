import SepVector
import Hypercube
import genericIO
import numpy as np


nz=100
nx=100
dz=10
dx=10
oz=0
ox=0
vc = 1500
dmax = 20

slow = SepVector.getSepVector(Hypercube.hypercube(ns=[nx,nz],os=[ox,oz],ds=[dx,dz]),storage='dataComplex')
slowNd = slow.getNdArray()
slow.set(vc)
num_scatter = int(nz*nx/4)

sign = -1
for i in range(num_scatter):
    x = np.random.randint(1,nx-1)
    z = np.random.randint(1,nz-1)
    dvel = sign*dmax*np.random.rand()
    slowNd[x,z] += dvel
    sign *= -1

num_gauss = 2
dev = [7,10]
perturb = [100,150]
center = [30,70]
x,z = np.meshgrid(np.arange(0,nz),np.arange(0,nz))
sign = -1
for i in range(num_gauss):
    mu = center[i]
    slowNd += sign * perturb[i] * np.exp((-(x-mu)**2 - (z-mu)**2)/(2*dev[i]**2))
    sign *= -1

genericIO.defaultIO.writeVector("Vel/random_true.H",slow)
