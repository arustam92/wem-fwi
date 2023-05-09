import SepVector
import Hypercube
import genericIO
import WEM
import json
import scipy.io

# modeling
nz = 512
nx = 1024
dz = 0.5
dx = 0.5
oz=0
ox=0

vel = SepVector.getSepVector(Hypercube.hypercube(ns=[nx,nz],
													os=[ox,oz],
													ds=[dx,dz]))
velNd = vel.getNdArray()

velNd[:] = scipy.io.loadmat("VelocityMaps200717.mat")['bodyTgt']

genericIO.defaultIO.writeVector("Vel/bodyTrue.H", vel)
