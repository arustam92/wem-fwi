#! /usr/bin/python3

import sys
import getopt

import SepVector
import Hypercube
import genericIO
import numpy as np	

opts,args = getopt.getopt(sys.argv[1:],"n1:")
for par in args:
	if "n1" in par:
		n1 = int(par.split("n1=")[1])
	elif "o1" in par:
		o1 = float(par.split("o1=")[1])
	elif "d1" in par:
		d1 = float(par.split("d1=")[1])
	elif "n2" in par:
		n2 = int(par.split("n2=")[1])
	elif "o2" in par:
		o2 = float(par.split("o2=")[1])
	elif "d2" in par:
		d2 = float(par.split("d2=")[1])
	elif "vc" in par:
		vc = float(par.split("vc=")[1])
	elif "layers" in par:
		layers = int(par.split("layers=")[1])
	elif "vr" in par:
		vr = list(map(float,(par.split("vr=")[1].split(','))))
	elif "angles" in par:
		angles = list(map(float,(par.split("angles=")[1].split(','))))
	elif "out" in par:
		out = par.split("out=")[1]

######################################################
thick = n1/(layers+1)
bounds = np.zeros((layers+3,n2))

sigs = len(angles)
offset = np.linspace(n2/4,3*n2/4,sigs)


x1 = np.linspace(0,n2/2-1,n2/2)
x2 = np.linspace(n2/2,n2-1,n2/2)
# x1 = x1.astype(int,copy=False)
# x2 = x2.astype(int,copy=False)

k = 4*np.tan(np.pi*np.asarray(angles)/180)
# sig1 = thick/(1+np.exp(-k[0]*(x1-off1)))
# sig2 = thick/(1+np.exp(k[1]*(x2-off2)))
dx = thick*np.sin(np.pi*np.asarray(angles)/180)

x = np.zeros(int(n2/sigs))
bounds[layers+1][:] = layers*thick
bounds[layers+2][:] = n1-1
for i in range(layers):
	sign = -1
	for j in range(sigs):
		x = np.linspace(j*n2/sigs,(j+1)*n2/sigs-1,n2/sigs)
		x = x.astype(int,copy=False)
		bounds[i+1][x] = thick/(1+np.exp(sign*k[j]*(x-offset[j]-sign*i*dx[j]))) + i*thick
		sign *= -1
		# bounds[i+1][x2] = thick/(1+np.exp(k[1]*(x2-off2-i*dx[1]))) + i*thick

bounds = bounds.astype(int,copy=False)
#########################################################
vel = SepVector.getSepVector(Hypercube.hypercube(ns=[n1,n2],os=[o1,o2],ds=[d1,d2]))
vel.set(vc)
velNP = vel.getNdArray()
for i in range (1,layers+1):
	for ix in range(n2):
		velNP[ix][bounds[i][ix]:bounds[i+1][ix]] = vr[i-1]

genericIO.defaultIO.writeVector(out,vel)



