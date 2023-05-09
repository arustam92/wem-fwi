import SepVector
import Hypercube
import genericIO
import Operator
import numpy as np
from scipy.ndimage import gaussian_filter

size0 = 141
size = int(size0/2)+1
overlap = 30

f1 = 45
n1 = 400

true = genericIO.defaultIO.getVector("Image2/ext_MarmMig_3-100_tau.H")
trueNd = np.transpose(true.getNdArray(),(2,0,1))

realizations = 100
npatch = int(1 + ((n1-1) - size0)/overlap)
total = 2 * npatch * realizations # 2 because of flipping

nt = int(trueNd.shape[0]/2)+1
dataset = np.zeros((2,total,nt,size,size))

for i in range(realizations):
    print("Realization %s" % i)
    image = genericIO.defaultIO.getVector("Image2/mig_rand_%s_tau.H" % i)
    imageNd = np.transpose(image.getNdArray(),(2,0,1))
    for j in range(npatch):
        ox = j*overlap + f1
        endx = ox + size0
        if endx < (n1 + f1):
            im = imageNd[::2,::2,ox:endx:2]
            tr = trueNd[::2,::2,ox:endx:2]
            print(npatch*i + j)
            ind = 2*(npatch*i + j)
            dataset[0,ind,:,:,:] = im
            dataset[1,ind,:,:,:] = tr
            dataset[0,ind + 1,:,:,:] = im[:,:,::-1]
            dataset[1,ind + 1,:,:,:] = tr[:,:,::-1]

np.save("jupyter/Patch/train_set_extMarm-30over.npy",dataset)
