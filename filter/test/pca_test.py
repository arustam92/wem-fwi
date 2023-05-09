import sys
import numpy as np

sys.path.append('../python')
import PCA

A = np.random.rand(3,10)

out = PCA.PCA(A)
