import numpy as np
import holoview as hv
hv.extension('matplotlib')

class HoloMap:

  def __init__(self, vec, axis=0, skip=1,  *args, **kwargs):
    
    ax = vec.getHyper().axes
    left = 
    d = {}
    for i in range(0, arr.shape[axis], skip):
      image = hv.Image(np.take(vec[:], i, axis=axis), bounds=bounds)
