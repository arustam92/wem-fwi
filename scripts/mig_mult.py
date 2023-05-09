import os
import sys

wave = "/home/user/project//projects/marmousiWEM/Wav/c_MarmWave_3-40.H"
data = "/home/user/project//projects/marmousiWEM/Dat/c_datWEM_3-40.H"
geom = "/home/user/project//projects/marmousiWEM/Par/geom.json"

realizations = 100
for i in range(realizations):
    print("Realization %s" % i)
    vel = str("Vel2/rand_%s.H" % i)
    out = str("Image2/mig_rand_%s.H" % i)
    os.system('python scripts/mig.py vel=%s wave=%s geom=%s data=%s out=%s ext=1' % (vel,wave,geom,data,out))
    os.system('python scripts/f2tau.py "Image2/mig_rand_%s.H" "Image2/mig_rand_%s_tau.H"' % (i,i))
