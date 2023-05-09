import subprocess
import sys
import re

for par in sys.argv[1:]:
    if ('in=' in par):
        input = par.split('in=')[1]
    if ('tagout=' in par):
        out = par.split('tagout=')[1]
    if ('bands=' in par):
        bands = []
        bands = par.split('bands=')[1].split(',')

for band in bands:
    fmin = band.split('-')[0]
    fmax = band.split('-')[1]
    subprocess.check_output("Window3d < %s min1=%s max1=%s > %s_%s.H" % (input,fmin,fmax,out,band), shell=True)
