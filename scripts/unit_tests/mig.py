import SepVector
import Hypercube
import genericIO
import ExtWEM
import json
import Operator

par = {
    "ns" : 20,
    "osx" : 0,
    "osz" : 0,
    "dsx" : 50,

    "nr" : 100,
    "orx" : 0,
    "orz" : 0,
    "drx" : 10,

    "prop" : 'ssf',
    "nref" : 1,
    "tap" : 10,
    "fmin" : 1,
    "fmax" : 15,
    "nfreq" : 1,
    "ng" : 1,

    "illum" : 1
}
parObj = genericIO.pythonParams(par).getCpp()

slow = genericIO.defaultIO.getVector("extBg.H")
wave = genericIO.defaultIO.getVector("c_wave_1-15.H")
data = genericIO.defaultIO.getVector("c_unit_test_spike_data_1-15.H")

# create a WEM propagator
rect=[2,5]
repeat=1
smoothOp = Operator.Smooth(slow,slow,rect,repeat)
bornOp = ExtWEM.ExtBorn(slow,data,wave,parObj)
wemOp = ExtWEM.ExtWEM(slow,data,parObj,wave)

image = slow.clone()
image.zero()
bornOp.adjoint(False,image,data)
genericIO.defaultIO.writeVector("mig.H" ,image)
