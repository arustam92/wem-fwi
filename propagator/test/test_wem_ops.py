import genericIO
import SepVector
import Hypercube
import unittest
import Operator
import wem_ops as wem
import pyOperator as Op
import numpy as np
import WEM
import newExtWEM
from pyVector import superVector

ns = [200,100]
os = [0,0]
ds = [10,10]
onepass = False

if onepass:
  orz = 900
else:
  orz = 0

par = {
    "ns" : 10,
    "osx" : 1000,
    "dsx" : ds[0],
    "osz" : 50,

    "nr" : 101,
    "orx" : 0,
    "drx" : ds[0],
    "orz" : orz,

    "eps" : 0.01,

    "prop" : "ssf",
    "nref" : 1,
    "tap" : 100,
    "pad" : 100,
    "fmin" : 1,
    "fmax" : 0,
    "nfreq" : 10,
    "ngs" : 3,
    "ngr" : 3,
    "ntaylor" : 1,

    "onepass" : onepass,
}
parObj = genericIO.pythonParams(par)
slow = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
slow.set(1/1000**2 - 1e-8j)
slow.getNdArray()[50:,:] = 1/2000**2 - 1e-8j

den = slow.clone()
den.set(1)

model2d = superVector(slow, den)

ref = wem.RefSampler(slow,par["nref"])
nt = 1001
ot = 0
dt = 0.004
t = np.linspace(ot,(nt-1)*dt,nt) - .5
f0 = 5

nf = int((par["fmax"]-par["fmin"]) * (nt-1) * dt + 1)
if nf < 0: 
  nf = par["nfreq"]
ns3 = ns
ds3 = ds
os3 = os
ns3.append(nf)
os3.append(1)
ds3.append(1) 
print(ns3)
slow3d = SepVector.getSepVector(Hypercube.hypercube(ns=ns3,ds=ds3,os=os3),storage='dataComplex')
slow3d.getNdArray()[:] = slow.getNdArray()[:]

den3d = slow3d.clone()
den3d.set(1)
model3d = superVector(slow3d, den3d)

wave = SepVector.getSepVector(Hypercube.hypercube(ns=[nt],ds=[dt],os=[ot]))
wave.getNdArray()[:] = (1-2*(np.pi*f0*t)**2)*np.exp(-((np.pi*f0*t)**2))

class TestInjection1D_2D(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200,100]
    os = [0,0]
    ds = [1,1]
    data = SepVector.getSepVector(ns=ns,ds=ds,os=os,storage='dataComplex')
    model = SepVector.getSepVector(ns=[ns[1]],ds=[ds[1]],os=[os[1]],storage='dataComplex')
    x = np.linspace(0,ns[0]/2,ns[0]).astype(int).tolist()
    z = np.linspace(0,0,ns[0]).astype(int).tolist()
    self.op = Operator.Injection(model,data,z,x,ng=1,tap=1)
    self.op.setStep(20)

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestInjection2D_3D(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [100,100,50]
    os = [0,0,0]
    ds = [1,1,1]

    dns = [200,200]
    dos = [0,0]
    dds = [1,1]
    model = SepVector.getSepVector(ns=ns,ds=ds,os=os,storage='dataComplex')
    data = SepVector.getSepVector(ns=dns,ds=dds,os=dos,storage='dataComplex')
    x = np.linspace(0,dns[0]-1,ns[1]).astype(int).tolist()
    z = np.linspace(0,0,ns[1]).astype(int).tolist()
    self.op = Operator.Injection(model,data,z,x,ng=0,tap=0)
    self.op.setStep(20)
    self.op.setShot(20)

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestSSF(unittest.TestCase):
  def setUp(self):
    print('SSF')
    model = SepVector.getSepVector(Hypercube.hypercube(ns=[ns[0]],ds=[ds[0]],os=[os[0]]),storage='dataComplex')
    data = model.clone()
    self.op = wem.SSF(model,data,slow,5,50,parObj,ref)

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestPhShift(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)
    ns = [200]
    os = [0]
    ds = [1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    self.op = wem.Phshift(model,data,model.getHyper(),ds[0],5,1/1000,1)

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestSplitStep(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200]
    os = [0]
    ds = [1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    loc = np.random.randint(0,200,size=50).tolist()
    self.op = wem.SplitStep(model,data,slow,5,loc,1/1000**2,50)

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestScatter(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    model = SepVector.getSepVector(Hypercube.hypercube(ns=[ns[0]],ds=[ds[0]],os=[os[0]]),storage='dataComplex')
    data = model.clone()
    self.op = wem.Scatter(model,data,slow,par["ntaylor"],1,50,parObj)
  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestDown(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    model = slow.clone()
    data = model.clone()
    self.op = wem.Down(model,data,slow,parObj,ref,1)

  # def test_ref(self):
    # self.assertEqual(complex(ref.getRefSlow(10,0)), 1/1000**2 - 1e-10j)

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestUp(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    model = slow.clone()
    data = model.clone()
    self.op = wem.Up(model,data,slow,parObj,ref,1)

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestDReflect(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    model = slow.clone()
    data = model.clone()
    self.op = wem.dReflect(model,data,slow)
  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestIC(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200]
    os = [0]
    ds = [1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    bg = slow.clone().rand()
    self.op = wem.IC(model,data,bg,10)
  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestIC_2(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200,100]
    os = [0,0]
    ds = [1,1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    self.op = wem.IC(model,data,slow,10)
  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestLinDown(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    model = slow.clone()
    data = model.clone()
    bg = model.clone()
    bg.set(1.+1j)
    down = wem.Down(model,data,slow,parObj,ref,1)
    self.op = wem.LinDown(model,data,slow,parObj,bg,down,1)
  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestReflect(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    model = slow.clone()
    data = slow.clone()
    self.op = wem.Reflect(model,data,model2d)
  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestDReflect(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    model = model2d.clone()
    data = slow.clone()
    self.op = wem.dReflect(model,data,model2d)
  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestBorn(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    model = slow.clone()
    data = SepVector.getSepVector(ns=[nt//2,par["nr"],par['ns']],
                          ds=[1,par["drx"],par['dsx']],
                          os=[1,par["orx"],par['osx']],storage='dataComplex')
    self.op = WEM.Born(model2d,data,wave,parObj)
    self.op.set_background(model2d)
  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestExtWEM(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    self.data = SepVector.getSepVector(ns=[nt//2,par["nr"],par['ns']],
                          ds=[1,par["drx"],par['dsx']],
                          os=[1,par["orx"],par['osx']],storage='dataComplex')
    self.wem = WEM.WEM(model2d,self.data,wave,parObj)
    self.ewem = newExtWEM.extWEM(model3d,self.data,wave,parObj)

  def test_equal(self):
    self.wem.forward(False, model2d, self.data)
    d1 = self.data.getNdArray().copy()
    self.ewem.forward(False, model3d, self.data)
    d2 = self.data.getNdArray().copy()
    self.assertTrue((d1==d2).all())

class TestExtBorn(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)   
    model = model3d.clone()
    data = SepVector.getSepVector(
      Hypercube.hypercube(ns=[par["nr"],nt//2,par['ns']],
                          ds=[par["drx"],1,par['dsx']],
                          os=[par["orx"],0,par['osx']]),storage='dataComplex')
    self.op = newExtWEM.extBorn(model,data,wave,parObj)

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

if __name__ == '__main__':
    unittest.main()