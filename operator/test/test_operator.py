import genericIO
import SepVector
import Hypercube
import unittest
import Operator
import pyOperator as Op
import numpy as np

class TestSpline3D(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200,100,50]
    os = [0,0,1]
    ds = [1,1,1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    self.op = Operator.Spline3D(model,data,type='CR-spline')

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestFFT_out(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200,100,50]
    os = [0,0,1]
    ds = [1,1,1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    self.op = Operator.FFT(model,data,1,1,mode="out")

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestTaper(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200,100,50]
    os = [0,0,1]
    ds = [1,1,1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    self.op = Operator.Taper(model,data,tap=[10],axes=[0],zero=[0])

  def test_dot_product(self):
    self.op.dotTest(verbose=True)


class TestLanczosInterpolation3d(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200,100,50]
    os = [0,0,1]
    ds = [1,1,1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    self.op = Operator.LanczosInterpolation3D(model,data,a=[1,1,3],taper=[0,0,0.4])

  def test_dot_product(self):
    self.op.dotTest(verbose=True)

class TestHyperbolicPenalty(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200,100,50]
    os = [0,0,1]
    ds = [1,1,1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    nlOp = Operator.HyperbolicPenalty(model,data,l=1,tau=1e-12)
    linOp = Operator.Softclip(model,data,l=1,tau=1e-12)
    self.op = Op.NonLinearOperator(nlOp,linOp,set_background_func=linOp.setBg)
    bg = model.clone().set(1+1j)
    self.op.lin_op.setBg(bg)

  def test_dot_product(self):
    self.op.lin_op.dotTest(verbose=True)

  def test_linearization(self):
    pass

class TestDerivative(unittest.TestCase):
  def setUp(self):
    print('.............. %s ................' % self.__class__.__name__)    
    ns = [200,100,50]
    os = [0,0,1]
    ds = [1,1,1]
    model = SepVector.getSepVector(Hypercube.hypercube(ns=ns,ds=ds,os=os),storage='dataComplex')
    data = model.clone()
    self.op = Operator.Derivative(model,data,which=1,order='exact',f0=1,dt=0.004,tmax=6,alpha=0,beta=1,mode='constant')

  def test_dot_product(self):
    self.op.dotTest(verbose=True)


if __name__ == '__main__':
    unittest.main()