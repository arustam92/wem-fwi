import py_wem_ops as pyOps
import pyOperator as Op
import numpy as np


class SplitStep (Op.Operator):
	def __init__(self,model,data,slow,freq,loc,s0,depth):
		self.cppMode = pyOps.SplitStep(slow.cppMode)
		self.cppMode.setFreq(freq)
		self.cppMode.setLocation(loc)
		self.cppMode.setSlow(s0)
		self.cppMode.setDepth(depth)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class Scatter (Op.Operator):
	def __init__(self,model,data,slow,ntaylor,freq,depth,par):
		self.cppMode = pyOps.Scatter(slow.cppMode, ntaylor, par.cppMode)
		self.cppMode.setFreq(freq)
		self.cppMode.setDepth(depth)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class SSF (Op.Operator):
	def __init__(self,model,data,slow,freq,depth,par,ref):
		self.cppMode = pyOps.SSF(slow.cppMode, par.cppMode, ref.cppMode)
		self.cppMode.setFreq(freq)
		self.cppMode.setDepth(depth)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class Reflect (Op.Operator):
	def __init__(self,model,data,slow):
		mod = [m.cppMode for m in slow.vecs]
		self.cppMode = pyOps.Reflect(mod)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class RefSampler (Op.Operator):
	def __init__(self,slow,nref):
		self.cppMode = pyOps.RefSampler(slow.cppMode,nref)

	def getRefSlow(self, iz, iref):
		return self.cppMode.getRefSlow(iz,iref)

class Phshift (Op.Operator):
	def __init__(self,model,data,hyper,dz,freq,s0,iref):
		self.cppMode = pyOps.Phshift(dz,hyper.cppMode,iref)
		self.cppMode.setFreq(freq)
		self.cppMode.setSlow(s0)
		self.cppMode.setRef(iref)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class Down (Op.Operator):
	def __init__(self,model,data,slow,par,ref,freq):
		self.cppMode = pyOps.Down(slow.cppMode,par.cppMode,ref.cppMode)
		self.cppMode.setFreq(freq)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class Up (Op.Operator):
	def __init__(self,model,data,slow,par,ref,freq):
		self.cppMode = pyOps.Up(slow.cppMode,par.cppMode,ref.cppMode)
		self.cppMode.setFreq(freq)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class IC (Op.Operator):
	def __init__(self,model,data,bg,depth):
		self.cppMode = pyOps.IC(bg.cppMode)
		self.cppMode.setDepth(depth)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class LinDown (Op.Operator):
	def __init__(self,model,data,slow,par,bg,oneway,freq):
		self.cppMode = pyOps.LinDown(slow.cppMode,par.cppMode,bg.cppMode,oneway.cppMode)
		self.cppMode.setFreq(freq)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)
	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class dReflect (Op.Operator):
	def __init__(self,model,data,slow):
		mod = [m.cppMode for m in slow.vecs]
		self.cppMode = pyOps.dReflect(mod)
		self.setDomainRange(model,data)
	def forward(self,add,model,data):
		mod = [m.cppMode for m in model.vecs]
		self.cppMode.forward(mod,data.cppMode,add)
	def adjoint(self,add,model,data):
		mod = [m.cppMode for m in model.vecs]
		self.cppMode.adjoint(mod,data.cppMode,add)

class Injection (Op.Operator):
	def __init__(self,model,data,z,x,ng=0,tap=0):
		self.cppMode = pyOp.Injection(z,x,ng,tap)
		self.setDomainRange(model,data)

	def setStep(self,step):
		self.cppMode.setStep(step)

	def setShot(self,shot):
		self.cppMode.setShot(shot)

	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)