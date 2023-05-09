import pyWEM
import pyOperator as Op
import genericIO
import numpy as np
from pyVector import superVector


class Propagator(Op.Operator):
	def __init__(self,model,data,wave,par):
		if not isinstance(model, superVector):
			raise TypeError("Model vector should be a superVector")
		self.wave = wave
		self.par = par
		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		mod = [m.cppMode for m in model]
		self.cppMode.forward(mod, data.cppMode,add)

	def adjoint(self,add,model,data):
		mod = [m.cppMode for m in model]
		self.cppMode.adjoint(mod,data.cppMode,add)

	@classmethod
	def from_subspace(cls, subdomain, subrange, wave, par, *args, **kw):
		pars = par.pars
		pars["ns"] = subrange.getHyper().getAxis(3).n
		pars["osx"] = subrange.getHyper().getAxis(3).o
		subpar = genericIO.pythonParams(pars)
		return cls(subdomain, subrange, wave, subpar, *args, **kw)
	
	def __getstate__(self):
		state = self.__dict__.copy()
		del state['cppMode']
		return state
	
class WEM (Propagator):
	def __init__(self,model,data,wave,par):
		super().__init__(model, data, wave, par)
		self.cppMode = pyWEM.WEM(wave.cppMode,model[0].getHyper().cppMode, par.cppMode)

	def setBgRefl(self,slow):
		self.cppMode.setBgRefl(slow.cppMode)

	def __setstate__(self, state):
		self.__dict__ = state.copy()
		self.cppMode = pyWEM.WEM(self.wave.cppMode,self.domain[0].getHyper().cppMode, self.par.cppMode)

	def wavefield(self,model,wavefield):
		self.cppMode.wavefield(model.getCpp(),wavefield.getCpp(),False)

class Born (Propagator):
	def __init__(self,model,data,wave,par,mode=None):
		super().__init__(model, data, wave, par)
		self.mode = mode
		mod = [m.cppMode for m in model]
		self.cppMode = pyWEM.Born(wave.cppMode,mod, par.cppMode)

	def forward(self,add,model,data):
		if self.mode == 'real':
			for m in model: m[:].imag = 0
		Propagator.forward(self, add, model, data)

	def adjoint(self,add,model,data):
		Propagator.adjoint(self, add, model, data)
		if self.mode == 'real':
			for m in model: m[:].imag = 0
	
	@classmethod
	def from_subspace(cls, subdomain, subrange, wave, par, mode=None):
		return super().from_subspace(subdomain, subrange, wave, par, mode=mode)

	def set_background(self,slow):
		mod = [m.cppMode for m in slow]
		self.cppMode.setBgSlow(mod)

	def __setstate__(self, state):
		self.__dict__ = state.copy()
		mod = [m.cppMode for m in self.domain]
		self.cppMode = pyWEM.Born(self.wave.cppMode, mod, self.par.cppMode)




class BornRefl (Born):
	def __init__(self,model,slowness,data,wave,par):
		self.cppMode = pyWEM.BornRefl(wave.cppMode,slowness.cppMode,par.cppMode)
		self.setDomainRange(model,data)
		self.prec = None
		self.pad = None

	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

	def setBgRefl(self,slow):
		self.cppMode.setBgRefl(slow.cppMode)

	def add_prec(self,op):
		self.prec = op
		self.big_model = op.range.clone()
	def add_pad(self,op):
		self.pad = op
		self.pad_model = op.range.clone()

	def set_background(self,slow):
		if (self.prec is not None):
			self.prec.forward(False,slow,self.big_model)
			slow = self.big_model
		if (self.pad is not None):
			self.pad.forward(False,slow,self.pad_model)
			slow = self.pad_model
		self.cppMode.setBgSlow(slow.cppMode)

class BornTomo (Born):
	def __init__(self,model,data,wave,par):
		self.cppMode = pyWEM.BornTomo(wave.cppMode,model.cppMode,par.cppMode)
		self.setDomainRange(model,data)
		self.pad = None

	def add_pad(self,op):
		self.pad = op
		self.pad_refl = op.range.clone()

	def setBgRefl(self,refl):
		if (self.pad is not None):
			self.pad.forward(False,refl,self.pad_refl)
			refl = self.pad_refl
		self.cppMode.setBgRefl(refl.cppMode)

	def set_background(self,slow):
		self.cppMode.setBgSlow(slow.cppMode)
