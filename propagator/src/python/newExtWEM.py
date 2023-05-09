import pyWEM
import WEM
import pyOperator as Op
import numpy as np
import SepVector
import Hypercube

class extWEM (WEM.Propagator):
	def __init__(self,model,data,wave,par):
		super().__init__(model, data, wave, par)
		self.cppMode = pyWEM.extWEM(wave.cppMode,model.vecs[0].getHyper().cppMode,par.cppMode)

	def wavefield(self,model,data):
		self.cppMode.wavefield(model.cppMode,data.cppMode,True)

	def setBgRefl(self,slow):
		self.cppMode.setBgRefl(slow.cppMode)

	def __setstate__(self, state):
		self.__dict__ = state.copy()
		self.cppMode = pyWEM.extWEM(self.wave.cppMode,self.domain[0].getHyper().cppMode, self.par.cppMode)


class extBorn (WEM.Propagator):
	def __init__(self,model,data,wave,par):
		super().__init__(model, data, wave, par)
		mod = [m.cppMode for m in model.vecs]
		self.cppMode = pyWEM.extBorn(wave.cppMode,mod,par.cppMode)

	def set_background(self,slow):
		print("setting bg")
		mod = [m.cppMode for m in slow.vecs]
		self.cppMode.setBgSlow(mod)

	def __setstate__(self, state):
		self.__dict__ = state.copy()
		mod = [m.cppMode for m in self.domain.vecs]
		self.cppMode = pyWEM.extBorn(self.wave.cppMode,mod, self.par.cppMode)
