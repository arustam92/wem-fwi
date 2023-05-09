import pyWEM
import pyOperator as Op
import numpy as np


class WEM (Op.Operator):
	def __init__(self,model,data,par,wave,src=None,rec=None,linear=False,verbose=False):
		message = 'Create propagator:'
		if src!=None and rec!=None:
			message += 'irregular,'
			if linear==False:
				message += 'nonlinear'
				self.cppMode = pyWEM.WEM(wave.cppMode,model.cppMode,par,src.cppMode,rec.cppMode)
				self.setDomainRange(model,data)
			elif linear==True:
				message += 'linear'
				self.cppMode = pyWEM.WEM(model.cppMode,par,src.cppMode,rec.cppMode)
				self.setDomainRange(wave,data)
		elif linear==False:
			message += 'reguular, nonlinear'
			self.cppMode = pyWEM.WEM(wave.cppMode,model.cppMode,par)
			self.setDomainRange(model,data)
		elif linear==True:
			message += 'reguular, linear'
			self.cppMode = pyWEM.WEM(model.cppMode,par)
			self.setDomainRange(wave,data)
		if(verbose): print(message)

	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

	def setFminFmax(self,fmin,fmax):
		self.cppMode.setFminFmax(fmin,fmax)

class Born (Op.Operator):
	def __init__(self,model,data,wave,par,src=None,rec=None):
		if src!=None and rec!=None:
			self.cppMode = pyWEM.Born(wave.cppMode,model.cppMode,par,src.cppMode,rec.cppMode)
		else:
			self.cppMode = pyWEM.Born(wave.cppMode,model.cppMode,par)
		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

	def setBgSlow(self,slow):
		self.cppMode.setBgSlow(slow.cppMode)
