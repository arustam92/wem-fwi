import WEM
import genericIO
import SepVector
import Hypercube
import pyOperator as Op
import numpy as np

class ExtWEM (Op.Operator):
	def __init__(self,model,data,par,wave):
		"""
			model -- [nx,nz,nw]
			data -- [nw,nr,ns] in frequency domain (complex)
		"""
		self.wemOp = []
		self.data = []
		self.model = []
		ax3 = Hypercube.axis(n=1,d=1,o=0)	# just to make data 3d
		dw = wave.getHyper().getAxis(1).d
		o = wave.getHyper().getAxis(1).o
		for i in range(model.getNdArray().shape[0]):
			self.data.append(SepVector.getSepVector(axes=[ax3,data.getHyper().getAxis(2),data.getHyper().getAxis(3)],storage='dataComplex'))
			self.model.append(SepVector.getSepVector(axes=[model.getHyper().getAxis(1),model.getHyper().getAxis(2)]))
			ow = o + i*dw
			wave1f = SepVector.getSepVector(ns=[1],ds=[dw],os=[ow],storage='dataComplex')
			wavNd = wave1f.getNdArray()
			wavNd[:] = wave.getNdArray()[i]
			# modNd = self.model[i].getNdArray()
			# datNd = self.data[i].getNdArray()

			# modNd[:,:] = model.getNdArray()[i,:,:]	# copy to all the slices

			# datNd = data.getNdArray()[i,:,:]	# copy to all the slices
			# modNd = 0
			# print(model.getNdArray()[i,:,:])

			wem = WEM.WEM(self.model[i],self.data[i],par,wave1f)
			# wem.setFminFmax(0, 1)	# because we dont need zero frequency
			self.wemOp.append(wem)

		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		# prange?
		datNd = data.getNdArray()
		modNd = model.getNdArray()
		for i in range(len(self.wemOp)):
			mod1fNd = self.model[i].getNdArray()
			mod1fNd[:] = modNd[i,:,:]
			self.wemOp[i].forward(add,self.model[i],self.data[i])
			datNd[:,:,i] = np.squeeze(self.data[i].getNdArray())

class ExtBorn(Op.Operator):
	def __init__(self,model,data,wave,par):
		"""
			model -- [nx,nz,nw]
			data -- [nw,nr,ns] in frequency domain (complex)
		"""
		self.bornOp = []
		self.data = []
		self.bg = []
		self.model = []
		ax3 = Hypercube.axis(n=1,d=1,o=0)	# just to make data 3d
		dw = wave.getHyper().getAxis(1).d
		o = wave.getHyper().getAxis(1).o
		for i in range(model.getNdArray().shape[0]):
			self.data.append(SepVector.getSepVector(axes=[ax3,data.getHyper().getAxis(2),data.getHyper().getAxis(3)],storage='dataComplex'))
			self.bg.append(SepVector.getSepVector(axes=[model.getHyper().getAxis(1),model.getHyper().getAxis(2)]))
			self.model.append(SepVector.getSepVector(axes=[model.getHyper().getAxis(1),model.getHyper().getAxis(2)]))
			ow = o + i*dw
			wave1f = SepVector.getSepVector(ns=[1],ds=[dw],os=[ow],storage='dataComplex')
			wavNd = wave1f.getNdArray()
			wavNd[:] = wave.getNdArray()[i]

			modNd = self.bg[i].getNdArray()
			modNd[:,:] = model.getNdArray()[i,:,:]

			born = WEM.Born(self.bg[i],self.data[i],wave1f,par)
			# born.setBgSlow(self.model[i])
			self.bornOp.append(born)

		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		datNd = data.getNdArray()
		modNd = model.getNdArray()
		for i in range(len(self.bornOp)):
			mod1fNd = self.model[i].getNdArray()
			mod1fNd[:] = modNd[i,:,:]
			self.bornOp[i].forward(add,self.model[i],self.data[i])
			datNd[:,:,i] = np.squeeze(self.data[i].getNdArray())

	def adjoint(self,add,model,data):
		datNd = data.getNdArray()
		modNd = model.getNdArray()
		for i in range(len(self.bornOp)):
			dat1fNd = self.data[i].getNdArray()
			dat1fNd[:,:,0] = np.ascontiguousarray(datNd[:,:,i])
			self.bornOp[i].adjoint(add,self.model[i],self.data[i])
			modNd[i,:,:] = (self.model[i].getNdArray())


	def setBgSlow(self,slow):
		slowNd = slow.getNdArray()
		for i in range(len(self.bornOp)):
			bgNd = self.bg[i].getNdArray()
			bgNd[:,:] = slowNd[i,:,:]
			self.bornOp[i].setBgSlow(self.bg[i])
