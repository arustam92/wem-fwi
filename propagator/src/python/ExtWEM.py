import WEM
import genericIO
import SepVector
import Hypercube
import pyOperator as Op
import numpy as np
import time
import Operator

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
			self.model.append(SepVector.getSepVector(axes=[model.getHyper().getAxis(1),model.getHyper().getAxis(2),model.getHyper().getAxis(3)],storage='dataComplex'))
			ow = o + i*dw
			wave1f = SepVector.getSepVector(ns=[1],ds=[dw],os=[ow],storage='dataComplex')
			wavNd = wave1f.getNdArray()
			wavNd[:] = wave.getNdArray()[i]

			wem = WEM.WEM(self.model[i],self.data[i],par,wave1f,ext=True)
			self.wemOp.append(wem)

		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		# prange?
		print('Extended WEM...')
		if (not add): data.scale(0)
		t0 = time.time()
		datNd = data.getNdArray()
		modNd = model.getNdArray()
		temp = np.zeros((model.getNdArray().shape[1],model.getNdArray().shape[2],model.getNdArray().shape[3]))
		for i in range(len(self.wemOp)):
			mod1fNd = self.model[i].getNdArray()
			mod1fNd[:] = modNd[i,:,:,:]
			temp += mod1fNd[:].real / len(self.wemOp)

			mod1fNd = self.model[i].getNdArray()
			mod1fNd[:] = temp[:]
			self.wemOp[i].forward(False,self.model[i],self.data[i])
			datNd[:,:,i] += np.squeeze(self.data[i].getNdArray())
			# self.wemOp[i].setBgRefl(self.model[i])
			# self.wemOp[i].forward(add,self.model[i],self.data[i])
			# datNd[:,:,i] += .5 * np.squeeze(self.data[i].getNdArray())

		print('%s s' % (time.time()-t0))

class ExtBorn(Op.Operator):
	def __init__(self,model,data,wave,par):
		"""
			model -- [nx,nz,nw]
			data -- [nw,nr,ns] in frequency domain (complex)
		"""
		self.bornRefl = []
		self.bornTomo = []
		self.data = []
		self.bg = []
		self.model = []
		self.bgRefl = []
		ax3 = Hypercube.axis(n=1,d=1,o=0)	# just to make data 3d
		dw = wave.getHyper().getAxis(1).d
		o = wave.getHyper().getAxis(1).o
		for i in range(model.getNdArray().shape[0]):
			self.data.append(SepVector.getSepVector(axes=[ax3,data.getHyper().getAxis(2),data.getHyper().getAxis(3)],storage='dataComplex'))
			self.bg.append(SepVector.getSepVector(axes=[model.getHyper().getAxis(1),model.getHyper().getAxis(2),model.getHyper().getAxis(3)],storage='dataComplex'))
			self.model.append(SepVector.getSepVector(axes=[model.getHyper().getAxis(1),model.getHyper().getAxis(2),model.getHyper().getAxis(3)],storage='dataComplex'))
			ow = o + i*dw
			wave1f = SepVector.getSepVector(ns=[1],ds=[dw],os=[ow],storage='dataComplex')
			wavNd = wave1f.getNdArray()
			wavNd[:] = wave.getNdArray()[i]

			modNd = self.bg[i].getNdArray()
			modNd[:,:] = model.getNdArray()[i,:,:]

			bornRefl = WEM.BornRefl(self.bg[i],self.data[i],wave1f,par,ext=True)
			bornTomo = WEM.BornTomo(self.bg[i],self.data[i],wave1f,par,ext=True)
			self.bornRefl.append(bornRefl)
			self.bornTomo.append(bornTomo)

		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		print('Extended fwd Born...')
		if (not add): data.scale(0)

		t0 = time.time()
		datNd = data.getNdArray()
		modNd = model.getNdArray()
		temp = np.zeros((model.getNdArray().shape[1],model.getNdArray().shape[2],model.getNdArray().shape[3]))
		for i in range(len(self.bornRefl)):
			mod1fNd = self.model[i].getNdArray()
			mod1fNd[:] = modNd[i,:,:,:]
			temp += mod1fNd[:].real / len(self.bornRefl)
			self.bornRefl[i].forward(False,self.model[i],self.data[i])
			datNd[:,:,i] += np.squeeze(self.data[i].getNdArray())
		for i in range(len(self.bornTomo)):
			mod1fNd = self.model[i].getNdArray()
			mod1fNd[:] = temp[:]
			self.bornTomo[i].forward(False,self.model[i],self.data[i])
			datNd[:,:,i] += np.squeeze(self.data[i].getNdArray())
		print('%s s' % (time.time()-t0))



	def adjoint(self,add,model,data):
		print('Extended adj Born...')
		if (not add): model.scale(0)

		t0 = time.time()
		datNd = data.getNdArray()
		modNd = model.getNdArray()
		temp = np.zeros((model.getNdArray().shape[1],model.getNdArray().shape[2],model.getNdArray().shape[3]))

		for i in range(len(self.bornRefl)):
			dat1fNd = self.data[i].getNdArray()
			dat1fNd[:,:,0] = np.ascontiguousarray(datNd[:,:,i])
			self.bornRefl[i].adjoint(False,self.model[i],self.data[i])
			modNd[i,:,:,:] += (self.model[i].getNdArray())

		for i in range(len(self.bornTomo)):
			dat1fNd = self.data[i].getNdArray()
			dat1fNd[:,:,0] = np.ascontiguousarray(datNd[:,:,i])
			self.bornTomo[i].adjoint(False,self.model[i],self.data[i])
			temp += self.model[i].getNdArray().real
		# Spread Tomo over all frequencies
		modNd[:] += temp
		print('%s s' % (time.time()-t0))



	def setBgSlow(self,slow):
		print('setting background')
		temp = np.zeros((slow.getNdArray().shape[1],slow.getNdArray().shape[2],slow.getNdArray().shape[3]))
		slowNd = slow.getNdArray()
		for i in range(len(self.bornRefl)):
			bgNd = self.bg[i].getNdArray()
			bgNd[:] = slowNd[i,:,:]
			temp += bgNd[:].real / len(self.bornRefl)
			self.bornRefl[i].setBgRefl(self.bg[i])
			self.bornTomo[i].setBgRefl(self.bg[i])
		for i in range(len(self.bornTomo)):
			bgNd = self.bg[i].getNdArray()
			bgNd[:] = temp[:]
			self.bornRefl[i].setBgSlow(self.bg[i])
			self.bornTomo[i].setBgSlow(self.bg[i])
