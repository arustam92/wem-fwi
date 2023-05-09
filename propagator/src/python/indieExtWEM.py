import WEM
import genericIO
import SepVector
import Hypercube
import pyOperator as Op
import numpy as np
import time
import Operator

class ExtWEM (Op.Operator):
	def __init__(self,model,data,par,wave,preOp=None):
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
		for i in range(len(self.wemOp)):
			mod1fNd = self.model[i].getNdArray()
			mod1fNd[:] = modNd[i,:,:,:]
			self.wemOp[i].setBgRefl(self.model[i])
			self.wemOp[i].forward(False,self.model[i],self.data[i])
			datNd[:,:,i] += np.squeeze(self.data[i].getNdArray())

		print('%s s' % (time.time()-t0))

class ExtBorn(Op.Operator):
	def __init__(self,model,data,wave,par,preOp=None):
		"""
			model -- [nx,nz,nw]
			data -- [nw,nr,ns] in frequency domain (complex)
		"""

		self.born = []
		# self.bornTomo = []
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

			born = WEM.Born(self.bg[i],self.data[i],wave1f,par)
			self.born.append(born)
			# bornTomo = WEM.BornTomo(self.bg[i],self.data[i],wave1f,par,ext=True)
			# self.bornTomo.append(bornTomo)

		self.setDomainRange(model,data)


	def forward(self,add,model,data):
		print('Extended fwd Born...')
		if (not add): data.scale(0)

		t0 = time.time()

		datNd = data.getNdArray()
		modNd = model.getNdArray()
		for i in range(len(self.born)):
			mod1fNd = self.model[i].getNdArray()
			mod1fNd[:] = modNd[i,:,:,:]
			self.born[i].forward(False,self.model[i],self.data[i])
			datNd[:,:,i] += np.squeeze(self.data[i].getNdArray())
			# self.bornTomo[i].forward(False,self.model[i],self.data[i])
			# datNd[:,:,i] += np.squeeze(self.data[i].getNdArray())
		print('%s s' % (time.time()-t0))



	def adjoint(self,add,model,data):
		print('Extended adj Born...')
		if (not add): model.scale(0)

		t0 = time.time()
		datNd = data.getNdArray()
		modNd = model.getNdArray()

		for i in range(len(self.born)):
			dat1fNd = self.data[i].getNdArray()
			dat1fNd[:,:,0] = np.ascontiguousarray(datNd[:,:,i])
			self.born[i].adjoint(False,self.model[i],self.data[i])
			modNd[i,:,:] += (self.model[i].getNdArray())
			# self.bornTomo[i].adjoint(False,self.model[i],self.data[i])
			# modNd[i,:,:] += (self.model[i].getNdArray()).real
		# modNd[:,0,:,:] += np.conj(modNd[:,1,:,:])
		# modNd[:,1,:,:] = np.conj(modNd[:,0,:,:])

		print('%s s' % (time.time()-t0))



	def setBgSlow(self,slow):
		print('setting background')
		slowNd = slow.getNdArray()
		for i in range(len(self.born)):
			bgNd = self.bg[i].getNdArray()
			bgNd[:] = slowNd[i,:,:,:]
			self.born[i].setBgSlow(self.bg[i])
			# bgNd[:] = slowNd[i,:,:].real
			# self.bornTomo[i].setBgSlow(self.bg[i])
