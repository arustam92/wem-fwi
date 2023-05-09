import pyOp
import pyOperator as Op
import numpy as np
import SepVector
import Hypercube
from scipy.signal import fftconvolve, correlate, convolve, convolve2d, hilbert
from scipy.ndimage import gaussian_filter, correlate1d, convolve1d
from scipy.linalg import toeplitz



class StackAndSpread(Op.Operator):

	def __init__(self, model, data, axes=(1,2)):
		self.setDomainRange(model,data)
		self.ax = axes
		self.N = 1
		for i in self.ax:
			self.N *= model.getNdArray().shape[i]

	def forward(self,add,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		if (not add): data.scale(0)
		if len(modNd.shape) > 2:
			for i in range(modNd.shape[0]):
				datNd[i,:] += np.sum(modNd[i,:],axis=None) / self.N
		else:
			datNd[:] += np.sum(modNd[:],axis=None) / self.N

	def adjoint(self,add,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		if (not add): model.scale(0)
		if len(modNd.shape) > 2:
			for i in range(modNd.shape[0]):
				modNd[i,:] += np.sum(datNd[i,:],axis=None) / self.N
		else:
			modNd[:] += np.sum(datNd[:],axis=None) / self.N

class ScaleOffset(Op.Operator):

	def __init__(self, model, data, par, scale=1):
		self.setDomainRange(model,data)
		self.scale = 1
		self.par = par

	def forward(self,add,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		if (not add): data.scale(0)
		rx = np.linspace(self.par["orx"],(self.par["nr"]-1)*self.par["drx"],self.par["nr"])
		for i in range(datNd.shape[0]):
			sx = self.par["osx"] + self.par["dsx"]*i
			h = np.abs(sx-rx)
			datNd[i,:,:] += np.apply_along_axis(lambda x: x*h*self.scale, 0, modNd[i,:,:])

	def adjoint(self,add,model,data):
		self.forward(add,data,model)

class Interpolation(Op.Operator):
	"""docstring for Interpolation."""

	def __init__(self, model, data):
		self.setDomainRange(model,data)
		self.ds = data.getHyper().getAxis(3).d
		self.os = data.getHyper().getAxis(3).o
		self.ns = data.getHyper().getAxis(3).n
		self.mat = np.zeros((self.ns,model.getHyper().getAxis(3).n))

	def forward(self,add,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		if (not add): data.scale(0)
		temp = np.transpose(modNd[:],(1,2,0))
		for i in range(datNd.shape[0]):
			datNd[i,:,:] += np.dot(temp,self.mat[i,:])
		del temp

	def adjoint(self,add,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		if (not add): model.scale(0)
		temp = np.transpose(datNd[:],(1,2,0))
		matT = np.transpose(self.mat)
		for i in range(modNd.shape[0]):
			modNd[i,:,:] += np.dot(temp,matT[i,:])
		del temp, matT

class NearestNeighbor(Interpolation):
	def __init__(self,model,data,points=[0], axis=0):
		self.setDomainRange(model,data)
		ds = data.getHyper().getAxis(4).d
		os = data.getHyper().getAxis(4).o
		ns = data.getHyper().getAxis(4).n
		self.mat = np.zeros((ns,model.getHyper().getAxis(4).n))
		last = 0
		for i in range(len(points)):
			ipt = int((points[i]-os)/ds)
			self.mat[last:ipt,i] = 1
			last = ipt
		self.mat[last:ns,len(points)-1] = 1
		print(self.mat)

class LinearInterpolation(Interpolation):
	def __init__(self,model,data,points=[0], axis=0):
		super().__init__(model,data)
		self.setDomainRange(model,data)
		last = 0
		icur = int((points[0]-self.os)/self.ds)
		self.mat[last:icur,0] = 1
		for i in range(len(points)-1):
			icur = int((points[i]-self.os)/self.ds + 0.5)
			inext = int((points[i+1]-self.os)/self.ds + 0.5)
			for j in range(icur,inext):
				cur_pt = self.os + j*self.ds
				alpha = max(0, 1 - (cur_pt - points[i])/(points[i+1]-points[i]))
				self.mat[j,i] = alpha
				self.mat[j,i+1] = 1 - alpha
			last = inext
		self.mat[last:self.ns,len(points)-1] = 1
		# self.cppMode = pyOp.LinearInterpolation3D(data.cppMode,points)

	# def forward(self,add,model,data):
	# 	self.cppMode.forward(model.cppMode,data.cppMode,add)
	#
	# def adjoint(self,add,model,data):
	# 	self.cppMode.adjoint(model.cppMode,data.cppMode,add)


class CosineInterpolation(Interpolation):
	def __init__(self,model,data,points=[0], axis=0):
		super().__init__(model,data)
		self.setDomainRange(model,data)
		last = 0
		icur = int((points[0]-self.os)/self.ds)
		self.mat[last:icur,0] = 1
		for i in range(len(points)-1):
			icur = int((points[i]-self.os)/self.ds + 0.5)
			inext = int((points[i+1]-self.os)/self.ds + 0.5)
			for j in range(icur,inext):
				cur_pt = self.os + j*self.ds
				# alpha = max(0, 1 - (cur_pt - points[i])/(points[i+1]-points[i]))
				alpha = .5 * (1 + np.cos(np.pi * (cur_pt - points[i])/(points[i+1]-points[i])))
				self.mat[j,i] = alpha
				self.mat[j,i+1] = 1 - alpha
			last = inext
		self.mat[last:self.ns,len(points)-1] = 1

class LeakyIntegrationInverse(Op.Operator):
	"""docstring for LeakyIntegration."""

	def __init__(self, model, data, alpha=0.9):
		self.mat = np.eye(model.getNdArray().shape[0])
		for i in range(1,self.mat.shape[0]):
			for j in range(i):
				self.mat[i,j] = alpha**(i-j)
		self.mat[:] = np.linalg.inv(self.mat)

	def forward(self,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		res = np.matmul(self.mat,np.reshape(modNd,(modNd.shape[0],modNd.shape[1]*modNd.shape[2]*modNd.shape[3])))
		datNd[:] = np.reshape(res,datNd.shape)

class LeakyIntegration(Op.Operator):
	"""docstring for LeakyIntegration."""

	def __init__(self, model, data, alpha=0.9, alpha_final=False):
		if (alpha_final == False):
			self.alpha = alpha * np.ones(model.getNdArray().shape[0])
		else:
			self.alpha = np.linspace(alpha,alpha_final,model.getNdArray().shape[0])
		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		if (not add): data.scale(0)
		nw = modNd.shape[0]
		t = 0
		for i in range(nw):
			datNd[i,:,:] += modNd[i,:,:] + self.alpha[i] * t
			t = datNd[i,:,:]

	def adjoint(self,add,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		if (not add): model.scale(0)
		t = 0
		for i in range(datNd.shape[0]-1,-1,-1):
			modNd[i,:,:] += datNd[i,:,:] + self.alpha[i] * t
			t = modNd[i,:,:]

class DoubleLeakyIntegration(LeakyIntegration):

	def __init__(self, model, data, alpha=0.9, alpha_final=False):
		super().__init__(model,data,alpha=alpha,alpha_final=alpha_final)


	def forward(self,add,model,data):
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		temp = data.clone()
		tempNd = temp.getNdArray()
		super().adjoint(False,temp,model)
		tempNd[1:tempNd.shape[0],:,:] = np.apply_along_axis(lambda m : m * (1 - self.alpha[1:tempNd.shape[0]]**2), axis=0,arr=tempNd[1:tempNd.shape[0],:,:])
		super().forward(add,temp,data)
		del temp

	def adjoint(self,add,model,data):
		self.forward(add,data,model)

class FFT (Op.Operator):
	def __init__(self,model,data,rank,axis,mode=None):
		if mode:
			self.cppMode = pyOp.FFT(model.cppMode,mode,rank,axis)
		else:
			self.cppMode = pyOp.FFT(model.cppMode,data.cppMode,rank,axis)
		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class Derivative(Op.Operator):
	def __init__(self,model,data,which=1,order=4,f0=None,dt=None, tmax=None,alpha=0,beta=None, mode='reflect'):
		if tmax is None:
			raise Exception("Must provide tmax!")
		if f0 is None:
			raise Exception("Must provide f0!")
		self.shift = 0
		self.alpha = alpha
		if beta is None:
			self.beta = 1
		else:
			self.beta = beta

		if which==1:
			if order == 1:
				self.filter = np.array([-1,1])
			if order == 2:
				self.filter = np.array([-1/2,0,1/2])
			if order == 4:
				self.filter = np.array([1/12,-2/3,0,2/3,-1/12])
			if order == 8:
				self.filter = np.array([1/280,-4/105,1/5,-4/5,0,4/5,-1/5,4/105,-1/280])
			if order == 'exact':
				if dt is None:
					raise Exception("Must provide dt!")
				sc = 1
				t = np.arange(-sc*tmax,sc*tmax,dt)
				w = 2*np.pi*np.linspace(-1/(2*dt),1/(2*dt),t.size)
				# flip as in forward we actually perform correlation and not convolution
				w = np.flip(w,axis=0)
				# tmax = sc*tmax
				f1 = (np.cos(w*tmax/2)/w) * tmax
				f2 = -2*(np.sin(w*tmax/2)/w**2)
				f1 = np.where(w==0,0,f1)
				f2 = np.where(w==0,0,f2)
				self.filter = (f1 + f2)[::2*sc]
		if which==2:
			if order == 2:
				self.filter = np.array([1,-2,1])
			if order == 4:
				self.filter = np.array([-1/12,4/3,-5/2,4/3,-1/12])
		self.setDomainRange(model,data)
		self.mode = mode

	def forward(self,add,model,data):
		if (not add): data.scale(0)
		for mod, dat in zip(model, data):
			modNd = mod.getNdArray()
			datNd = dat.getNdArray()
			sum = np.sum(modNd.imag,axis=0) / modNd.shape[0]
			datNd.real[:] += -self.beta*np.apply_along_axis(lambda x: correlate1d(x,self.filter,mode=self.mode,origin=self.shift),0,modNd.imag)
			datNd.imag[:] += self.beta*np.apply_along_axis(lambda x: correlate1d(x,self.filter,mode=self.mode,origin=self.shift),0,modNd.real)
			# penalizing the zero-time-lag of imaginary part
			datNd.imag[:] += self.alpha * sum

	def adjoint(self,add,model,data):
		if (not add): model.scale(0)
		for mod, dat in zip(model, data):
			modNd = mod.getNdArray()
			datNd = dat.getNdArray()
			sum = np.sum(datNd.imag,axis=0) / datNd.shape[0]
			modNd.real[:] += self.beta*np.apply_along_axis(lambda x: convolve1d(x,self.filter,mode=self.mode,origin=self.shift),0,datNd.imag)
			modNd.imag[:] += -self.beta*np.apply_along_axis(lambda x: convolve1d(x,self.filter,mode=self.mode,origin=self.shift),0,datNd.real)
			# penalizing the zero-time-lag of imaginary part
			modNd.imag[:] += self.alpha * sum


class Laplacian2D(Op.Operator):
	"""docstring for Laplacian2D."""

	def __init__(self, model, data):
		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		if (not add): data.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		for iw in range(1,modNd.shape[0]-1):
			for iz in range(1,modNd.shape[1]-1):
				datNd[iw,iz,:] += modNd[iw-1,iz,:] + modNd[iw+1,iz,:] + \
								modNd[iw,iz-1,:] + modNd[iw,iz+1,:] - 4*modNd[iw,iz,:]

	def adjoint(self,add,model,data):
		self.forward(add,data,model)



class Smooth(Op.Operator):
	"""docstring for Smooth"""
	def __init__(self,model,data,l,rep=1,mod='all'):
		self.setDomainRange(model,data)
		self.box = np.array(l)
		if mod == 'all':
			print('using this')
			self.func = self.smooth_all
		else:
			if len(model.shape) > 3: self.func = self.smooth3d
			else: self.func = self.smooth2d
		self.repeat = rep
		self.tap = self.box[2]

	def forward(self,add,model,data):
		if (not add): data.scale(0)
		modelN = model.getNdArray()
		dataN = data.getNdArray()
		self.func(modelN,dataN)

	def adjoint(self,add,model,data):
		if (not add): model.scale(0)
		modelN = model.getNdArray()
		dataN = data.getNdArray()
		self.func(dataN,modelN)

	def smooth2d(self,input,output):
		for i in range(self.repeat):
			sig = np.copy(self.box)
			output[:,self.tap::,:] += input[:,self.tap::,:]
			for j in range(self.tap):
				sig[2] = self.box[2] - j
				t = gaussian_filter(input[:,j:j+1,:].real,sigma=tuple(sig))
				output[:,j:j+1,:].real += t
				t = gaussian_filter(input[:,j:j+1,:].imag,sigma=tuple(sig))
				output[:,j:j+1,:].imag += t

	def smooth3d(self,input,output):
		for i in range(self.repeat):
			sig = np.copy(self.box)
			output[:,:,self.tap::,:] += input[:,:,self.tap::,:]
			for iw in range(output.shape[0]):
				for j in range(self.tap):
					sig[2] = self.box[2] - j
					t = gaussian_filter(input[iw,:,j:j+1,:].real,sigma=tuple(sig))
					output[iw,:,j:j+1,:].real += t
					t = gaussian_filter(input[iw,:,j:j+1,:].imag,sigma=tuple(sig))
					output[iw,:,j:j+1,:].imag += t


	def smooth_all(self,input,output):
		# for iw in range(input.shape[0]):
			# t = input
			# for i in range(self.repeat):
		t = gaussian_filter(input[:].real,sigma=tuple(self.box))
		output[:].real += t
		t = gaussian_filter(input[:].imag,sigma=tuple(self.box))
		output[:].imag += t

	def smooth_all_2d(self,input,output):
		t = input
		for i in range(self.repeat):
			output[:] = fftconvolve(t[:],self.box,mode='same')
			t = output

class Taper(Op.Operator):
	"""docstring for Taper."""

	def __init__(self, model, data, tap=[0], axes=[0], zero=[0], shift=[(0,0)]):
		self.setDomainRange(model,data)
		self.filter = []
		self.axes = axes
		for axis in range(len(axes)):
			ind = np.linspace(0,tap[axis],tap[axis])
			filter = np.ones(model.getNdArray().shape[axes[axis]])
			lcos = .5 * (1 - np.cos(np.pi/tap[axis] * ind))
			filter[zero[axis]+shift[axis][0]:tap[axis]+zero[axis]+shift[axis][0]] = lcos
			filter[filter.size-tap[axis]-zero[axis]-shift[axis][1]:filter.size-zero[axis]-shift[axis][1]] = lcos[::-1]
			filter[shift[axis][0]:zero[axis]+shift[axis][0]] = 0
			filter[filter.size-zero[axis]-shift[axis][1]:filter.size-shift[axis][1]] = 0
			self.filter.append(filter)

	def forward(self, add, model, data):
		if (not add): data.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		for axis in range(len(self.filter)):
			datNd[:] += np.apply_along_axis(lambda m: m * self.filter[axis], axis=self.axes[axis], arr=modNd)

	def adjoint(self, add, model, data):
		if (not add): model.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		for axis in range(len(self.filter)):
			modNd[:] += np.apply_along_axis(lambda m: m * self.filter[axis], axis=self.axes[axis], arr=datNd)



class LanczosInterpolation3D (Op.Operator):
	def __init__(self,model,data,a=[3,3,3],taper=[0,0,0]):
		self.a = a
		self.taper = taper
		self.cppMode = pyOp.LanczosInterpolation3D(model.cppMode,data.cppMode,a,taper)
		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

	def __getstate__(self):
		state = self.__dict__.copy()
		del state['cppMode']
		return state
	
	def __setstate__(self, state):
		self.__dict__ = state.copy()
		self.cppMode = pyOp.LanczosInterpolation3D(self.domain.cppMode, self.range.cppMode, self.a, self.taper)


class LanczosInterpolation2D (Op.Operator):
	def __init__(self,model,data,a=[3,3,3]):
		self.cppMode = pyOp.LanczosInterpolation2D(model.cppMode,data.cppMode,a)
		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class Spline (Op.Operator):
	def __init__(self,model,data,a=1,b=0,type=None):
		if type == 'CR-spline':
			self.a = 0.5
			self.b = 0
		if type == 'B-spline':
			self.a = 0
			self.b = 1
		if type == 'MN-spline':
			self.a = 1/3
			self.b = 1/3

		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

class Spline2D (Spline):
	def __init__(self,model,data,a=1,b=0,type=None):
		super().__init__(model,data,a,b,type=type)
		self.cppMode = pyOp.Spline2D(model.cppMode,data.cppMode,a,b)

	def __getstate__(self):
		state = self.__dict__.copy()
		del state['cppMode']
		return state
	
	def __setstate__(self, state):
		self.__dict__ = state.copy()
		self.cppMode = pyOp.Spline2D(self.domain.cppMode, self.range.cppMode, self.a, self.b)

class Spline3D (Spline):
	def __init__(self,model,data,a=1,b=0,type=None,taper=[0,0,0]):
		super().__init__(model,data,a,b,type=type)
		self.taper = taper
		self.cppMode = pyOp.Spline3D(model.cppMode,data.cppMode,a,b,taper)

	def __getstate__(self):
		state = self.__dict__.copy()
		del state['cppMode']
		return state
	
	def __setstate__(self, state):
		self.__dict__ = state.copy()
		self.cppMode = pyOp.Spline3D(self.domain.cppMode, self.range.cppMode, self.a, self.b, self.taper)

class Pad(Op.Operator):
	"""docstring for Pad"""
	def __init__(self, model, data, beg=[], end=[], mode='edge'):
		self.beg = beg
		self.end = end
		self.data = data
		self.padding = tuple(np.vstack([beg, end]).T)
		self.mode = mode
		self.setDomainRange(model, data)

	def get_padded(self):
		return self.data
		
	@classmethod
	def from_params(cls, model, beg, end, mode='edge'):
		# construct a data vector
		ax = model.getHyper().axes
		beg = np.flip(np.array(beg))
		end = np.flip(np.array(end))
		os = [-a.d * beg[i] for i,a in enumerate(ax)]
		ds = [a.d for a in ax]
		ns = [a.n + end[i] for i,a in enumerate(ax)]
		data = SepVector.getSepVector(ns=ns, ds=ds, os=os, storage=model.kw["storage"])
		return cls(model, data, beg=np.flip(beg), end=np.flip(end), mode=mode)

	def forward(self,add,model,data):
		if not add: data.scale(0)
		data[:] += np.pad(model[:], self.padding, mode=self.mode)

	def adjoint(self,add,model,data):
		if (not add): model.scale(0)
		# for now can only do along the slowest axis
		model[:,...] += data[self.beg[0]:self.end[0],...]
		model[0,...] += np.sum(data[:self.beg[0],...], axis=0)
		model[-1,...] += np.sum(data[-self.end[0]:,...], axis=0)


class Mask(Op.Operator):
	def __init__(self,model,data, mask):
		self.mask = mask
		self.setDomainRange(model,data)

	def forward(self,add,model,data):
		if (not add): data.scale(0)
		for m, d, mask in zip(model, data, self.mask):
			modNd = m.getNdArray()
			datNd = d.getNdArray()
			datNd[:] += modNd[:] * mask[:]

	def adjoint(self,add,model,data):
		self.forward(add,data,model)

class Hilbert(Op.Operator):
	"""docstring for Hilbert."""

	def __init__(self, model, data, npad=None, sign=-1):
		self.real = np.zeros(model.shape,dtype=np.float32)
		self.npad = npad
		self.sign = sign
		self.setDomainRange(model,data)

	def forward(self, add, model, data):
		if (not add): data.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		self.real[:] = modNd.real[:]
		temp = hilbert(self.real,N=self.npad,axis=0)[0:modNd.shape[0],:]
		datNd.real[:] += temp.real[:]
		datNd.imag[:] += self.sign*temp.imag[:]

	def adjoint(self, add, model, data):
		if (not add): model.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		self.real[:] = datNd.imag[:]
		temp = hilbert(self.real,N=self.npad,axis=0)[0:datNd.shape[0],:]
		modNd.real[:] += datNd.real[:] - self.sign*temp.imag[:]
		modNd.imag[:] += 0
		# modNd.imag[:] = 0

class Hilbert2(Op.Operator):
	"""docstring for Hilbert."""

	def __init__(self, model, data, npad=None):
		self.real = np.zeros(model.shape,dtype=np.float32)
		self.npad = npad
		self.setDomainRange(model,data)

	def forward(self, add, model, data):
		if (not add): data.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		abs = np.abs(modNd)
		self.real[:] = np.cos(np.arctan2(modNd.imag[:],modNd.real[:]))
		temp = hilbert(self.real,N=self.npad,axis=0)[0:modNd.shape[0],:]
		datNd.real[:] += modNd.real[:]
		datNd.imag[:] += abs * temp.imag[:]

	def adjoint(self, add, model, data):
		if (not add): model.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		abs = np.abs(datNd)
		# self.real[:] = datNd.imag[:]
		self.real[:] = np.sin(np.arctan2(datNd.imag[:],datNd.real[:]))
		temp = hilbert(self.real,N=self.npad,axis=0)[0:datNd.shape[0],:]
		modNd.real[:] += datNd.real[:] - abs*temp.imag[:]
		modNd.imag[:] += 0

class HilbertC(Op.Operator):
	"""docstring for Hilbert."""

	def __init__(self, model, data, npad=0):
		self.cppMode = pyOp.Hilbert(npad)
		self.setDomainRange(model,data)

	def forward(self, add, model, data):
		self.cppMode.forward(model.cppMode,data.cppMode,add)

	def adjoint(self, add, model, data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)


def h(x, l, tau):
	return l*x - np.sqrt(l*l*x*x + tau**2)

def dh(x, l, tau):
	return l - l*l*x/np.sqrt(l*l*x*x + tau**2)

class HyperbolicPenalty(Op.Operator):
	"""docstring for HyperbolicPenalty."""

	def __init__(self, model, data, l=1, tau=.001):
		self.setDomainRange(model,data)
		self.l = l / 2
		self.tau = tau

	def forward(self, add, model, data):
		if not add: data.scale(0)
		# slownesss
		modNd = model[0].getNdArray()
		datNd = data[0].getNdArray()
		datNd.real[:] += modNd.real[:]
		datNd.imag[:] += h(modNd.imag, self.l, self.tau)
		# density 
		modNd = model[1].getNdArray()
		datNd = data[1].getNdArray()
		datNd[:] += modNd[:]

class Softclip(Op.Operator):
	"""docstring for Softclip."""

	def __init__(self, model, data, l=1, tau=.001):
		self.setDomainRange(model,data)
		self.l = l / 2
		self.tau = tau
		self.bg = np.copy(model[0].getNdArray().imag)

	def forward(self, add, model, data):
		if not add: data.scale(0)
		# slowness
		modNd = model[0].getNdArray()
		datNd = data[0].getNdArray()
		datNd.real[:] += modNd.real[:]
		sc = dh(self.bg, self.l, self.tau)
		datNd.imag[:] += modNd.imag[:]*sc
		# density 
		modNd = model[1].getNdArray()
		datNd = data[1].getNdArray()
		datNd[:] += modNd[:]
		

	def adjoint(self, add, model, data):
		self.forward(add,data,model)

	def set_background(self,bg):
		self.bg[:] = bg[0].getNdArray().imag

class Softmax(Op.Operator):
	def __init__(self, model, data, xmax, l=1, tau=.001):
		self.setDomainRange(model,data)
		self.l = l / 2
		self.tau = tau
		self.xmax = xmax

	def forward(self, add, model, data):
		if not add: data.scale(0)
		# slownesss
		modNd = model[0].getNdArray()
		datNd = data[0].getNdArray()
		datNd.imag[:] += modNd.imag[:]
		datNd.real[:] += h(modNd.real - self.xmax, self.l, self.tau) + self.xmax
		# density 
		modNd = model[1].getNdArray()
		datNd = data[1].getNdArray()
		datNd[:] += modNd[:]

class Softmin(Op.Operator):
	def __init__(self, model, data, xmin, l=1, tau=.001):
		self.setDomainRange(model,data)
		self.l = l / 2
		self.tau = tau
		self.xmin = xmin

	def forward(self, add, model, data):
		if not add: data.scale(0)
		# slownesss
		modNd = model[0].getNdArray()
		datNd = data[0].getNdArray()
		datNd.imag[:] += modNd.imag[:]
		datNd.real[:] += -h(-(modNd.real - self.xmin), self.l, self.tau) + self.xmin
		# density 
		modNd = model[1].getNdArray()
		datNd = data[1].getNdArray()
		datNd[:] += modNd[:]

class dSoftmax(Op.Operator):
	"""docstring for Softclip."""

	def __init__(self, model, data, xmax, l=1, tau=.001):
		self.setDomainRange(model,data)
		self.l = l / 2
		self.tau = tau
		self.xmax = xmax
		self.bg = np.copy(model[0].getNdArray().real)

	def forward(self, add, model, data):
		if not add: data.scale(0)
		# slowness
		modNd = model[0].getNdArray()
		datNd = data[0].getNdArray()
		datNd.imag[:] += modNd.imag[:]
		sc = dh(self.bg - self.xmax, self.l, self.tau)
		datNd.real[:] += modNd.real[:]*sc
		# density 
		modNd = model[1].getNdArray()
		datNd = data[1].getNdArray()
		datNd[:] += modNd[:]

	def adjoint(self, add, model, data):
		self.forward(add,data,model)

	def set_background(self,bg):
		self.bg[:] = bg[0].getNdArray().real

class dSoftmin(Op.Operator):
	"""docstring for Softclip."""

	def __init__(self, model, data, xmin, l=1, tau=.001):
		self.setDomainRange(model,data)
		self.l = l / 2
		self.tau = tau
		self.xmin = xmin
		self.bg = np.copy(model[0].getNdArray().real)

	def forward(self, add, model, data):
		if not add: data.scale(0)
		# slowness
		modNd = model[0].getNdArray()
		datNd = data[0].getNdArray()
		datNd.imag[:] += modNd.imag[:]
		sc = dh(-(self.bg - self.xmin), self.l, self.tau)
		datNd.real[:] += modNd.real[:]*sc
		# density 
		modNd = model[1].getNdArray()
		datNd = data[1].getNdArray()
		datNd[:] += modNd[:]

	def adjoint(self, add, model, data):
		self.forward(add,data,model)

	def set_background(self,bg):
		self.bg[:] = bg[0].getNdArray().real


def SoftMinMax(model, data, xmin, xmax, l=1, tau=[0.001, 0.001]):
	maxOp = Op.NonLinearOperator(Softmax(model, data, xmax, l=l, tau=tau[0]),dSoftmax(model, data, xmax, l=l, tau=tau[0]))
	minOp = Op.NonLinearOperator(Softmin(model, data, xmin, l=l, tau=tau[1]),dSoftmin(model, data, xmin, l=l, tau=tau[1]))
	return Op.CombNonlinearOp(maxOp,minOp)


class ComplexExp(Op.Operator):
	"""docstring for Softclip."""

	def __init__(self, model, data):
		self.setDomainRange(model,data)

	def forward(self, add, model, data):
		if not add: data.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		datNd[:] += np.exp(modNd)

class LinComplexExp(Op.Operator):
	"""docstring for Softclip."""

	def __init__(self, model, data):
		self.bg = np.copy(model.getNdArray())
		self.setDomainRange(model,data)

	def forward(self, add, model, data):
		if not add: data.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		datNd[:] += modNd[:] * np.exp(self.bg)

	def adjoint(self, add, model, data):
		if not add: model.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		modNd[:] += datNd[:] * np.exp(np.conj(self.bg))

	def set_background(self,bg):
		self.bg[:] = bg.getNdArray()

class Decon(Op.Operator):
	"""docstring for Decon."""

	def __init__(self, model, data, wave, eps=1e-2):
		self.setDomainRange(model,data)
		waveNd = wave.getNdArray()
		abs = np.abs(np.fft.fft(waveNd,norm='ortho')[:model.shape[2]])
		self.w = 1 / (abs + eps * np.amax(abs))

	def forward(self, add, model, data):
		if not add: data.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		datNd[:] += np.apply_along_axis(lambda m: m * self.w, axis=-1, arr=modNd)

	def adjoint(self, add, model, data):
		if not add: model.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		modNd[:] += np.apply_along_axis(lambda m: m * self.w, axis=-1, arr=datNd)

class WindowImag(Op.Operator):

	def __init__(self, model, data):
		self.setDomainRange(model,data)

	def forward(self, add, model, data):
		if not add: data.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		datNd.imag[:] += modNd[:].imag

	def adjoint(self, add, model, data):
		if not add: model.scale(0)
		modNd = model.getNdArray()
		datNd = data.getNdArray()
		modNd.imag[:] += datNd[:].imag

