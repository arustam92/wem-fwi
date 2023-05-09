import WEM
import genericIO
import pyOperator as Op
import numpy as np

class SuperWEM (Op.Operator):
	def __init__(self,model,data_list,par,wave,src_list,rec_list):
		"""	mainly for medical imaging,
			src, rec -- list of the src/rec files in .H format
		"""
		self.wemOp = []
		for i in range(len(src_list)):
			""" deal with coordinates first """
			src = genericIO.defaultIO.getVector(src_list[i])
			rec = genericIO.defaultIO.getVector(rec_list[i])
			coef = fit_line(src.getNdArray())
			angle = np.math.atan(coef[0])
			sx,sz = rotate(src.getNdArray()[0],src.getNdArray()[1],angle,coef[1])
			rx,rz = rotate(rec.getNdArray()[0],rec.getNdArray()[1],angle,coef[1])
			shift = Operator.Shift(sz,rz,0)

			""" deal with the models now """



			shift.forward(False,data_list[i],data_list[i])



			del src,rec
		# 	fit src, rec to line -> a,b
		# 	create local model + rotate
		# 	rotate src and rec
		# 	shift using src/rec
		#
		# 	self.wemOp.append(WEM.WEM(model[i],data[i],par,wave,src[i],rec[i]))


	# def forward(self,add,model,data):
	# 	update local_models
	# 	for i in range():
	# 		self.wemOp[i].forward(add,local_model[i],data[i])

	def adjoint(self,add,model,data):
		self.cppMode.adjoint(model.cppMode,data.cppMode,add)

# class Born (Op.Operator):
	# def __init__(self,model,data,wave,par,src=None,rec=None):
	# 	if src!=None and rec!=None:
	# 		self.cppMode = pyWEM.Born(wave.cppMode,model.cppMode,par,src.cppMode,rec.cppMode)
	# 	else:
	# 		self.cppMode = pyWEM.Born(wave.cppMode,model.cppMode,par)
	# 	self.setDomainRange(model,data)
	#
	# def forward(self,add,model,data):
	# 	self.cppMode.forward(model.cppMode,data.cppMode,add)
	#
	# def adjoint(self,add,model,data):
	# 	self.cppMode.adjoint(model.cppMode,data.cppMode,add)
	#
	# def setBgSlow(self,slow):
	# 	self.cppMode.setBgSlow(slow.cppMode)


def fit_line(xz_arr):
	return np.polyfit(xz_arr[0],xz_arr[1],1)

def rotate(x,z,angRad,b=0):
	xr = np.cos(angRad)*x + np.sin(angRad)*(z - b)
	zr = -np.sin(angRad)*x + np.cos(angRad)*(z - b)
	return xr,zr
