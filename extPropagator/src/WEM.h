#pragma once
#include <Operator.h>
#include <float1DReg.h>
#include <float2DReg.h>
#include <float3DReg.h>
#include <Reflect.h>
#include <Phshift.h>
#include <Injection.h>
#include <FFT1.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <complex3DReg.h>
#include <Acquisition.h>
#include <omp.h>

#include <OneWay.h>

namespace SEP {

class WEM : public Operator<float1DReg,float3DReg>, Acquisition
{
public:
	WEM(std::shared_ptr<float3DReg> slow, std::shared_ptr<paramObj> par,
			std::shared_ptr<float2DReg> scoord = nullptr, std::shared_ptr<float2DReg> rcoord = nullptr);

	void forward(std::shared_ptr<float1DReg> model, std::shared_ptr<float3DReg> data, bool add);
	void adjoint(std::shared_ptr<float1DReg> model, std::shared_ptr<float3DReg> data, bool add);

	// nonlinear w.r.t. slowness
	WEM(std::shared_ptr<float1DReg> wave, std::shared_ptr<float3DReg> modspace, std::shared_ptr<paramObj> par,
			std::shared_ptr<float2DReg> scoord = nullptr, std::shared_ptr<float2DReg> rcoord = nullptr);
	void forward(std::shared_ptr<float3DReg> model, std::shared_ptr<float3DReg> data, bool add);
	// for complex-valued data (multi-scale approach)
	WEM(std::shared_ptr<complex1DReg> wave, std::shared_ptr<float2DReg> modspace, std::shared_ptr<paramObj> par,
			std::shared_ptr<float2DReg> scoord = nullptr, std::shared_ptr<float2DReg> rcoord = nullptr);
	void forward(std::shared_ptr<float3DReg> model, std::shared_ptr<complex3DReg> data, bool add);

	void setFminFmax(int _fmin, int _fmax) {
		fmin = _fmin;
		fmax = _fmax;
	}

	void setWfld(std::shared_ptr<complex3DReg> wavefield) {_full_wfld_w = wavefield;}
	std::shared_ptr<float3DReg> getWfld() {
		_full_wfld_t = std::make_shared<float3DReg> (_full_wfld_w->getHyper());
		FFT1 fft (_full_wfld_t,_full_wfld_w,1,1);
		fft.adjoint(_full_wfld_t,_full_wfld_w,0);
		return _full_wfld_t;
	}
	void snap_wfld(std::shared_ptr<complex2DReg> wfld,int i) {
		std::vector<axis> axes = _full_wfld_w->getHyper()->getAxes();
		for (size_t iz = 0; iz < axes[2].n; iz++) {
			for (size_t ix = 0; ix < axes[1].n; ix++) {
				(*_full_wfld_w->_mat)[iz][ix][i] += (*wfld->_mat)[iz][ix];
			}
		}
	}

	std::shared_ptr<complex3DReg> _full_wfld_w;
	std::shared_ptr<float3DReg> _full_wfld_t;

private:

	std::shared_ptr<Reflect> reflect;
	std::shared_ptr<float3DReg> _slow; // for linear op
	std::shared_ptr<complex1DReg> _wave_f; //for nonlinear op
	std::shared_ptr<RefSampler> ref;

	std::shared_ptr<paramObj> param;

	boost::multi_array<float,1> freq;

	bool wantwfld;
	int fmax, fmin;
	float fMIN, fMAX;

	int nref;
	int nshots;

};

}
