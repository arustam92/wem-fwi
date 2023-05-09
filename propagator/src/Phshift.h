#pragma once

#include <Operator.h>
#include <complex1DReg.h>
#include <float1DReg.h>
#include <cmath>
#include <hypercube.h>

namespace SEP {

class Scatter;

constexpr double pi = 4*std::atan(1.);

class Phshift : public Operator<complex1DReg,complex1DReg>
{
public:
	Phshift (float dz, std::shared_ptr<hypercube> hyp, int nref, float eps=0.04);

	void forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);

	inline void setRef(int iref) {_iref = iref;}
	inline void setFreq(float freq) {
		w = 2*pi*freq;
		// float vmin = 1500;
		// float sin = .5;
		// float kmax = 2*freq*sin / vmin;
		// nk = k.size() * kmax / kN;
		// ik0 = (kN - kmax) / dk;
	}
	inline void setSlow(std::complex<float> slow1d) {
		_s = slow1d;
	}
	inline void swapKz() {
		for (int i=0; i<kz.size(); ++i) {
			std::copy_n(kz[i]->data(),k.size(),kz_prev[i]->data());
		}
	}

	inline void lookupTable();

	void reflect(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
		if(!add) data->scale(0);
		std::complex<float> R;
		for (int ik = 0; ik < k.size(); ++ik) {
			R = ((*kz_prev[_iref])[ik] - (*kz[_iref])[ik])/((*kz_prev[_iref])[ik] + (*kz[_iref])[ik]);
			(*data->_mat)[ik] = R * (*model->_mat)[ik];
		}
	}

protected:
	float _eps;
	float w;
	float _dz;
	float1D k;
	complex1D _sqrt;
	// std::shared_ptr<complex1D> kz_prev, kz;
	std::shared_ptr<float1D> _amp;
	float1D epsilon;
	int ntable;
	float dtable, dk, kN;
	int _kstep_, nk, ik0;

	float min,max,std_dev;
	std::vector<std::shared_ptr<complex1D>> kz_prev, kz;
	complex1D sinc;
	float xmin, xmax;

	std::complex<float> _s;
	int _iref;
};



}
