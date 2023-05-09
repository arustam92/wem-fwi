#pragma once
#include <Operator.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <FFT1.h>

namespace SEP {

class Transmission : public Operator<complex2DReg,complex2DReg>
{
public:
	Transmission(std::shared_ptr<complex2DReg> slow) {

		_slow = slow;
		nz = slow->getHyper()->getAxis(2).n;
		nx = slow->getHyper()->getAxis(1).n;
		t.resize(boost::extents[nz][nx]);

		// zero-incidence reflectivity
		std::complex<float> x;
		for (long long iz=1; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {
				t[iz][ix]  = 2.f * std::sqrt((*slow->_mat)[iz][ix])) /
		  						 (std::sqrt((*slow->_mat)[iz-1][ix])+std::sqrt((*slow->_mat)[iz][ix]));
		  }
		}

	};

	std::shared_ptr<complex2DReg>& getSlow() {return _slow;}


	void forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add) {

		if (!add) data->scale(0.);
	// tbb::parallel_for(tbb::blocked_range<int>())
		for (long long iz=0; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {
		  	(*data->_mat)[iz][ix] = t[iz][ix]*(*model->_mat)[iz][ix];
		  }
		}
	}

	void adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add) {
		if(!add) model->scale(0.);
		for (long long iz=0; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {
		  	(*model->_mat)[iz][ix] = std::conj(t[iz][ix])*(*data->_mat)[iz][ix];
		  }
		}
	}

private:
	std::shared_ptr<complex2DReg> _slow;
	complex2D t;
	int nz, nx;
};
}
