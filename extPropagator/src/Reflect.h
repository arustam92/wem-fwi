#pragma once
#include <Operator.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <FFT1.h>

namespace SEP {

class Reflect : public Operator<complex2DReg,complex2DReg>
{
public:
	Reflect(std::shared_ptr<float2DReg> slow) {
		r = slow->clone();
		r->scale(0.);
		nz = slow->getHyper()->getAxis(2).n;
		nx = slow->getHyper()->getAxis(1).n;

		// zero-incidence reflectivity
		for (long long iz=1; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {
		  	(*r->_mat)[iz][ix] = ((*slow->_mat)[iz-1][ix]-(*slow->_mat)[iz][ix])/
		  						 ((*slow->_mat)[iz-1][ix]+(*slow->_mat)[iz][ix]);
		  }
		}

	};


	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

		if (!add) data->scale(0.);
	// tbb::parallel_for(tbb::blocked_range<int>())
		for (long long iz=0; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {
		  	(*data->_mat)[iz][ix] = (*r->_mat)[iz][ix]*(*model->_mat)[iz][ix];
		  }
		}
	}

	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if(!add) model->scale(0.);
		for (long long iz=0; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {
		  	(*model->_mat)[iz][ix] = (*r->_mat)[iz][ix]*(*data->_mat)[iz][ix];
		  }
		}
	}

private:
	std::shared_ptr<float2DReg> r;
	int nz, nx;
};
}
