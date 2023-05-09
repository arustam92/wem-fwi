#pragma once
#include "complex2DReg.h"
#include <Operator.h>

namespace SEP {

class dTransmission : public Operator<complex2DReg,complex2DReg>
{
public:
	dTransmission(const std::shared_ptr<complex2DReg> slow) : _slow(slow) {
		nz = _slow->getHyper()->getAxis(2).n;
		nx = _slow->getHyper()->getAxis(1).n;
		s1s2sq = std::make_shared<complex2DReg>(slow->getHyper()->getAxis(1),slow->getHyper()->getAxis(2));
		for (int iz=1; iz<nz; ++iz) {
			for (int ix=0; ix<nx; ++ix) {
				(*s1s2sq->_mat)[iz][ix] = std::pow(std::sqrt((*_slow->_mat)[iz-1][ix])+
																					 std::sqrt((*_slow->_mat)[iz][ix]),-2);
			}
		}
	};

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

		if(!add) data->scale(0.);

		std::complex<float> a, b;

		for (long long iz=1; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {
				a = 2.f*(*_slow->_mat)[iz][ix] * (*s1s2sq->_mat)[iz][ix];
				b = 2.f*(*_slow->_mat)[iz-1][ix] * (*s1s2sq->_mat)[iz][ix];
				(*data->_mat)[iz][ix] += -a * (*model->_mat)[iz-1][ix] + b * (*model->_mat)[iz][ix];
			}
		}
		for (long long ix=0; ix < nx; ix++) {
			// a = float(2)*(*_slow->_mat)[0][ix] * (*s1s2sq->_mat)[0][ix];
			// b = float(2)*(*_slow->_mat)[1][ix] * (*s1s2sq->_mat)[0][ix];
			// (*data->_mat)[0][ix] += a * (*model->_mat)[1][ix] - b * (*model->_mat)[0][ix];
			// (*data->_mat)[0][ix] = 0;
		}
	};

	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

		if(!add) model->scale(0.);

		std::complex<float> a, b;

		for (long long iz=1; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {
				a = std::conj(2.f*(*_slow->_mat)[iz][ix] * (*s1s2sq->_mat)[iz][ix]);
				b = std::conj(2.f*(*_slow->_mat)[iz-1][ix] * (*s1s2sq->_mat)[iz][ix]);
				(*model->_mat)[iz-1][ix] += -a *	(*data->_mat)[iz][ix];
				(*model->_mat)[iz][ix] += b * (*data->_mat)[iz][ix];
		  }
		}
		for (long long ix=0; ix < nx; ix++) {
			// a = std::conj(float(2)*(*_slow->_mat)[nz-2][ix] * (*s1s2sq->_mat)[nz-2][ix]);
			// b = std::conj(float(2)*(*_slow->_mat)[1][ix] * (*s1s2sq->_mat)[0][ix]);
			// (*model->_mat)[nz-1][ix] += a *	(*data->_mat)[nz-2][ix];
			// (*model->_mat)[0][ix] += -b * (*data->_mat)[0][ix];
			// (*model->_mat)[0][ix] = 0;
			// (*model->_mat)[nz-1][ix] = 0;
		}
	};

private:
	const std::shared_ptr<complex2DReg> _slow;
	std::shared_ptr<complex2DReg> s1s2sq;
	int nz,nx;
};

}
