#pragma once
#include "complex2DReg.h"
#include <Operator.h>
#include <algorithm>

namespace SEP {

class dReflect : public Operator<complex2DReg,complex2DReg>
{
public:
	dReflect(std::vector<std::shared_ptr<complex2DReg>> model) : _model(model) {
		
		nz = _model[0]->getHyper()->getAxis(2).n;
		nx = _model[0]->getHyper()->getAxis(1).n;
		s1s2sq = model[0]->clone();
		r1r2sq = model[0]->clone();
		
		for (int iz=1; iz<nz; ++iz) 
			for (int ix=0; ix<nx; ++ix) {
			std::complex<float> v = std::pow(
																	(*model[1]->_mat)[iz][ix] * std::sqrt((*model[0]->_mat)[iz-1][ix])+
																	(*model[1]->_mat)[iz-1][ix] * std::sqrt((*model[0]->_mat)[iz][ix]),-2);

				(*s1s2sq->_mat)[iz][ix] = ((*model[1]->_mat)[iz-1][ix] * (*model[1]->_mat)[iz][ix]) / v;
				(*r1r2sq->_mat)[iz][ix] = 2.f * ((*model[0]->_mat)[iz-1][ix] * (*model[0]->_mat)[iz][ix]) / v;
			}
	};

	void forward(std::vector<std::shared_ptr<complex2DReg>> model, std::shared_ptr<complex2DReg> data, bool add) {

		if(!add) data->scale(0.);

		typedef std::complex<float> cmplx;
		cmplx a, b, c;
		

		for (int iz=1; iz < nz-1; iz++) {
		  for (int ix=0; ix < nx; ix++) {
				// slowness
				c = std::sqrt((*_model[0]->_mat)[iz][ix]/(*_model[0]->_mat)[iz-1][ix]);
				a = c * (*s1s2sq->_mat)[iz][ix];
				b = 1.f/c * (*s1s2sq->_mat)[iz][ix];
				(*data->_mat)[iz][ix] += a * (*model[0]->_mat)[iz-1][ix] - b * (*model[0]->_mat)[iz][ix];
				// density
				a = (*_model[1]->_mat)[iz-1][ix] * (*r1r2sq->_mat)[iz][ix];
				b = (*_model[1]->_mat)[iz][ix] * (*r1r2sq->_mat)[iz][ix];
				(*data->_mat)[iz][ix] += a * (*model[1]->_mat)[iz][ix] - b * (*model[1]->_mat)[iz-1][ix];
			}
		}
		for (int ix=0; ix < nx; ix++) {
			// a = c * (*s1s2sq->_mat)[iz][ix];
			// b = 1.f / (*_slow->_mat)[0][ix];
			// (*data->_mat)[0][ix] += a * (*model->_mat)[1][ix] - b * (*model->_mat)[0][ix];
			// (*data->_mat)[0][ix] = - b * (*model->_mat)[0][ix];
		}
	};

	void adjoint(std::vector<std::shared_ptr<complex2DReg>> model, std::shared_ptr<complex2DReg> data, bool add) {

		if(!add) std::for_each(model.begin(), model.end(), [](std::shared_ptr<complex2DReg>& v) { v->scale(0); });

		typedef std::complex<float> cmplx;
		cmplx a, b, c;

		for (int iz=1; iz < nz-1; iz++) {
		  for (int ix=0; ix < nx; ix++) {
				// slowness
				c = std::sqrt((*_model[0]->_mat)[iz][ix]/(*_model[0]->_mat)[iz-1][ix]);
				a = std::conj(c * (*s1s2sq->_mat)[iz][ix]);
				b = std::conj(1.f/c * (*s1s2sq->_mat)[iz][ix]);
				(*model[0]->_mat)[iz-1][ix] += a *	(*data->_mat)[iz][ix];
				(*model[0]->_mat)[iz][ix] += -b * (*data->_mat)[iz][ix];
				// density
				a = std::conj((*_model[1]->_mat)[iz-1][ix] * (*r1r2sq->_mat)[iz][ix]);
				b = std::conj((*_model[1]->_mat)[iz][ix] * (*r1r2sq->_mat)[iz][ix]);
				(*model[1]->_mat)[iz][ix] += a * (*data->_mat)[iz][ix];
				(*model[1]->_mat)[iz-1][ix] += -b * (*data->_mat)[iz][ix];
		  }
		}
		for (int ix=0; ix < nx; ix++) {
			// b = 1.f / (*_slow->_mat)[0][ix];
			// (*model->_mat)[0][ix] += -std::conj(b) * (*data->_mat)[0][ix];
		}
	};

protected:
	const std::vector<std::shared_ptr<complex2DReg>> _model;
	std::shared_ptr<complex2DReg> s1s2sq, r1r2sq;
	int nz,nx;
};

class dTransmission : public dReflect
{
public:
	dTransmission(const std::vector<std::shared_ptr<complex2DReg>> slow) : dReflect(slow) {
		s1s2sq->scale(-1.f);
	};
};

}
