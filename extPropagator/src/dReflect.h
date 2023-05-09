#pragma once
#include <Operator.h>

namespace SEP {

class dReflect : public Operator<float2DReg,float2DReg>
{
public:
	dReflect(std::shared_ptr<float2DReg> slow) : _slow(slow) {
		ind = _slow->getHyper()->getAxis(2).n-1;
		nz = _slow->getHyper()->getAxis(2).n;
		nx = _slow->getHyper()->getAxis(1).n;
		s1s2sq = std::make_shared<float2DReg>(slow->getHyper());
		for (int iz=1; iz<nz; ++iz) {
			for (int ix=0; ix<nx; ++ix) {
				(*s1s2sq->_mat)[iz][ix] = 1./std::pow((*_slow->_mat)[iz-1][ix]+(*_slow->_mat)[iz][ix],2);
			}
		}
	};

	void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<float2DReg> data, bool add) {

		if(!add) data->scale(0.);

		for (long long iz=1; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {

			(*data->_mat)[iz-1][ix] += 2.*(*_slow->_mat)[iz][ix] * (*s1s2sq->_mat)[iz][ix] *
										(*model->_mat)[iz-1][ix] -

		  						  	   2.*(*_slow->_mat)[iz-1][ix] * (*s1s2sq->_mat)[iz][ix] *
		  						  	   (*model->_mat)[iz][ix];
			}

		}
		for (long long ix=0; ix < nx; ix++) {

			(*data->_mat)[ind][ix] = 0;
		}
	};

	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<float2DReg> data, bool add) {

		if(!add) model->scale(0.);

		for (long long iz=1; iz < nz; iz++) {
		  for (long long ix=0; ix < nx; ix++) {

			(*model->_mat)[iz-1][ix] += 2.*(*_slow->_mat)[iz][ix] * (*s1s2sq->_mat)[iz][ix] *
										(*data->_mat)[iz-1][ix];
			(*model->_mat)[iz][ix] += -2.*(*_slow->_mat)[iz-1][ix]* (*s1s2sq->_mat)[iz][ix]  *
		  						  	   (*data->_mat)[iz-1][ix];
		  }
			for (long long ix=0; ix < nx; ix++) {

				(*model->_mat)[ind][ix] = 0;
			}
		}
	};

private:
	std::shared_ptr<float2DReg> _slow;
	std::shared_ptr<float2DReg> s1s2sq;
	int ind;
	int nz,nx;
};

}
