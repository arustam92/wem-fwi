#pragma once
#include <Operator.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <FFT1.h>

namespace SEP {

class Reflect : public Operator<complex2DReg,complex2DReg>
{
public:
	Reflect(std::vector<std::shared_ptr<complex2DReg>> model) {

		_model = model; 
		nz = model[0]->getHyper()->getAxis(2).n;
		nx = model[0]->getHyper()->getAxis(1).n;
		r.resize(boost::extents[nz][nx]);

		for (int iz=1; iz < nz; iz++) {
		  for (int ix=0; ix < nx; ix++) {
		  	r[iz][ix]  = ((*model[1]->_mat)[iz][ix]*std::sqrt((*model[0]->_mat)[iz-1][ix])-
											(*model[1]->_mat)[iz-1][ix]*std::sqrt((*model[0]->_mat)[iz][ix]))/
		  						 ((*model[1]->_mat)[iz][ix]*std::sqrt((*model[0]->_mat)[iz-1][ix])+
									 (*model[1]->_mat)[iz-1][ix]*std::sqrt((*model[0]->_mat)[iz][ix]));
		  }
		}

	};

	std::vector<std::shared_ptr<complex2DReg>>& getBg() {return _model;}


	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

		if (!add) data->scale(0.);
		for (int iz=0; iz < nz; iz++) 
		  for (int ix=0; ix < nx; ix++) 
		  	(*data->_mat)[iz][ix] = r[iz][ix]*(*model->_mat)[iz][ix];
	}

	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if(!add) model->scale(0.);
		for (int iz=0; iz < nz; iz++) {
		  for (int ix=0; ix < nx; ix++) {
		  	(*model->_mat)[iz][ix] = std::conj(r[iz][ix])*(*data->_mat)[iz][ix];
		  }
		}
	}

protected:
	std::vector<std::shared_ptr<complex2DReg>> _model;
	complex2D r;
	int nz, nx;
};

class Transmission : public Reflect
{
public:
	Transmission(std::vector<std::shared_ptr<complex2DReg>> slow) : Reflect(slow) {

		// for (int iz=1; iz < nz; iz++) {
		//   for (int ix=0; ix < nx; ix++) {
		//   	r[iz][ix]  = 2.f * std::sqrt((*slow->_mat)[iz][ix]) /
		//   						 (std::sqrt((*slow->_mat)[iz-1][ix])+std::sqrt((*slow->_mat)[iz][ix]));
		//   }
		// }

	};
};

}
