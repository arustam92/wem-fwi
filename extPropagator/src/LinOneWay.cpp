#include <LinOneWay.h>
#include <ic.h>
// #include "illumination.h"

using namespace SEP;

// model are the scatterers
// here the scattering is done using Taylor expansion (see 3dseis)

LinOneWay::LinOneWay(std::shared_ptr<float2DReg> slow,std::shared_ptr<paramObj> par,
	std::shared_ptr<complex2DReg> bg_wfld, std::shared_ptr<OneWay> oneway) {

	_slow = slow;
	_bg_wfld = bg_wfld;

	nz = slow->getHyper()->getAxis(2).n;
	nx = slow->getHyper()->getAxis(1).n;

	ntaylor = par->getInt("ntaylor",1);
	sc = std::make_shared<Scatter> (slow,ntaylor);
	_dz = slow->getHyper()->getAxis(2).d;

	_wfld_sc = std::make_shared<complex2DReg>(slow->getHyper());
	_wfld_prev = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
	_wfld_next = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));

	prop = oneway->getProp();
	ic = std::make_shared<IC>(bg_wfld);
	// if (par->getBool("illum",false)) illum = std::make_shared<Illumination>(bg_wfld->getHyper());

	bounds = oneway->getBounds();

}

void LinDown::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if(!add) data->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->_mat->data(),nx,data->_mat->data());

	for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		for (int ix=0; ix<nx; ix++) {
			(*_wfld_prev->_mat)[ix] = (*_bg_wfld->_mat)[iz][ix] * (*model->_mat)[iz][ix];
		}

		// loop over taylor expansion
		sc->setDepth(iz);
		sc->forward(_wfld_prev,_wfld_next,false);

		prop->setDepth(iz);
		prop->forward(_wfld_next,_wfld_prev,false);

		std::copy_n(_wfld_prev->_mat->data(),nx,data->_mat->data()+(iz+1)*nx);


	}
}

void LinDown::adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if(!add) model->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->_mat->data(),nx,_wfld_sc->_mat->data()+bounds[1]*nx);

	for (int iz=bounds[1]; iz>bounds[0]; iz--) {

		std::copy_n(data->_mat->data()+iz*nx,nx,_wfld_prev->_mat->data());

		prop->setDepth(iz-1);
		prop->adjoint(_wfld_next,_wfld_prev,false);

		// loop over taylor expansion
		sc->setDepth(iz-1);
		sc->adjoint(_wfld_prev,_wfld_next,false);

		std::copy_n(_wfld_prev->_mat->data(),nx,_wfld_sc->_mat->data()+(iz-1)*nx);

	}

	ic->adjoint(model,_wfld_sc,true);
}

void LinUp::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if(!add) data->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->_mat->data(),nx,data->_mat->data()+bounds[0]*nx);

	for (int iz=bounds[0]; iz>bounds[1]; iz--) {

		for (int ix=0; ix<nx; ix++) {
			(*_wfld_prev->_mat)[ix] = (*_bg_wfld->_mat)[iz][ix] * (*model->_mat)[iz][ix];
		}

		// loop over taylor expansion
		sc->setDepth(iz);
		sc->forward(_wfld_prev,_wfld_next,false);

		prop->setDepth(iz);
		prop->forward(_wfld_next,_wfld_prev,false);

		std::copy_n(_wfld_prev->_mat->data(),nx,data->_mat->data()+(iz-1)*nx);


	}
}

void LinUp::adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if(!add) model->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->_mat->data(),nx,_wfld_sc->_mat->data());

	for (int iz=bounds[1]; iz<bounds[0]; iz++) {

		std::copy_n(data->_mat->data()+iz*nx,nx,_wfld_prev->_mat->data());

		prop->setDepth(iz+1);
		prop->adjoint(_wfld_next,_wfld_prev,false);

		// loop over taylor expansion
		sc->setDepth(iz+1);
		sc->adjoint(_wfld_prev,_wfld_next,false);

		std::copy_n(_wfld_prev->_mat->data(),nx,_wfld_sc->_mat->data()+(iz+1)*nx);

	}

	ic->adjoint(model,_wfld_sc,true);
}
