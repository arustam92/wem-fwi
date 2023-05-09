
#include <OneWay.h>
#include <cstring>
#include <algorithm>
#include <PSPI.h>
#include <SSF.h>

using namespace SEP;

OneWay::OneWay(const std::shared_ptr<float2DReg> slow, const std::shared_ptr<paramObj> par,
			std::shared_ptr<RefSampler> ref) {

	nz = slow->getHyper()->getAxis(2).n;
	nx = slow->getHyper()->getAxis(1).n;
	oz = 0;
	orec = 0;

	// if (par->getFloat("osz",0) < slow->getHyper()->getAxis(2).o) oz = 0;
	// else oz = (par->getFloat("osz",0) - slow->getHyper()->getAxis(2).o)/slow->getHyper()->getAxis(2).d;
	// if (par->getFloat("orz",0) < slow->getHyper()->getAxis(2).o) orec = 0;
	// else orec = (par->getFloat("orz",0) - slow->getHyper()->getAxis(2).o)/slow->getHyper()->getAxis(2).d;
	_wfld_prev = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
	_wfld_temp = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
	_wfld_next = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));

	if (par->getString("prop")=="pspi") {
		prop = std::make_shared<PSPI>(slow,par,ref);
	}
	else {
		prop = std::make_shared<SSF>(slow,par,ref);
	}

}


void Down::forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if (!add) data->scale(0.);
	// boundaries
		std::copy_n(model->_mat->data(),nx,data->_mat->data());

	for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_prev	->_mat->data());
		std::copy_n(model->_mat->data()+(iz+1)*nx,nx,_wfld_next->_mat->data());

		prop->setDepth(iz);
		prop->forward(_wfld_prev,_wfld_temp,0);

		_wfld_next->scaleAdd(_wfld_temp,1.,1.);

		std::copy_n(_wfld_next->_mat->data(),nx,data->_mat->data()+(iz+1)*nx);

	}

}

void Down::adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	// boundaries
	if(!add) model->scale(0.);
	std::copy_n(data->_mat->data()+(nz-1)*nx,nx,model->_mat->data()+(nz-1)*nx);

		for (int iz=bounds[1]; iz>bounds[0]; iz--) {

			std::copy_n(data->_mat->data()+iz*nx,nx,_wfld_temp->_mat->data());
			std::copy_n(data->_mat->data()+(iz-1)*nx,nx,_wfld_next->_mat->data());

			prop->setDepth(iz-1);
			prop->adjoint(_wfld_prev,_wfld_temp,0);

			_wfld_next->scaleAdd(_wfld_prev,1.,1.);

			std::copy_n(_wfld_next->_mat->data(),nx,model->_mat->data()+(iz-1)*nx);
	}
}

void Up::forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if (!add) data->scale(0.);
	// boundaries
	std::copy_n(model->_mat->data()+(nz-1)*nx,nx,data->_mat->data()+(nz-1)*nx);

	for (int iz=bounds[0]; iz>bounds[1]; iz--) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_prev->_mat->data());
		std::copy_n(model->_mat->data()+(iz-1)*nx,nx,_wfld_next->_mat->data());

		prop->setDepth(iz);
		prop->forward(_wfld_prev,_wfld_temp,0);

		_wfld_next->scaleAdd(_wfld_temp,1.,1.);
		std::copy_n(_wfld_next->_mat->data(),nx,data->_mat->data()+(iz-1)*nx);
	}
}

void Up::adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if (!add) model->scale(0.);
	// boundaries
	std::copy_n(data->_mat->data(),nx,model->_mat->data());

	for (int iz=bounds[1]; iz<bounds[0]; iz++) {

		std::copy_n(data->_mat->data()+iz*nx,nx,_wfld_temp->_mat->data());
		std::copy_n(data->_mat->data()+(iz+1)*nx,nx,_wfld_next->_mat->data());

		prop->setDepth(iz+1);
		prop->adjoint(_wfld_prev,_wfld_temp,0);

		_wfld_next->scaleAdd(_wfld_prev,1.,1.);
		std::copy_n(_wfld_next->_mat->data(),nx,model->_mat->data()+(iz+1)*nx);
	}

}
