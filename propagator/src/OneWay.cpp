
#include <OneWay.h>
#include <cstring>
#include <algorithm>
#include <PSPI.h>
#include <SSF.h>
#include "StatPhase.h"

using namespace SEP;

OneWay::OneWay(const std::shared_ptr<complex2DReg>& slow, const std::shared_ptr<paramObj>& par, const std::shared_ptr<RefSampler>& ref) {

	nz = slow->getHyper()->getAxis(2).n;
	nx = slow->getHyper()->getAxis(1).n;
	oz = 0;
	orec = 0;

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


void Down::forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add) {

	if (!add) data->scale(0.);
	// boundaries
		std::copy_n(model->_mat->data(),nx,data->_mat->data());

	for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_prev	->_mat->data());

		prop->setDepth(iz);
		prop->forward(_wfld_prev,_wfld_temp,0);

		std::transform(_wfld_temp->_mat->data(),_wfld_temp->_mat->data()+nx,
										model->_mat->data()+(iz+1)*nx, data->_mat->data()+(iz+1)*nx,
										[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );

	}

}

void Down::adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add) {

	// boundaries
	if(!add) model->scale(0.);
	std::copy_n(data->_mat->data()+(nz-1)*nx,nx,model->_mat->data()+(nz-1)*nx);

		for (int iz=bounds[1]; iz>bounds[0]; iz--) {

			std::copy_n(data->_mat->data()+iz*nx,nx,_wfld_temp->_mat->data());

			prop->setDepth(iz-1);
			prop->adjoint(_wfld_prev,_wfld_temp,0);

			std::transform(_wfld_prev->_mat->data(),_wfld_prev->_mat->data()+nx,
											data->_mat->data()+(iz-1)*nx, model->_mat->data()+(iz-1)*nx,
											[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );
	}
}

void Up::forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add) {

	if (!add) data->scale(0.);
	// boundaries
	std::copy_n(model->_mat->data()+(nz-1)*nx,nx,data->_mat->data()+(nz-1)*nx);

	for (int iz=bounds[0]; iz>bounds[1]; iz--) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_prev->_mat->data());

		prop->setDepth(iz);
		prop->forward(_wfld_prev,_wfld_temp,0);

		std::transform(_wfld_temp->_mat->data(),_wfld_temp->_mat->data()+nx,
										model->_mat->data()+(iz-1)*nx, data->_mat->data()+(iz-1)*nx,
										[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );
	}
}

void Up::adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add) {

	if (!add) model->scale(0.);
	// boundaries
	std::copy_n(data->_mat->data(),nx,model->_mat->data());

	for (int iz=bounds[1]; iz<bounds[0]; iz++) {

		std::copy_n(data->_mat->data()+iz*nx,nx,_wfld_temp->_mat->data());

		prop->setDepth(iz+1);
		prop->adjoint(_wfld_prev,_wfld_temp,0);

		std::transform(_wfld_prev->_mat->data(),_wfld_prev->_mat->data()+nx,
										data->_mat->data()+(iz+1)*nx, model->_mat->data()+(iz+1)*nx,
										[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );
	}

}

void DownSource::forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add) {

	if (!add) data->scale(0.);
	// boundaries
		std::copy_n(model->_mat->data(),nx,data->_mat->data());

		for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_prev	->_mat->data());

		src_filter->setDepth(iz);
		src_filter->forward(_wfld_prev, _wfld_temp, 0);

		std::transform(_wfld_temp->_mat->data(),_wfld_temp->_mat->data()+nx,
										model->_mat->data()+iz*nx,
										[](const std::complex<float> &i) {return i; } );

	}

	Down::forward(model, data, add);

}

void DownSource::adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add) {
	
	Down::adjoint(model, data, add);

	// boundaries
		std::copy_n(model->_mat->data(),nx,data->_mat->data());

		for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_prev	->_mat->data());

		src_filter->setDepth(iz);
		src_filter->adjoint(_wfld_temp, _wfld_prev, 0);

		std::transform(_wfld_temp->_mat->data(),_wfld_temp->_mat->data()+nx,
										model->_mat->data()+iz*nx,
										[](const std::complex<float> &i) {return i; } );

	}

}

void UpSource::forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add) {

	if (!add) data->scale(0.);
	// boundaries
		std::copy_n(model->_mat->data(),nx,data->_mat->data());

		for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_prev	->_mat->data());

		src_filter->setDepth(iz);
		src_filter->forward(_wfld_prev, _wfld_temp, 0);

		std::transform(_wfld_temp->_mat->data(),_wfld_temp->_mat->data()+nx,
										model->_mat->data()+iz*nx,
										[](const std::complex<float> &i) {return i; } );

	}

	Up::forward(model, data, add);

}

void UpSource::adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add) {

	Up::adjoint(model, data, add);

	// boundaries
		std::copy_n(model->_mat->data(),nx,data->_mat->data());

		for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_prev	->_mat->data());

		src_filter->setDepth(iz);
		src_filter->adjoint(_wfld_temp, _wfld_prev, 0);

		std::transform(_wfld_temp->_mat->data(),_wfld_temp->_mat->data()+nx,
										model->_mat->data()+iz*nx,
										[](const std::complex<float> &i) {return i; } );

	}


}
