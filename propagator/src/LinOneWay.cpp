#include <LinOneWay.h>
#include <ic.h>
#include <algorithm>
#include <functional>

using namespace SEP;

// model are the scatterers
// here the scattering is done using Taylor expansion (see 3dseis)

LinOneWay::LinOneWay(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,
			const std::shared_ptr<complex2DReg>& bg_wfld, const std::shared_ptr<OneWay>& oneway) : _slow(slow),_bg_wfld(bg_wfld) {

	nz = slow->getHyper()->getAxis(2).n;
	nx = slow->getHyper()->getAxis(1).n;

	ntaylor = par->getInt("ntaylor",1);
	sc = std::make_shared<Scatter> (slow,ntaylor,par);
	_dz = slow->getHyper()->getAxis(2).d;

	// _wfld_sc = std::make_shared<complex2DReg>(slow->getHyper());
	_wfld_prev = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
	_wfld_next = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));

	prop = oneway->getProp();
	ic = std::make_shared<IC>(bg_wfld);
	// if (par->getBool("illum",false)) illum = std::make_shared<Illumination>(bg_wfld->getHyper());

	bounds = oneway->getBounds();

}

void LinDown::forward(const std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if(!add) data->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->_mat->data(),nx,data->_mat->data());

	for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		std::copy_n(model->_mat->data()+iz*nx,nx,_wfld_next->_mat->data());

		ic->setDepth(iz);
		ic->forward(_wfld_next,_wfld_prev,false);

		// loop over taylor expansion
		sc->setDepth(iz);
		sc->forward(_wfld_prev,_wfld_next,false);

		prop->setDepth(iz);
		prop->forward(_wfld_next,_wfld_prev,false);

		std::transform(_wfld_prev->_mat->data(),_wfld_prev->_mat->data()+nx,
										data->_mat->data()+(iz+1)*nx, data->_mat->data()+(iz+1)*nx,
										[]( std::complex<float> &i,  std::complex<float> &j) {return i+j; } );
	}
}

void LinDown::adjoint(std::shared_ptr<complex2DReg> model, const std::shared_ptr<complex2DReg> data, bool add) {

	if(!add) model->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->_mat->data(),nx,model->_mat->data()+bounds[1]*nx);

	for (int iz=bounds[1]; iz>bounds[0]; iz--) {

		std::copy_n(data->_mat->data()+iz*nx,nx,_wfld_prev->_mat->data());

		prop->setDepth(iz-1);
		prop->adjoint(_wfld_next,_wfld_prev,false);

		// loop over taylor expansion
		sc->setDepth(iz-1);
		sc->adjoint(_wfld_prev,_wfld_next,false);

		ic->setDepth(iz-1);
		ic->adjoint(_wfld_next,_wfld_prev,false);

		std::transform(_wfld_next->_mat->data(),_wfld_next->_mat->data()+nx,
										model->_mat->data()+(iz-1)*nx, model->_mat->data()+(iz-1)*nx,
										[]( std::complex<float> &partial,  std::complex<float> &full) {return partial+full;	} );

	}

}

void LinUp::forward(const std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if(!add) data->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->getVals(),nx,data->getVals()+bounds[0]*nx);

	for (int iz=bounds[0]; iz>bounds[1]; iz--) {

		std::copy_n(model->getVals()+iz*nx,nx,_wfld_next->getVals());

		ic->setDepth(iz);
		ic->forward(_wfld_next,_wfld_prev,false);

		// loop over taylor expansion
		sc->setDepth(iz);
		sc->forward(_wfld_prev,_wfld_next,false);

		prop->setDepth(iz);
		prop->forward(_wfld_next,_wfld_prev,false);

		std::transform(_wfld_prev->getVals(),_wfld_prev->getVals()+nx,
										data->getVals()+(iz-1)*nx, data->getVals()+(iz-1)*nx,
										[]( std::complex<float> &i,  std::complex<float> &j) {return i+j; } );

	}
}

void LinUp::adjoint(std::shared_ptr<complex2DReg> model, const std::shared_ptr<complex2DReg> data, bool add) {

	if(!add) model->scale(0.);

	_wfld_prev->scale(0.);

	for (int iz=bounds[1]; iz<bounds[0]; iz++) {

		std::copy_n(data->getVals()+iz*nx,nx,_wfld_prev->getVals());

		prop->setDepth(iz+1);
		prop->adjoint(_wfld_next,_wfld_prev,false);

		// loop over taylor expansion
		sc->setDepth(iz+1);
		sc->adjoint(_wfld_prev,_wfld_next,false);

		ic->setDepth(iz+1);
		ic->adjoint(_wfld_next,_wfld_prev,false);

		std::transform(_wfld_next->getVals(),_wfld_next->getVals()+nx,
										model->getVals()+(iz+1)*nx, model->getVals()+(iz+1)*nx,
										[]( std::complex<float> &partial,  std::complex<float> &full) {return partial+full;	} );

	}
}