#include <LinTrans.h>
#include <ic.h>
#include <algorithm>
#include <functional>

using namespace SEP;

void LinTransDown::forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add) {

	if(!add) data->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->getVals(),nx,data->getVals());

	auto tmp = tmp_wfld->clone();
	ic->forward(model,tmp,false);
	dtrans->forward(tmp,tmp_wfld,false);

	for (int iz=bounds[0]; iz<bounds[1]; iz++) {

		std::copy_n(tmp_wfld->getVals()+iz*nx,nx,_wfld_next->getVals());

		// ic->setDepth(iz);
		// ic->forward(_wfld_next,_wfld_prev,false);

		prop->setDepth(iz);
		prop->forward(_wfld_next,_wfld_prev,false);

		std::transform(_wfld_prev->getVals(),_wfld_prev->getVals()+nx,
										data->getVals()+(iz+1)*nx, data->getVals()+(iz+1)*nx,
										[](const std::complex<float> &i, const std::complex<float> &j) {return i + j; } );
	}

}

void LinTransDown::adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add) {

	if(!add) model->scale(0.);

	_wfld_prev->scale(0.);
	std::copy_n(_wfld_prev->getVals(),nx,model->getVals()+bounds[1]*nx);

	for (int iz=bounds[1]; iz>bounds[0]; iz--) {

		std::copy_n(data->getVals()+iz*nx,nx,_wfld_prev->getVals());

		prop->setDepth(iz-1);
		prop->adjoint(_wfld_next,_wfld_prev,false);

		// ic->setDepth(iz-1);
		// ic->adjoint(_wfld_next,_wfld_prev,false);

		std::transform(_wfld_next->getVals(),_wfld_next->getVals()+nx, tmp_wfld->getVals()+(iz-1)*nx,
										[](const std::complex<float> &partial) {return partial;	} );

	}
	auto tmp = tmp_wfld->clone();
	dtrans->adjoint(tmp,tmp_wfld,false);
	ic->adjoint(model,tmp,add);


}

void LinTransUp::forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add) {

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

		trans->setDepth(iz);
		trans->forward(_wfld_prev,_wfld_prev,1);

		std::transform(_wfld_prev->getVals(),_wfld_prev->getVals()+nx,
										data->getVals()+(iz-1)*nx, data->getVals()+(iz-1)*nx,
										[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );

	}
}

void LinTransUp::adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add) {

	if(!add) model->scale(0.);

	_wfld_prev->scale(0.);

	for (int iz=bounds[1]; iz<bounds[0]; iz++) {

		std::copy_n(data->getVals()+iz*nx,nx,_wfld_prev->getVals());

		trans->setDepth(iz+1);
		trans->adjoint(_wfld_prev,_wfld_prev,1);

		prop->setDepth(iz+1);
		prop->adjoint(_wfld_next,_wfld_prev,false);

		// loop over taylor expansion
		sc->setDepth(iz+1);
		sc->adjoint(_wfld_prev,_wfld_next,false);

		ic->setDepth(iz+1);
		ic->adjoint(_wfld_next,_wfld_prev,false);

		std::transform(_wfld_next->getVals(),_wfld_next->getVals()+nx,
										model->getVals()+(iz+1)*nx, model->getVals()+(iz+1)*nx,
										[](const std::complex<float> &partial, const std::complex<float> &full) {return partial+full;	} );

	}
}
