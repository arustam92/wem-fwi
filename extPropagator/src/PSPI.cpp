#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <PSPI.h>

using namespace SEP;

void PSPI::forward(std::shared_ptr<complex1DReg> model,
					std::shared_ptr<complex1DReg> data, bool add) {

	  if(!add) data->scale(0.);

	  fft_out->forward(model,model_kx,0);

		for (int iref=0; iref < nref; ++iref) {

			ps->setSlow(ref->getRefSlow(_iz,iref));
			ps->forward(model_kx,_wfld_ref,0);

			fft_in->adjoint(_wfld_ref,_wfld_ref,1);
			taper->forward(_wfld_ref,_wfld_ref,1);

			select->setLocation(ref->getRefLoc(_iz,iref));
			select->forward(_wfld_ref,data,1);
		}

}

void PSPI::adjoint(std::shared_ptr<complex1DReg> model,
					std::shared_ptr<complex1DReg> data, bool add) {

		if(!add) {model->scale(0.);model_kx->scale(0.);}

		for (int iref=0; iref < nref; ++iref) {

			select->setLocation(ref->getRefLoc(_iz,iref));
			select->adjoint(_wfld_ref,data,0);

			taper->adjoint(_wfld_ref,_wfld_ref,1);
			fft_in->forward(_wfld_ref,_wfld_ref,1);

			ps->setSlow(ref->getRefSlow(_iz,iref));
			ps->adjoint(model_kx,_wfld_ref,1);
		}

		fft_out->adjoint(model,model_kx,1);

}
