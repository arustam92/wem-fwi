#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <PSPI.h>

using namespace SEP;

void PSPI::forward(std::shared_ptr<complex1DReg> model,
					std::shared_ptr<complex1DReg> data, bool add) {

	  if(!add) data->scale(0.);

		pad->forward(model,model_pad,0);
	  fft_in->forward(model_pad,model_pad,1);

		for (int iref=0; iref < nref; ++iref) {

			ps->setSlow(ref->getRefSlow(_iz,iref));
			ps->forward(model_pad,_wfld_ref_pad,0);

			fft_in->adjoint(_wfld_ref_pad,_wfld_ref_pad,1);
			pad->adjoint(_wfld_ref, _wfld_ref_pad, 0);
			taper->forward(_wfld_ref,_wfld_ref,1);

			select->setLocation(ref->getRefLoc(_iz,iref));
			select->forward(_wfld_ref,data,1);
		}

}

void PSPI::adjoint(std::shared_ptr<complex1DReg> model,
					std::shared_ptr<complex1DReg> data, bool add) {

		if(!add) model->scale(0.);

		model_pad->scale(0.);

		for (int iref=0; iref < nref; ++iref) {

			select->setLocation(ref->getRefLoc(_iz,iref));
			select->adjoint(_wfld_ref,data,0);

			taper->adjoint(_wfld_ref,_wfld_ref,1);
			pad->forward(_wfld_ref,_wfld_ref_pad,0);
			fft_in->forward(_wfld_ref_pad,_wfld_ref_pad,1);

			ps->setSlow(ref->getRefSlow(_iz,iref));
			ps->adjoint(model_pad,_wfld_ref_pad,1);
		}

		fft_in->adjoint(model_pad,model_pad,1);
		pad->adjoint(model, model_pad, 1);

}
