#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <StatPhase.h>

using namespace SEP;

void StatPhase::forward(std::shared_ptr<complex1DReg> model,
					std::shared_ptr<complex1DReg> data, bool add) {

	  if(!add) data->scale(0.);

		constexpr std::complex<float> I = {0, 1};
		constexpr std::complex<float> pi4 = {0, pi/4};
		float k1, k2;

		std::complex<float> phi1, phi2;
		float sc;

	  fft_out->forward(model,model_kx,0);

		for (int i=0; i < nx; ++i) {

			std::complex<float> ws = w*_slow_slice[i];

			k1 = ws.real() * scale1[i];
			k2 = -k1;

			if (k1 < 0) k1 += (nx-1)*dk;
			else k2 += (nx-1)*dk;

			int ik1 = (k1 - k[0]) / dk;
			int ik2 = (k2 - k[0]) / dk;

			phi1 = I * ws * scale21[i] - pi4;
			phi2 = I * ws * scale22[i] - pi4;
			sc = std::sqrt(std::abs(ws * scale3[i]));
			(*data->_mat)[i] += (*model_kx->_mat)[ik1] * std::exp(phi1) * sc;
			(*data->_mat)[i] += (*model_kx->_mat)[ik2] * std::exp(phi2) * sc;

		}
}

void StatPhase::adjoint(std::shared_ptr<complex1DReg> model,
					std::shared_ptr<complex1DReg> data, bool add) {

		// if(!add) {model->scale(0.);model_kx->scale(0.);}
		//
		// for (int iref=0; iref < nref; ++iref) {
		//
		// 	select->setLocation(ref->getRefLoc(_iz,iref));
		// 	select->adjoint(_wfld_ref,data,0);
		//
		// 	ss->setLocation(ref->getRefLoc(_iz,iref));
		// 	ss->setSlow(ref->getRefSlow(_iz,iref));
		// 	ss->adjoint(_wfld_ref,_wfld_ref,1);
		//
		// 	taper->adjoint(_wfld_ref,_wfld_ref,1);
		// 	fft_in->forward(_wfld_ref,_wfld_ref,1);
		//
		// 	ps->setSlow(ref->getRefSlow(_iz,iref));
		// 	ps->adjoint(model_kx,_wfld_ref,1);
		// }
		//
		// fft_out->adjoint(model,model_kx,1);

}
