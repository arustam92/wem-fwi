#include "BornTomo.h"
#include "illumination.h"
#include "Weight.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

using namespace SEP;

void BornTomo::fwd_propagate(const std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data) {

	int nf = _refl_->getHyper()->getAxis(3).n;

	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
		[=] (const tbb::blocked_range<int> &r) {

		for (int i=r.begin(); i<r.end(); i++) {

	for (int is=0; is<nshots; ++is) {

		std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ng",0)));
		std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX()));
		std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
		std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
		std::shared_ptr<Down> down (new Down(_slow,param,ref));
		std::shared_ptr<Up> up (new Up (_slow,param,ref));

		std::unique_ptr<Born_up> born_up (new Born_up(_slow,param,inj_rec,up,_bg_wfld,_sc_wfld));

		std::shared_ptr<complex2DReg> reflOneFreq = std::make_shared<complex2DReg>(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2));
		sliceModel(_refl_,reflOneFreq,i);
		std::shared_ptr<Operator<complex2DReg,complex2DReg>> reflect = std::make_shared<Weight<complex2DReg>>(reflOneFreq);

		std::unique_ptr<Born_down> born_down(new Born_down (_slow,param,inj_rec,down,up,reflect,_bg_wfld,_sc_wfld));

			inj_rec->setShot(is);
			inj_src->setShot(is);

				// setting up
				inj_src->setStep(index[i]);
				inj_rec->setStep(index[i]);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);
				born_down->setFreq(freq[i]);
				born_up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(getWavelet(),_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// down-scattering
				born_down->forward(model,data,1);

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				// up-scattering
				born_up->forward(model,data,1);
			}
		// });
		}
	});
}



void BornTomo::adj_propagate(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data) {

	int nf = _refl_->getHyper()->getAxis(3).n;

	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
		[=] (const tbb::blocked_range<int> &r) {

		for (int i=r.begin(); i<r.end(); i++) {

			std::shared_ptr<complex2DReg> reflOneFreq = std::make_shared<complex2DReg>(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2));
			sliceModel(_refl_,reflOneFreq,i);
			std::shared_ptr<Operator<complex2DReg,complex2DReg>> reflect = std::make_shared<Weight<complex2DReg>>(reflOneFreq);

			for (int is=0; is<nshots; ++is) {
				std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ng",0)));
				std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX()));
				std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
				std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
				std::shared_ptr<Down> down (new Down(_slow,param,ref));
				std::shared_ptr<Up> up (new Up (_slow,param,ref));

				std::unique_ptr<Born_up> born_up (new Born_up(_slow,param,inj_rec,up,_bg_wfld,_sc_wfld));

				std::unique_ptr<Born_down> born_down(new Born_down (_slow,param,inj_rec,down,up,reflect,_bg_wfld,_sc_wfld));


			inj_rec->setShot(is);
			inj_src->setShot(is);

				// setting up
				inj_src->setStep(index[i]);
				inj_rec->setStep(index[i]);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);
				born_down->setFreq(freq[i]);
				born_up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(getWavelet(),_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// illum->accumulate(_bg_wfld);

				// down-scattering
				born_down->adjoint(model,data,1);
				//
				//
				// // up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);
				//
				// illum->accumulate(_bg_wfld);
				// up-scattering
				born_up->adjoint(model,data,1);
			}
		}
	});
}

void BornTomo::sliceModel(std::shared_ptr<complex3DReg> input,std::shared_ptr<complex2DReg> output, int ifreq) {
	// auto view4d = (c4Dview)(*input->_mat)[indices[ifreq][range()][range()][range()]];
	// (*output->_mat) = view4d;

	std::copy_n(input->_mat->data()+ifreq*nz*nx,nz*nx,output->_mat->data());
}
void BornTomo::sliceModel(std::shared_ptr<complex2DReg> input,std::shared_ptr<complex3DReg> output, int ifreq) {
	// auto view4d = (c4Dview)(*output->_mat)[indices[ifreq][range()][range()][range()]];
	// view4d = (*input->_mat);

	// std::copy_n(input->_mat->data(),nz*nx,output->_mat->data()+ifreq*nz*nx);
	std::transform(input->_mat->data(),input->_mat->data()+nx*nz,
									output->_mat->data()+ifreq*nz*nx, output->_mat->data()+ifreq*nz*nx,
									[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );
}
