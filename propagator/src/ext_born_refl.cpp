#include "BornRefl.h"
#include "illumination.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

using namespace SEP;

void BornRefl::fwd_propagate(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {

	int nf = model->getHyper()->getAxis(3).n;

	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
		[=] (const tbb::blocked_range<int> &r) {

		for (int i=r.begin(); i<r.end(); i++) {

			std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2));
			sliceModel(model,modOneFreq,i);

	// tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
	// 	[=](const tbb::blocked_range<int> &rs) {

	for (int is=0; is<nshots; ++is) {

			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ng",0)));
			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX()));
			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<Down> down (new Down(_slow,param,ref));
			std::shared_ptr<Up> up (new Up (_slow,param,ref));
			std::shared_ptr<IC> ic = std::make_shared<IC>(_bg_wfld);

			inj_rec->setShot(is);
			inj_src->setShot(is);

				// setting up
				inj_src->setStep(index[i]);
				inj_rec->setStep(index[i]);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(getWavelet(),_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// linearizetion of reflectivity
				ic->forward(modOneFreq,_sc_wfld,0);

				up->forward(_sc_wfld,_sc_wfld,1);

				inj_rec->adjoint(data,_sc_wfld,1);
			}
		// });
		}
	});
}



void BornRefl::adj_propagate(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {

	int nf = model->getHyper()->getAxis(3).n;

	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
		[=] (const tbb::blocked_range<int> &r) {

			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<Down> down (new Down(_slow,param,ref));
			std::shared_ptr<Up> up (new Up (_slow,param,ref));
			std::shared_ptr<IC> ic = std::make_shared<IC>(_bg_wfld);

		for (int i=r.begin(); i<r.end(); i++) {

			std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(_slow->getHyper());
			std::shared_ptr<Illumination> illum = std::make_shared<Illumination>(_slow->getHyper());

			// tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
			// 	[=](const tbb::blocked_range<int> &rs) {

	for (int is=0; is<nshots; ++is) {

			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ng",0)));
			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX()));

			inj_rec->setShot(is);
			inj_src->setShot(is);

				// setting up
				inj_src->setStep(index[i]);
				inj_rec->setStep(index[i]);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(getWavelet(),_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// illum->accumulate(_bg_wfld);

				inj_rec->forward(data,_sc_wfld,0);
				up->adjoint(_sc_wfld,_sc_wfld,1);
				ic->adjoint(modOneFreq, _sc_wfld,1);
			}
		// });
		// illum->compensate(modOneFreq);
		sliceModel(modOneFreq,model,i);
		}
	});
}

void BornRefl::sliceModel(std::shared_ptr<complex3DReg> input,std::shared_ptr<complex2DReg> output, int ifreq) {
	std::copy_n(input->_mat->data()+ifreq*nz*nx,nz*nx,output->_mat->data());
}
void BornRefl::sliceModel(std::shared_ptr<complex2DReg> input,std::shared_ptr<complex3DReg> output, int ifreq) {
	std::transform(input->_mat->data(),input->_mat->data()+nx*nz,
									output->_mat->data()+ifreq*nz*nx, output->_mat->data()+ifreq*nz*nx,
									[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );
}
