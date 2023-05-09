#include "extWEM.h"
#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

using namespace SEP;

void extWEM::full_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data) {

	auto slow_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
	pad_slow->forward(model[0],slow_pad,0);
	auto den_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
	pad_slow->forward(model[1],den_pad,0);

	tbb::parallel_for(tbb::blocked_range2d<int>(0,getNshots(),0,freq.size()),
		[=](const tbb::blocked_range2d<int> &r) {

		std::vector<std::shared_ptr<complex2DReg>> modOneFreq;
		// slowness
		modOneFreq.push_back(std::make_shared<complex2DReg>(slow_pad->getHyper()->getAxis(1),slow_pad->getHyper()->getAxis(2)));
		// density
		modOneFreq.push_back(std::make_shared<complex2DReg>(slow_pad->getHyper()->getAxis(1),slow_pad->getHyper()->getAxis(2)));
		std::shared_ptr<complex2DReg> _wfld = std::make_shared<complex2DReg>(modOneFreq[0]->getHyper());
		
	for (int i=r.cols().begin(); i < r.cols().end(); i++) {
		sliceModel(slow_pad,modOneFreq[0],i);
		sliceModel(den_pad,modOneFreq[1],i);
		std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(modOneFreq[0],nref);
		std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(modOneFreq);
		std::shared_ptr<Down> down = std::make_shared<Down>(modOneFreq[0],param,ref);
		std::shared_ptr<Up> up = std::make_shared<Up>(modOneFreq[0],param,ref);

		for (int is=r.rows().begin(); is < r.rows().end(); ++is) {
			std::shared_ptr<Injection> inj_src = std::make_shared<Injection> (getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0));
			std::shared_ptr<Injection> inj_rec = std::make_shared<Injection>(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0));

				inj_src->setShot(is);
				inj_rec->setShot(is);

					inj_src->setStep(index[i]);
					inj_src->forward(_wave_f,_wfld,0);

					down->setFreq(freq[i]);
					down->forward(_wfld,_wfld,1);
					inj_rec->setStep(index[i]);
					inj_rec->adjoint(data,_wfld,1);

					reflect->forward(_wfld,_wfld,1);

					up->setFreq(freq[i]);
					up->forward(_wfld,_wfld,1);

					inj_rec->setStep(index[i]);
					inj_rec->adjoint(data,_wfld,1);
				}
			}
		},tbb::static_partitioner());
	// });
}

void extWEM::full_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex4DReg> data) {

	// int nf = model->getHyper()->getAxis(3).n;
	// auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
	// pad_model->forward(model,model_pad,0);

	// tbb::parallel_for(tbb::blocked_range<int>(0,freq.size()),
	// 	[=](const tbb::blocked_range<int> &r) {

	// // tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
	// // 	[=](const tbb::blocked_range<int> &rs) {

	// 	std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));

	// 	std::shared_ptr<complex2DReg> _wfld = std::make_shared<complex2DReg>(modOneFreq->getHyper());

	// for (int i=r.begin(); i < r.end(); i++) {
	// 	sliceModel(model_pad,modOneFreq,i);
	// 	std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(modOneFreq,nref);
	// 	std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(modOneFreq);
	// 	std::shared_ptr<Down> down = std::make_shared<Down>(modOneFreq,param,ref);
	// 	std::shared_ptr<Up> up = std::make_shared<Up>(modOneFreq,param,ref);

	// 	for (int is=0; is < nshots; ++is) {
	// 		std::shared_ptr<Injection> inj_src = std::make_shared<Injection> (getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("tap",0));
	// 		std::shared_ptr<Injection> inj_rec = std::make_shared<Injection>(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("tap",0));


	// 			inj_src->setShot(is);
	// 			inj_rec->setShot(is);

	// 				inj_src->setStep(index[i]);
	// 				inj_rec->setStep(index[i]);

	// 				inj_src->forward(_wave_f,_wfld,0);

	// 				down->setFreq(freq[i]);
	// 				down->forward(_wfld,_wfld,1);
	// 				inj_rec->adjoint(data,_wfld,1);
	// 				//
	// 				// reflect->forward(_wfld,_wfld,1);
	// 				//
	// 				// up->setFreq(freq[i]);
	// 				// up->forward(_wfld,_wfld,1);
	// 				//
	// 				// inj_rec->adjoint(data,_wfld,1);
	// 			}
	// 		}
	// 	// });
	// });
}

void extWEM::onepass_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data) {

	// int nf = model->getHyper()->getAxis(3).n;
	// auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
	// pad_slow->forward(model,model_pad,0);

	// tbb::parallel_for(tbb::blocked_range2d<int>(0,getNshots(),0,freq.size()),
	// 	[=](const tbb::blocked_range2d<int> &r) {

	// // tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
	// // 	[=](const tbb::blocked_range<int> &rs) {

	// 	std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));

	// 	std::shared_ptr<complex2DReg> _wfld = std::make_shared<complex2DReg>(modOneFreq->getHyper());

	// for (int i=r.cols().begin(); i < r.cols().end(); i++) {
		
	// 	sliceModel(model_pad,modOneFreq,i);
	// 	std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(modOneFreq,nref);
	// 	std::shared_ptr<DownSource> down_src = std::make_shared<DownSource>(modOneFreq,param,ref);

	// 	for (int is=r.rows().begin(); is < r.rows().end(); ++is) {
	// 		std::shared_ptr<Injection> inj_src = std::make_shared<Injection> (getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0));
	// 		std::shared_ptr<Injection> inj_rec = std::make_shared<Injection>(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0));

	// 			inj_src->setShot(is);
	// 			inj_rec->setShot(is);

	// 				inj_src->setStep(index[i]);
	// 				inj_src->forward(_wave_f,_wfld,0);

	// 				inj_src->setStep(index[i]);
	// 				inj_src->forward(_wave_f,_wfld,0);
	// 				down_src->setFreq(freq[i]);
	// 				down_src->forward(_wfld,_wfld,1);

	// 				inj_rec->setStep(index[i]);
	// 				inj_rec->adjoint(data,_wfld,1);
	// 			}
	// 		}
	// 	},tbb::static_partitioner());
	// // });
}

void extWEM::partial_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data) {

	// int nf = model->getHyper()->getAxis(3).n;
	// auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
	// pad_model->forward(model,model_pad,0);

	// auto mod2d = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
	// stack->forward(model_pad,mod2d,false);

	// tbb::parallel_for(tbb::blocked_range<int>(0,freq.size()),
	// 	[&](const tbb::blocked_range<int> &r) {

	// // tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
	// // 	[=](const tbb::blocked_range<int> &rs) {

	// 	std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));

	// 	std::shared_ptr<complex2DReg> _wfld = std::make_shared<complex2DReg>(modOneFreq->getHyper());

	// for (int i=r.begin(); i < r.end(); i++) {
	// 	sliceModel(model_pad,modOneFreq,i);
	// 	std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(mod2d,nref);
	// 	std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(modOneFreq);
	// 	std::shared_ptr<DownSource> down_src = std::make_shared<DownSource>(mod2d,param,ref);
	// 	std::shared_ptr<UpSource> up_src = std::make_shared<UpSource>(mod2d,param,ref);

	// 	for (int is=0; is < nshots; ++is) {
	// 		std::shared_ptr<Injection> inj_src = std::make_shared<Injection> (getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("tap",0));
	// 		std::shared_ptr<Injection> inj_rec = std::make_shared<Injection>(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("tap",0));


	// 			inj_src->setShot(is);
	// 			inj_rec->setShot(is);

	// 				inj_src->setStep(index[i]);
	// 				inj_src->forward(_wave_f,_wfld,0);

	// 				down_src->setFreq(freq[i]);
	// 				down_src->forward(_wfld,_wfld,1);
	// 				inj_rec->setStep(index[i]);
	// 				inj_rec->adjoint(data,_wfld,1);

	// 				reflect->forward(_wfld,_wfld,1);

	// 				up_src->setFreq(freq[i]);
	// 				up_src->forward(_wfld,_wfld,1);

	// 				inj_rec->setStep(index[i]);
	// 				inj_rec->adjoint(data,_wfld,1);
	// 			}
	// 		}
	// 	});
	// // });
}
