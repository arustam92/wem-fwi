#include "extBornFull.h"
#include "illumination.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

using namespace SEP;

void extBornFull::full_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data) {

	auto dslow_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
	pad_model->forward(model[0],dslow_pad,0);
	auto dden_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
	pad_model->forward(model[1],dden_pad,0);

	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
		[=] (const tbb::blocked_range<int> &r) {

		std::vector<std::shared_ptr<complex2DReg>> modOneFreq;
		std::vector<std::shared_ptr<complex2DReg>> bgOneFreq;
		// slowness
		modOneFreq.push_back(std::make_shared<complex2DReg>(dslow_pad->getHyper()->getAxis(1),dslow_pad->getHyper()->getAxis(2)));
		bgOneFreq.push_back(std::make_shared<complex2DReg>(dslow_pad->getHyper()->getAxis(1),dslow_pad->getHyper()->getAxis(2)));
		// density
		modOneFreq.push_back(std::make_shared<complex2DReg>(dslow_pad->getHyper()->getAxis(1),dslow_pad->getHyper()->getAxis(2)));
		bgOneFreq.push_back(std::make_shared<complex2DReg>(dslow_pad->getHyper()->getAxis(1),dslow_pad->getHyper()->getAxis(2)));
		
		std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0)));
		std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));
		std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));

		for (int i=r.begin(); i<r.end(); i++) {

			sliceModel(getBgSlow(),bgOneFreq[0],i);
			sliceModel(getBgDensity(),bgOneFreq[1],i);
			std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(bgOneFreq[0],nref);
			std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(bgOneFreq);
			sliceModel(dslow_pad,modOneFreq[0],i);
			sliceModel(dden_pad,modOneFreq[1],i);

			std::shared_ptr<OneWay> down (new Down(bgOneFreq[0],param,ref));
			std::shared_ptr<OneWay> up (new Up (bgOneFreq[0],param,ref));

	// tbb::parallel_for(tbb::blocked_range<int>(0,getNshots()),
	// 	[=](const tbb::blocked_range<int> &rs) {

	for (int is=0; is<getNshots(); ++is) {

			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0)));

			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));
			std::unique_ptr<LinDown> lin_down(new LinDown(bgOneFreq[0],param,_bg_wfld,down));
			std::unique_ptr<LinUp> lin_up(new LinUp(bgOneFreq[0],param,_bg_wfld,up));

			inj_rec->setShot(is);
			inj_src->setShot(is);

				// setting up
				inj_src->setStep(index[i]);
				inj_rec->setStep(index[i]);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				lin_down->setFreq(freq[i]);
				lin_up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(getWavelet(),_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// forward scattering
				lin_down->forward(modOneFreq[0],_sc_wfld,0);
				down->forward(_sc_wfld,_sc_wfld,1);
				inj_rec->adjoint(data,_sc_wfld,1);

				// linearizetion of reflectivity
				born_refl->forward(modOneFreq,data,1);

				// down-scattering
				lin_down->forward(modOneFreq[0],_sc_wfld,0);
				down->forward(_sc_wfld,_sc_wfld,1);
				reflect->forward(_sc_wfld,_sc_wfld,1);
				up->forward(_sc_wfld,_sc_wfld,1);
				inj_rec->adjoint(data,_sc_wfld,1);

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				// up-scattering
				lin_up->forward(modOneFreq[0],_sc_wfld,0);
				up->forward(_sc_wfld,_sc_wfld,1);
				inj_rec->adjoint(data,_sc_wfld,1);
			}
		// });
		}
	},tbb::static_partitioner());
}

void extBornFull::full_adj(std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data) {

	auto dslow_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
	auto dden_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());

	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
		[=] (const tbb::blocked_range<int> &r) {

			std::vector<std::shared_ptr<complex2DReg>> modOneFreq;
			std::vector<std::shared_ptr<complex2DReg>> bgOneFreq;
			// slowness
			modOneFreq.push_back(std::make_shared<complex2DReg>(dslow_pad->getHyper()->getAxis(1),dslow_pad->getHyper()->getAxis(2)));
			bgOneFreq.push_back(std::make_shared<complex2DReg>(dslow_pad->getHyper()->getAxis(1),dslow_pad->getHyper()->getAxis(2)));
			// density
			modOneFreq.push_back(std::make_shared<complex2DReg>(dslow_pad->getHyper()->getAxis(1),dslow_pad->getHyper()->getAxis(2)));
			bgOneFreq.push_back(std::make_shared<complex2DReg>(dslow_pad->getHyper()->getAxis(1),dslow_pad->getHyper()->getAxis(2)));
		
			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0)));
			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));

		for (int i=r.begin(); i<r.end(); i++) {

			std::for_each(modOneFreq.begin(), modOneFreq.end(), [](std::shared_ptr<complex2DReg>& v) { v->scale(0); });
			sliceModel(getBgSlow(),bgOneFreq[0],i);
			sliceModel(getBgDensity(),bgOneFreq[1],i);
			std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(bgOneFreq[0],nref);
			std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(bgOneFreq);

			std::shared_ptr<OneWay> down (new Down(bgOneFreq[0],param,ref));
			std::shared_ptr<OneWay> up (new Up (bgOneFreq[0],param,ref));

			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));
			std::unique_ptr<LinDown> lin_down(new LinDown(bgOneFreq[0],param,_bg_wfld,down));
			std::unique_ptr<LinUp> lin_up(new LinUp(bgOneFreq[0],param,_bg_wfld,up));

			// tbb::parallel_for(tbb::blocked_range<int>(0,getNshots()),
			// 	[=](const tbb::blocked_range<int> &rs) {

	for (int is=0; is<getNshots(); ++is) {

			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0)));


			inj_rec->setShot(is);
			inj_src->setShot(is);

				// setting up
				inj_src->setStep(index[i]);
				inj_rec->setStep(index[i]);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				lin_down->setFreq(freq[i]);
				lin_up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(getWavelet(),_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// forward scattering
				inj_rec->forward(data,_sc_wfld,0);
				down->adjoint(_sc_wfld,_sc_wfld,1);
				lin_down->adjoint(modOneFreq[0],_sc_wfld,1);

				// linearizetion of reflectivity
				born_refl->adjoint(modOneFreq,data,1);

				// down-scattering
				inj_rec->forward(data,_sc_wfld,0);
				up->adjoint(_sc_wfld,_sc_wfld,1);
				// cmap[is]->forward(_sc_wfld, _sc_up_wfld);
				reflect->adjoint(_sc_wfld,_sc_wfld,1);
				down->adjoint(_sc_wfld,_sc_wfld,1);
				lin_down->adjoint(modOneFreq[0],_sc_wfld,1);

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				// illum->accumulate(_bg_wfld);
				// up-scattering
				inj_rec->forward(data,_sc_wfld,0);
				up->adjoint(_sc_wfld,_sc_wfld,1);
				lin_up->adjoint(modOneFreq[0],_sc_wfld,1);
			}
		// });
		// illum->compensate(modOneFreq);
		sliceModel(modOneFreq[0],dslow_pad,i);
		sliceModel(modOneFreq[1],dden_pad,i);
		}
	},tbb::static_partitioner());
	// stack->adjoint(model_pad,mod2d,1);
	pad_model->adjoint(model[0],dslow_pad,1);
	pad_model->adjoint(model[1],dden_pad,1);
}

// void extBornFull::onepass_fwd(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {

// 	int nf = model->getHyper()->getAxis(3).n;
// 	auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
// 	pad_model->forward(model,model_pad,0);

// 	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
// 		[=] (const tbb::blocked_range<int> &r) {

// 			std::shared_ptr<complex2DReg> slowOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
// 			std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2));
// 			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0)));
// 			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));
// 			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));

// 		for (int i=r.begin(); i<r.end(); i++) {
			
// 			sliceModel(getBgSlow(),slowOneFreq,i);
// 			std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(slowOneFreq,nref);
// 			sliceModel(model_pad,modOneFreq,i);

// 			std::shared_ptr<Down> down (new Down(slowOneFreq,param,ref));
// 			std::shared_ptr<DownSource> down_src (new DownSource(slowOneFreq,param,ref));
// 	// tbb::parallel_for(tbb::blocked_range<int>(0,getNshots()),
// 	// 	[=](const tbb::blocked_range<int> &rs) {

// 	for (int is=0; is<getNshots(); ++is) {

// 			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0)));
// 			std::unique_ptr<LinDown> lin_down(new LinDown(slowOneFreq,param,_bg_wfld,down));

// 			inj_rec->setShot(is);
// 			inj_src->setShot(is);

// 				// setting up
// 				inj_src->setStep(index[i]);
// 				inj_rec->setStep(index[i]);
// 				down->setFreq(freq[i]);
// 				down_src->setFreq(freq[i]);
// 				lin_down->setFreq(freq[i]);

// 				// down background wavefield
// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				lin_down->forward(modOneFreq,_sc_wfld,0);
// 				down->forward(_sc_wfld,_sc_wfld,1);
// 				inj_rec->adjoint(data,_sc_wfld,1);

// 				// lin_up->forward(modOneFreq,_sc_wfld,0);
// 				// up->forward(_sc_wfld,_sc_wfld,1);
// 				// inj_rec->adjoint(data,_sc_wfld,1);
// 			}
// 		// });
// 		}
// 	},tbb::static_partitioner());
// }

// void extBornFull::onepass_adj(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {

// 	int nf = model->getHyper()->getAxis(3).n;
// 	auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());

// 	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
// 		[=] (const tbb::blocked_range<int> &r) {

// 			std::shared_ptr<complex2DReg> slowOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
// 			std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(slowOneFreq->getHyper());

// 			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0)));
// 			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));
// 			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));

// 		for (int i=r.begin(); i<r.end(); i++) {

// 			modOneFreq->scale(0);
// 			sliceModel(getBgSlow(),slowOneFreq,i);
// 			std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(slowOneFreq,nref);

// 			std::shared_ptr<Down> down (new Down(slowOneFreq,param,ref));
// 			std::shared_ptr<DownSource> down_src (new DownSource(slowOneFreq,param,ref));
// 			std::unique_ptr<LinDown> lin_down(new LinDown(slowOneFreq,param,_bg_wfld,down));

// 			// tbb::parallel_for(tbb::blocked_range<int>(0,getNshots()),
// 			// 	[=](const tbb::blocked_range<int> &rs) {

// 	for (int is=0; is<getNshots(); ++is) {

// 			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0)));


// 			inj_rec->setShot(is);
// 			inj_src->setShot(is);

// 				// setting up
// 				inj_src->setStep(index[i]);
// 				inj_rec->setStep(index[i]);
// 				down->setFreq(freq[i]);
// 				down_src->setFreq(freq[i]);
// 				lin_down->setFreq(freq[i]);

// 				// down background wavefield
// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				// down-scattering
// 				inj_rec->forward(data,_sc_wfld,0);
// 				down->adjoint(_sc_wfld,_sc_wfld,1);
// 				lin_down->adjoint(modOneFreq,_sc_wfld,1);

// 				// inj_rec->forward(data,_sc_wfld,0);
// 				// up->adjoint(_sc_wfld,_sc_wfld,1);
// 				// lin_up->adjoint(modOneFreq,_sc_wfld,1);

// 			}
// 		// });
// 		// illum->compensate(modOneFreq);
// 		sliceModel(modOneFreq,model_pad,i);
// 		}
// 	},tbb::static_partitioner());
// 	pad_model->adjoint(model,model_pad,1);
// }

// void extBornFull::partial_fwd(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {

// 	int nf = model->getHyper()->getAxis(3).n;
// 	auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
// 	pad_model->forward(model,model_pad,0);

// 	auto mod2d = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
// 	stack->forward(model_pad,mod2d,false);

// 	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
// 		[&] (const tbb::blocked_range<int> &r) {

// 			std::shared_ptr<complex2DReg> slowOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
// 			std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2));
// 			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0)));
// 			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));
// 			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));

// 		for (int i=r.begin(); i<r.end(); i++) {

// 			sliceModel(getBgSlow(),slowOneFreq,i);
// 			std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(_slow2d,nref);
// 			std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(slowOneFreq);
// 			sliceModel(model_pad,modOneFreq,i);

// 			std::shared_ptr<DownSource> down_src (new DownSource(_slow2d,param,ref));
// 			std::shared_ptr<OneWay> down (new Down(_slow2d,param,ref));
// 			std::shared_ptr<OneWay> up (new Up (_slow2d,param,ref));
// 			std::shared_ptr<UpSource> up_src (new UpSource(_slow2d,param,ref));

// 	// tbb::parallel_for(tbb::blocked_range<int>(0,getNshots()),
// 	// 	[=](const tbb::blocked_range<int> &rs) {

// 	for (int is=0; is<getNshots(); ++is) {

// 			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0)));

// 			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));
// 			std::unique_ptr<LinDown> lin_down(new LinDown(_slow2d,param,_bg_wfld,down));
// 			std::unique_ptr<LinUp> lin_up(new LinUp(_slow2d,param,_bg_wfld,up));

// 			inj_rec->setShot(is);
// 			inj_src->setShot(is);

// 				// setting up
// 				inj_src->setStep(index[i]);
// 				inj_rec->setStep(index[i]);
// 				down->setFreq(freq[i]);
// 				up->setFreq(freq[i]);
// 				down_src->setFreq(freq[i]);
// 				up_src->setFreq(freq[i]);
// 				lin_down->setFreq(freq[i]);
// 				lin_up->setFreq(freq[i]);

// 				// down background wavefield
// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				// forward scattering
// 				lin_down->forward(mod2d,_sc_wfld,0);
// 				down->forward(_sc_wfld,_sc_wfld,1);
// 				inj_rec->adjoint(data,_sc_wfld,1);

// 				// linearizetion of reflectivity
// 				born_refl->forward(modOneFreq,data,1);

// 				// down-scattering
// 				lin_down->forward(mod2d,_sc_wfld,0);
// 				down->forward(_sc_wfld,_sc_wfld,1);
// 				reflect->forward(_sc_wfld,_sc_wfld,1);
// 				up->forward(_sc_wfld,_sc_wfld,1);
// 				inj_rec->adjoint(data,_sc_wfld,1);

// 				// up background wavefield
// 				reflect->forward(_bg_wfld,_bg_wfld,1);
// 				up_src->forward(_bg_wfld,_bg_wfld,1);

// 				// up-scattering
// 				lin_up->forward(mod2d,_sc_wfld,0);
// 				up->forward(_sc_wfld,_sc_wfld,1);
// 				inj_rec->adjoint(data,_sc_wfld,1);
// 			}
// 		// });
// 		}
// 	});
// }

// void extBornFull::partial_adj(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {

// 	int nf = model->getHyper()->getAxis(3).n;
// 	auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());

// 	std::shared_ptr<complex2DReg> mod2d = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));

// 	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
// 		[&] (const tbb::blocked_range<int> &r) {

// 			std::shared_ptr<complex2DReg> slowOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
// 			std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(slowOneFreq->getHyper());

// 			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0)));
// 			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));
// 			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));

// 		for (int i=r.begin(); i<r.end(); i++) {

// 			sliceModel(getBgSlow(),slowOneFreq,i);
// 			std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(_slow2d,nref);
// 			std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(slowOneFreq);

// 			std::shared_ptr<DownSource> down_src (new DownSource(_slow2d,param,ref));
// 			std::shared_ptr<OneWay> down (new Down(_slow2d,param,ref));
// 			std::shared_ptr<OneWay> up (new Up (_slow2d,param,ref));
// 			std::shared_ptr<UpSource> up_src (new UpSource(_slow2d,param,ref));

// 			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));
// 			std::unique_ptr<LinDown> lin_down(new LinDown(_slow2d,param,_bg_wfld,down));
// 			std::unique_ptr<LinUp> lin_up(new LinUp(_slow2d,param,_bg_wfld,up));

// 			// tbb::parallel_for(tbb::blocked_range<int>(0,getNshots()),
// 			// 	[=](const tbb::blocked_range<int> &rs) {

// 	for (int is=0; is<getNshots(); ++is) {

// 			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0)));


// 			inj_rec->setShot(is);
// 			inj_src->setShot(is);

// 				// setting up
// 				inj_src->setStep(index[i]);
// 				inj_rec->setStep(index[i]);
// 				down->setFreq(freq[i]);
// 				up->setFreq(freq[i]);
// 				down_src->setFreq(freq[i]);
// 				up_src->setFreq(freq[i]);

// 				lin_down->setFreq(freq[i]);
// 				lin_up->setFreq(freq[i]);

// 				// down background wavefield
// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				// forward scattering
// 				inj_rec->forward(data,_sc_wfld,0);
// 				down->adjoint(_sc_wfld,_sc_wfld,1);
// 				lin_down->adjoint(mod2d,_sc_wfld,1);

// 				// linearizetion of reflectivity
// 				born_refl->adjoint(modOneFreq,data,1);

// 				// down-scattering
// 				inj_rec->forward(data,_sc_wfld,0);
// 				up->adjoint(_sc_wfld,_sc_wfld,1);
// 				// cmap[is]->forward(_sc_wfld, _sc_up_wfld);
// 				reflect->adjoint(_sc_wfld,_sc_wfld,1);
// 				down->adjoint(_sc_wfld,_sc_wfld,1);
// 				lin_down->adjoint(mod2d,_sc_wfld,1);

// 				// up background wavefield
// 				reflect->forward(_bg_wfld,_bg_wfld,1);
// 				up_src->forward(_bg_wfld,_bg_wfld,1);

// 				// illum->accumulate(_bg_wfld);
// 				// up-scattering
// 				inj_rec->forward(data,_sc_wfld,0);
// 				up->adjoint(_sc_wfld,_sc_wfld,1);
// 				lin_up->adjoint(mod2d,_sc_wfld,1);
// 			}
// 		// });
// 		// illum->compensate(modOneFreq);
// 		sliceModel(modOneFreq,model_pad,i);
// 		}
// 	});
// 	stack->adjoint(model_pad,mod2d,1);
// 	pad_model->adjoint(model,model_pad,1);
// }

// void extBornFull::dru_fwd(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {

// 	int nf = model->getHyper()->getAxis(3).n;
// 	auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());
// 	pad_model->forward(model,model_pad,0);

// 	auto mod2d = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
// 	// stack->forward(model_pad,mod2d,false);

// 	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
// 		[=] (const tbb::blocked_range<int> &r) {

// 			std::shared_ptr<complex2DReg> slowOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
// 			std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2));
// 			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0)));
// 			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));
// 			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));

// 		for (int i=r.begin(); i<r.end(); i++) {

// 			sliceModel(getBgSlow(),slowOneFreq,i);
// 			std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(slowOneFreq,nref);
// 			std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(slowOneFreq);
// 			sliceModel(model_pad,modOneFreq,i);

// 			std::shared_ptr<OneWay> down (new Down(slowOneFreq,param,ref));
// 			std::shared_ptr<DownSource> down_src (new DownSource(slowOneFreq,param,ref));
// 			std::shared_ptr<OneWay> up (new Up (slowOneFreq,param,ref));
// 			std::shared_ptr<UpSource> up_src (new UpSource(slowOneFreq,param,ref));

// 	// tbb::parallel_for(tbb::blocked_range<int>(0,getNshots()),
// 	// 	[=](const tbb::blocked_range<int> &rs) {

// 	for (int is=0; is<getNshots(); ++is) {

// 			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0)));
// 			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,down,_bg_wfld,_sc_wfld));
// 			std::unique_ptr<LinDown> lin_down(new LinDown(slowOneFreq,param,_bg_wfld,down));
// 			std::unique_ptr<LinUp> lin_up(new LinUp(slowOneFreq,param,_bg_wfld,up));

// 			inj_rec->setShot(is);
// 			inj_src->setShot(is);

// 				// setting up
// 				inj_src->setStep(index[i]);
// 				inj_rec->setStep(index[i]);
// 				down_src->setFreq(freq[i]);
// 				up_src->setFreq(freq[i]);
// 				down->setFreq(freq[i]);
// 				up->setFreq(freq[i]);

// 				lin_down->setFreq(freq[i]);
// 				lin_up->setFreq(freq[i]);

// 				// down background wavefield
// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				// // forward scattering
// 				lin_down->forward(modOneFreq,_sc_wfld,0);
// 				down->forward(_sc_wfld,_sc_wfld,1);
// 				inj_rec->adjoint(data,_sc_wfld,1);

// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				up_src->forward(_bg_wfld,_bg_wfld,1);

// 				// linearizetion of reflectivity
// 				born_refl->forward(modOneFreq,data,1);

// 				// up-scattering
// 				lin_up->forward(modOneFreq,_sc_wfld,0);
// 				up->forward(_sc_wfld,_sc_wfld,1);
// 				reflect->forward(_sc_wfld,_sc_wfld,1);
// 				down->forward(_sc_wfld,_sc_wfld,1);
// 				inj_rec->adjoint(data,_sc_wfld,1);

// 				// down background wavefield
// 				reflect->forward(_bg_wfld,_bg_wfld,1);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				// down-scattering
// 				lin_down->forward(modOneFreq,_sc_wfld,0);
// 				down->forward(_sc_wfld,_sc_wfld,1);
// 				inj_rec->adjoint(data,_sc_wfld,1);
// 			}
// 		// });
// 		}
// 	},tbb::static_partitioner());
// }

// void extBornFull::dru_adj(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {

// 	int nf = model->getHyper()->getAxis(3).n;
// 	auto model_pad = std::make_shared<complex3DReg>(pad_model->getPaddedHyper());

// 	std::shared_ptr<complex2DReg> mod2d = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));

// 	tbb::parallel_for(tbb::blocked_range<int> (0,freq.size()),
// 		[=] (const tbb::blocked_range<int> &r) {

// 			std::shared_ptr<complex2DReg> slowOneFreq = std::make_shared<complex2DReg>(model_pad->getHyper()->getAxis(1),model_pad->getHyper()->getAxis(2));
// 			std::shared_ptr<complex2DReg> modOneFreq = std::make_shared<complex2DReg>(slowOneFreq->getHyper());

// 			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX(),param->getInt("ngr",0),param->getInt("pad",0)));
// 			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));
// 			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(getBgSlow()->getHyper()->getAxis(1),getBgSlow()->getHyper()->getAxis(2)));

// 		for (int i=r.begin(); i<r.end(); i++) {

// 			modOneFreq->scale(0);
// 			sliceModel(getBgSlow(),slowOneFreq,i);
// 			std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(slowOneFreq,nref);
// 			std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(slowOneFreq);

// 			std::shared_ptr<OneWay> down (new Down(slowOneFreq,param,ref));
// 			std::shared_ptr<DownSource> down_src (new DownSource(slowOneFreq,param,ref));
// 			std::shared_ptr<OneWay> up (new Up (slowOneFreq,param,ref));
// 			std::shared_ptr<UpSource> up_src (new UpSource(slowOneFreq,param,ref));

// 			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,down,_bg_wfld,_sc_wfld));
// 			std::unique_ptr<LinDown> lin_down(new LinDown(slowOneFreq,param,_bg_wfld,down));
// 			std::unique_ptr<LinUp> lin_up(new LinUp(slowOneFreq,param,_bg_wfld,up));

// 			// tbb::parallel_for(tbb::blocked_range<int>(0,getNshots()),
// 			// 	[=](const tbb::blocked_range<int> &rs) {

// 	for (int is=0; is<getNshots(); ++is) {

// 			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ngs",0),param->getInt("pad",0)));


// 			inj_rec->setShot(is);
// 			inj_src->setShot(is);

// 				// setting up
// 				inj_src->setStep(index[i]);
// 				inj_rec->setStep(index[i]);
// 				down->setFreq(freq[i]);
// 				up->setFreq(freq[i]);
// 				down_src->setFreq(freq[i]);
// 				up_src->setFreq(freq[i]);

// 				lin_down->setFreq(freq[i]);
// 				lin_up->setFreq(freq[i]);

// 				// down background wavefield
// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				// forward scattering
// 				inj_rec->forward(data,_sc_wfld,0);
// 				down->adjoint(_sc_wfld,_sc_wfld,1);
// 				lin_down->adjoint(modOneFreq,_sc_wfld,1);

// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				up_src->forward(_bg_wfld,_bg_wfld,1);

// 				// linearizetion of reflectivity
// 				born_refl->adjoint(modOneFreq,data,1);

// 				// up-scattering
// 				inj_rec->forward(data,_sc_wfld,0);
// 				down->adjoint(_sc_wfld,_sc_wfld,1);
// 				reflect->adjoint(_sc_wfld,_sc_wfld,1);
// 				up->adjoint(_sc_wfld,_sc_wfld,1);
// 				lin_up->adjoint(modOneFreq,_sc_wfld,1);

// 				// up background wavefield
// 				reflect->forward(_bg_wfld,_bg_wfld,1);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				// up-scattering
// 				inj_rec->forward(data,_sc_wfld,0);
// 				down->adjoint(_sc_wfld,_sc_wfld,1);
// 				lin_down->adjoint(modOneFreq,_sc_wfld,1);
// 			}
// 		// });
// 		// illum->compensate(modOneFreq);
// 		sliceModel(modOneFreq,model_pad,i);
// 		}
// 	},tbb::static_partitioner());
// 	// stack->adjoint(model_pad,mod2d,1);
// 	pad_model->adjoint(model,model_pad,1);
// }

// sukkk
