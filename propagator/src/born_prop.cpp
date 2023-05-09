#include "BornFull.h"
#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/global_control.h>
#include <execution>

using namespace SEP;

void BornFull::full_fwd(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data) {
	
	// tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);

	tbb::parallel_for(tbb::blocked_range2d<int>(0,getNshots(),0,freq.size()),
		[=](const tbb::blocked_range2d<int> &r) {

		auto _bg_wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
		auto _sc_wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());

		std::vector<std::shared_ptr<complex2DReg>> grid, mod;
		grid.push_back(std::make_shared<complex2DReg>(getBgSlow()->getHyper()));
		grid.push_back(std::make_shared<complex2DReg>(getBgSlow()->getHyper()));

		mod.push_back(std::make_shared<complex2DReg>(getBgSlow()->getHyper()));
		mod.push_back(std::make_shared<complex2DReg>(getBgSlow()->getHyper()));

		for (int is=r.rows().begin(); is<r.rows().end(); ++is) {

		auto scoord = cmap[is]->forward(getSrc()->getZ(is),getSrc()->getX(is));
		auto rcoord = cmap[is]->forward(getRec()->getZ(),getRec()->getX());
		
		cmap[is]->forward(getBgSlow(), grid[0]);
		cmap[is]->forward(getBgDensity(), grid[1]);

		int pad = param->getInt("pad",0);
		Pad padOp(grid[0]->getHyper(), pad, pad, true);
		std::vector<std::shared_ptr<complex2DReg>> grid_pad;
		grid_pad.push_back(std::make_shared<complex2DReg>(padOp.getPaddedHyper()));
		grid_pad.push_back(std::make_shared<complex2DReg>(padOp.getPaddedHyper()));

		padOp.forward(grid[0], grid_pad[0], 0);
		padOp.forward(grid[1], grid_pad[1], 0);

		std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(grid_pad);
		cmap[is]->scaleByJacobian(grid_pad[0]);
		cmap[is]->scaleByJacobian(grid_pad[0]);
		std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(grid_pad[0],nref);

		cmap[is]->forward(model[0],mod[0]);
		cmap[is]->forward(model[1],mod[1]);
		Pad padModel(mod[0]->getHyper(), pad, pad, false);
		std::vector<std::shared_ptr<complex2DReg>> grid_mod;
		grid_mod.push_back(std::make_shared<complex2DReg>(padModel.getPaddedHyper()));
		grid_mod.push_back(std::make_shared<complex2DReg>(padModel.getPaddedHyper()));

		padModel.forward(mod[0], grid_mod[0], 0);
		padModel.forward(mod[1], grid_mod[1], 0);

		_bg_wfld->setHyper(grid_pad[0]->getHyper());
		_sc_wfld->setHyper(grid_pad[0]->getHyper());

		std::unique_ptr<Injection> inj_src (new Injection(scoord->real,scoord->imag,param->getInt("ngs",0),param->getInt("pad",0)));
		std::shared_ptr<Injection> inj_rec (new Injection(rcoord->real,rcoord->imag,param->getInt("ngr",0),param->getInt("pad",0)));
		std::shared_ptr<OneWay> down (new Down(grid_pad[0],param,ref));
		std::shared_ptr<OneWay> up = std::make_shared<Up>(grid_pad[0],param,ref);
		
		std::unique_ptr<LinDown> lin_down (new LinDown (grid_pad[0],param,_bg_wfld,down));
		std::unique_ptr<LinUp> lin_up(new LinUp(grid_pad[0],param,_bg_wfld,up));
		std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));
		
		inj_rec->setShot(is);
		inj_src->setShot(is);
		
			for (int i=r.cols().begin(); i<r.cols().end(); i++) {

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
				lin_down->forward(grid_mod[0],_sc_wfld,0);
				down->forward(_sc_wfld,_sc_wfld,1);
				inj_rec->adjoint(data,_sc_wfld,1);

				// linearizetion of reflectivity
				born_refl->forward(grid_mod,data,1);

				// down-scattering
				lin_down->forward(grid_mod[0],_sc_wfld,0);
				down->forward(_sc_wfld,_sc_wfld,1);
				reflect->forward(_sc_wfld,_sc_wfld,1);
				up->forward(_sc_wfld,_sc_wfld,1);
				inj_rec->adjoint(data,_sc_wfld,1);

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				// up-scattering
				lin_up->forward(grid_mod[0],_sc_wfld,0);
				up->forward(_sc_wfld,_sc_wfld,1);
				inj_rec->adjoint(data,_sc_wfld,1);
			}
		}
	},tbb::static_partitioner());
}


void BornFull::full_adj(std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data) {

	tbb::parallel_for(
		tbb::blocked_range2d<int>(0,getNshots(),0,freq.size()),
		[=](const tbb::blocked_range2d<int> &r) {

		auto _bg_wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
		auto _sc_wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
		auto mm = std::make_shared<complex2DReg>(model[0]->getHyper());

		std::vector<std::shared_ptr<complex2DReg>> grid, mod, grid_mod;
		grid.push_back(std::make_shared<complex2DReg>(getBgSlow()->getHyper()));
		grid.push_back(std::make_shared<complex2DReg>(getBgSlow()->getHyper()));

		mod.push_back(std::make_shared<complex2DReg>(getBgSlow()->getHyper()));
		mod.push_back(std::make_shared<complex2DReg>(getBgSlow()->getHyper()));

		for (int is=r.rows().begin(); is<r.rows().end(); ++is) {
			
		auto scoord = cmap[is]->forward(getSrc()->getZ(is),getSrc()->getX(is));
		auto rcoord = cmap[is]->forward(getRec()->getZ(),getRec()->getX());
		
		cmap[is]->forward(getBgSlow(), grid[0]);
		cmap[is]->forward(getBgDensity(), grid[1]);

		int pad = param->getInt("pad",0);
		Pad padOp(grid[0]->getHyper(), pad, pad, true);
		std::vector<std::shared_ptr<complex2DReg>> grid_pad, grid_mod;
		grid_pad.push_back(std::make_shared<complex2DReg>(padOp.getPaddedHyper()));
		grid_pad.push_back(std::make_shared<complex2DReg>(padOp.getPaddedHyper()));

		grid_mod.push_back(std::make_shared<complex2DReg>(padOp.getPaddedHyper()));
		grid_mod.push_back(std::make_shared<complex2DReg>(padOp.getPaddedHyper()));

		padOp.forward(grid[0], grid_pad[0], 0);
		padOp.forward(grid[1], grid_pad[1], 0);

		std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(grid_pad);
		cmap[is]->scaleByJacobian(grid_pad[0]);
		cmap[is]->scaleByJacobian(grid_pad[0]);
		std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(grid_pad[0],nref);

		_bg_wfld->setHyper(grid_pad[0]->getHyper());
		_sc_wfld->setHyper(grid_pad[0]->getHyper());
		mod[0]->setHyper(grid[0]->getHyper());
		mod[1]->setHyper(grid[0]->getHyper());

		std::unique_ptr<Injection> inj_src (new Injection(scoord->real,scoord->imag,param->getInt("ngs",0),param->getInt("pad",0)));
		std::shared_ptr<Injection> inj_rec (new Injection(rcoord->real,rcoord->imag,param->getInt("ngr",0),param->getInt("pad",0)));
		std::shared_ptr<OneWay> down (new Down(grid_pad[0],param,ref));
		std::shared_ptr<OneWay> up = std::make_shared<Up>(grid_pad[0],param,ref);
		
		std::unique_ptr<LinDown> lin_down (new LinDown (grid_pad[0],param,_bg_wfld,down));
		std::unique_ptr<LinUp> lin_up(new LinUp(grid_pad[0],param,_bg_wfld,up));
		std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));
		
		inj_rec->setShot(is);
		inj_src->setShot(is);

			for (int i=r.cols().begin(); i<r.cols().end(); i++) {

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

				// linearizetion of reflectivity
				born_refl->adjoint(grid_mod,data,1);

				// forward scattering
				inj_rec->forward(data,_sc_wfld,0);
				down->adjoint(_sc_wfld,_sc_wfld,1);
				lin_down->adjoint(grid_mod[0],_sc_wfld,1);

				// down-scattering
				inj_rec->forward(data,_sc_wfld,0);
				up->adjoint(_sc_wfld,_sc_wfld,1);
				reflect->adjoint(_sc_wfld,_sc_wfld,1);
				down->adjoint(_sc_wfld,_sc_wfld,1);
				lin_down->adjoint(grid_mod[0],_sc_wfld,1);

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				// up-scattering
				inj_rec->forward(data,_sc_wfld,0);
				up->adjoint(_sc_wfld,_sc_wfld,1);
				lin_up->adjoint(grid_mod[0],_sc_wfld,1);
			}
			padOp.adjoint(mod[0],grid_mod[0],0);
			padOp.adjoint(mod[1],grid_mod[1],0);

			cmap[is]->inverse(mm, mod[0]);
			std::transform(
									model[0]->getVals(), model[0]->getVals()+model[0]->getHyper()->getN123(), 
									mm->getVals(), model[0]->getVals(),
									[](std::complex<float> v1, std::complex<float> v2){
										return v1 + v2;});
			
			cmap[is]->inverse(mm, mod[1]);
			std::transform(
									model[1]->getVals(), model[1]->getVals()+model[1]->getHyper()->getN123(), 
									mm->getVals(), model[1]->getVals(),
									[](std::complex<float> v1, std::complex<float> v2){
										return v1 + v2;});
		}
		}, tbb::static_partitioner()) ;
}

// void BornFull::onepass_fwd(const std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data) {

// 	tbb::parallel_for(
// 		tbb::blocked_range2d<int>(0,getNshots(),0,freq.size()),
// 		[=](const tbb::blocked_range2d<int> &r) {

// 		auto _bg_wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
// 		auto _sc_wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
// 		auto grid_mod = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
// 		auto grid_slow = std::make_shared<complex2DReg>(_slow->getHyper());
// 		auto mod = std::make_shared<complex2DReg>(model->getHyper());

// 		for (int is=r.rows().begin(); is < r.rows().end(); ++is) {
		
// 		auto scoord = cmap[is]->forward(getSrc()->getZ(is),getSrc()->getX(is));
// 		auto rcoord = cmap[is]->forward(getRec()->getZ(),getRec()->getX());
	
// 		cmap[is]->forward(_slow, grid_slow);
// 		int pad = param->getInt("pad",0);
// 		Pad padOp(grid_slow->getHyper(), pad, pad, true);
// 		auto grid_pad = std::make_shared<complex2DReg>(padOp.getPaddedHyper());
// 		padOp.forward(grid_slow, grid_pad, 0);
		
// 		cmap[is]->scaleByJacobian(grid_pad);
// 		cmap[is]->scaleByJacobian(grid_pad);
// 		std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(grid_pad,nref);

// 		cmap[is]->forward(model,mod);
// 		grid_mod->setHyper(grid_slow->getHyper());
// 		padOp.forward(mod, grid_mod, 0);

// 		_bg_wfld->setHyper(grid_pad->getHyper());
// 		_sc_wfld->setHyper(grid_pad->getHyper());
		
// 		std::unique_ptr<Injection> inj_src (new Injection(scoord->real,scoord->imag,param->getInt("ngs",0),param->getInt("pad",0)));
// 		std::shared_ptr<Injection> inj_rec (new Injection(rcoord->real,rcoord->imag,param->getInt("ngr",0),param->getInt("pad",0)));
// 		std::shared_ptr<Down> down (new Down(grid_pad,param,ref));
// 		std::shared_ptr<DownSource> down_src (new DownSource(grid_pad,param,ref));
// 		std::unique_ptr<LinDown> lin_down (new LinDown (grid_pad,param,_bg_wfld,down));

// 			inj_rec->setShot(is);
// 			inj_src->setShot(is);

// 		for (int i=r.cols().begin(); i < r.cols().end(); ++i) {

// 				// setting up
// 				inj_src->setStep(index[i]);
// 				inj_rec->setStep(index[i]);
// 				down->setFreq(freq[i]);
// 				down_src->setFreq(freq[i]);
// 				lin_down->setFreq(freq[i]);

// 				// down background wavefield
// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);

// 				lin_down->forward(grid_mod,_sc_wfld,0);
// 				down->forward(_sc_wfld,_sc_wfld,1);

// 				inj_rec->adjoint(data,_sc_wfld,1);
// 			}
// 		}

// 	},tbb::static_partitioner());
// }


// void BornFull::onepass_adj(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data) {
// 	// tbb::global_control c(tbb::global_control::max_allowed_parallelism, 1);
// 	auto m0 = model->clone();
// 	m0->scale(0);

// 	tbb::parallel_for(
// 		tbb::blocked_range2d<int>(0,getNshots(),0,freq.size()),
// 		[=](const tbb::blocked_range2d<int> &r) {

// 			auto _bg_wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
// 			auto _sc_wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
// 			auto grid_mod = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
// 			auto grid_slow = std::make_shared<complex2DReg>(_slow->getHyper());
// 			auto mod = std::make_shared<complex2DReg>(model->getHyper());
// 			auto mm = std::make_shared<complex2DReg>(model->getHyper());

// 		for (int is=r.rows().begin(); is < r.rows().end(); ++is) {

// 		auto scoord = cmap[is]->forward(getSrc()->getZ(is),getSrc()->getX(is));
// 		auto rcoord = cmap[is]->forward(getRec()->getZ(),getRec()->getX());
		
// 		cmap[is]->forward(_slow, grid_slow);
// 		int pad = param->getInt("pad",0);
// 		Pad padOp(grid_slow->getHyper(), pad, pad, true);
// 		auto grid_pad = std::make_shared<complex2DReg>(padOp.getPaddedHyper());
// 		padOp.forward(grid_slow, grid_pad, 0);
	
// 		cmap[is]->scaleByJacobian(grid_pad);
// 		cmap[is]->scaleByJacobian(grid_pad);
// 		std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(grid_pad,nref);
		
// 		_bg_wfld->setHyper(grid_pad->getHyper());
// 		_sc_wfld->setHyper(grid_pad->getHyper());
// 		grid_mod->setHyper(grid_pad->getHyper());
// 		grid_mod->zero();
// 		mod->setHyper(grid_slow->getHyper());

// 		std::unique_ptr<Injection> inj_src (new Injection(scoord->real,scoord->imag,param->getInt("ngs",0),param->getInt("pad",0)));
// 		std::shared_ptr<Injection> inj_rec (new Injection(rcoord->real,rcoord->imag,param->getInt("ngr",0),param->getInt("pad",0)));
// 		std::shared_ptr<Down> down (new Down(grid_pad,param,ref));
// 		std::shared_ptr<DownSource> down_src (new DownSource(grid_pad,param,ref));
// 		std::unique_ptr<LinDown> lin_down (new LinDown (grid_pad,param,_bg_wfld,down));

// 			inj_rec->setShot(is);
// 			inj_src->setShot(is);

// 		for (int i=r.cols().begin(); i < r.cols().end(); ++i) {

// 				// setting up
// 				inj_src->setStep(index[i]);
// 				inj_rec->setStep(index[i]);
// 				down->setFreq(freq[i]);
// 				down_src->setFreq(freq[i]);
// 				lin_down->setFreq(freq[i]);

// 				// down background wavefield
// 				inj_src->forward(getWavelet(),_bg_wfld,0);
// 				down_src->forward(_bg_wfld,_bg_wfld,1);
// 				//
// 				inj_rec->forward(data,_sc_wfld,0);
// 				down->adjoint(_sc_wfld,_sc_wfld,1);
// 				lin_down->adjoint(grid_mod,_sc_wfld,1);
// 			}
// 			padOp.adjoint(mod,grid_mod,0);
// 			cmap[is]->inverse(mm,mod);
// 			std::transform(
// 									model->getVals(), model->getVals()+model->getHyper()->getN123(), 
// 									mm->getVals(), model->getVals(),
// 									[](std::complex<float> v1, std::complex<float> v2){
// 										return v1 + v2;});
// 		}
// 		}, tbb::static_partitioner()) ;
// }
