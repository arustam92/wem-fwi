#include "WEM.h"
#include "Pad.h"
#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

using namespace SEP;

void WEM::full_fwd(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data) {

	tbb::parallel_for(tbb::blocked_range2d<int>(0,getNshots(),0,freq.size()),
		[=](const tbb::blocked_range2d<int> &r) {

		std::shared_ptr<complex2DReg> _wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
		std::vector<std::shared_ptr<complex2DReg>> grid;
		grid.push_back(std::make_shared<complex2DReg>(model[0]->getHyper()));
		grid.push_back(std::make_shared<complex2DReg>(model[0]->getHyper()));

	for (int is=r.rows().begin(); is < r.rows().end(); ++is) {

		auto scoord = cmap[is]->forward(getSrc()->getZ(is),getSrc()->getX(is));
		auto rcoord = cmap[is]->forward(getRec()->getZ(),getRec()->getX());
		
		cmap[is]->forward(model[0],grid[0]);
		cmap[is]->forward(model[1],grid[1]);
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
		_wfld->setHyper(grid_pad[0]->getHyper());

		std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(grid_pad[0],nref);

		std::shared_ptr<Injection> inj_src = std::make_shared<Injection> (scoord->real,scoord->imag,param->getInt("ngs",0),param->getInt("pad",0));
		std::shared_ptr<Injection> inj_rec = std::make_shared<Injection>(rcoord->real,rcoord->imag,param->getInt("ngr",0),param->getInt("pad",0));

		std::shared_ptr<Down> down (new Down(grid_pad[0],param,ref));
		std::shared_ptr<Up> up (new Up(grid_pad[0],param,ref));

			inj_src->setShot(is);
			inj_rec->setShot(is);

			for (int i=r.cols().begin(); i < r.cols().end(); ++i) {

				inj_src->setStep(index[i]);
				inj_src->forward(_wave_f,_wfld,0);
				down->setFreq(freq[i]);
				down->forward(_wfld,_wfld,1);
				
				// cmap[is]->inverse(ph_wfld,_wfld);
				inj_rec->setStep(index[i]);
				inj_rec->adjoint(data,_wfld,1);
	
				reflect->forward(_wfld,_wfld,1);
				
				// cmap[is]->inverse(ph_wfld,_wfld);
				//
				up->setFreq(freq[i]);
				up->forward(_wfld,_wfld,1);

				// ph_wfld = cmap[is]->inverse(_wfld);
				inj_rec->setStep(index[i]);
				inj_rec->adjoint(data,_wfld,1);

			}
		}
	},tbb::static_partitioner());
	// });
}

// void WEM::onepass_fwd(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data) {

// 	tbb::parallel_for(tbb::blocked_range2d<int>(0,getNshots(),0,freq.size()),
// 		[=](const tbb::blocked_range2d<int> &r) {

// 		std::shared_ptr<complex2DReg> _wfld = std::make_shared<complex2DReg>(pad_model->getPaddedHyper());
// 		std::shared_ptr<complex2DReg> grid = std::make_shared<complex2DReg>(model->getHyper());

// 		for (int is=r.rows().begin(); is < r.rows().end(); ++is) {

// 		auto scoord = cmap[is]->forward(getSrc()->getZ(is),getSrc()->getX(is));
// 		auto rcoord = cmap[is]->forward(getRec()->getZ(),getRec()->getX());
		
// 		cmap[is]->forward(model, grid);
// 		int pad = param->getInt("pad",0);
// 		Pad padOp(grid->getHyper(), pad, pad, true);
// 		auto grid_pad = std::make_shared<complex2DReg>(padOp.getPaddedHyper());
// 		padOp.forward(grid, grid_pad, 0);
		
// 		cmap[is]->scaleByJacobian(grid_pad);
// 		cmap[is]->scaleByJacobian(grid_pad);
// 		_wfld->setHyper(grid_pad->getHyper());

// 		std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(grid_pad,nref);
// 		std::shared_ptr<Injection> inj_src = std::make_shared<Injection> (scoord->real,scoord->imag,param->getInt("ngs",0),param->getInt("pad",0));
// 		std::shared_ptr<Injection> inj_rec = std::make_shared<Injection>(rcoord->real,rcoord->imag,param->getInt("ngr",0),param->getInt("pad",0));
// 		std::shared_ptr<DownSource> down_src = std::make_shared<DownSource>(grid_pad,param,ref);
		
// 		inj_src->setShot(is);
// 		inj_rec->setShot(is);

// 			for (int i=r.cols().begin(); i < r.cols().end(); ++i) {
// 				inj_src->setStep(index[i]);
// 				inj_src->forward(_wave_f,_wfld,0);

// 				down_src->setFreq(freq[i]);
// 				down_src->forward(_wfld,_wfld,1);
				
// 				inj_rec->setStep(index[i]);
// 				inj_rec->adjoint(data,_wfld,1);

// 			}
// 		}
// 	},tbb::static_partitioner());
// }

// void WEM::full_fwd(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex4DReg> data) {

// 	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
// 		[=](const tbb::blocked_range<int> &rs) {

// 	tbb::parallel_for(tbb::blocked_range<int>(0,freq.size()),
// 		[=](const tbb::blocked_range<int> &r) {

// 	for (int is=rs.begin(); is < rs.end(); ++is) {

// 		auto scoord = cmap[is]->forward(getSrc()->getZ(is),getSrc()->getX(is));
// 		// auto rcoord = cmap[is]->forward(getRec()->getZ(),getRec()->getX());
// 		auto grid = cmap[is]->forward(model);
// 		std::shared_ptr<Reflect> reflect = std::make_shared<Reflect>(grid);

// 		cmap[is]->scaleByJacobian(grid);

// 		std::shared_ptr<RefSampler> ref = std::make_shared<RefSampler>(grid,nref);
// 		std::shared_ptr<RefSampler> ref2 = std::make_shared<RefSampler>(model,nref);

// 		std::shared_ptr<Injection> inj_src = std::make_shared<Injection> (scoord->real,scoord->imag,param->getInt("ngs",0));

// 		std::shared_ptr<Injection> inj_rec = std::make_shared<Injection>(getRec()->getZ(),getRec()->getX());

// 		std::shared_ptr<complex2DReg> _wfld = std::make_shared<complex2DReg>(grid->getHyper()->getAxis(1),grid->getHyper()->getAxis(2));
// 		std::shared_ptr<Down> down = std::make_shared<Down>(grid,param,ref);
// 		std::shared_ptr<Up> up = std::make_shared<Up>(model,param,ref2);

// 			inj_src->setShot(is);
// 			inj_rec->setShot(is);

// 			for (int i=r.begin(); i < r.end(); ++i) {

// 				inj_src->setStep(index[i]);
// 				inj_src->forward(_wave_f,_wfld,0);
// 				down->setFreq(freq[i]);
// 				down->forward(_wfld,_wfld,1);
// 				auto ph_wfld = cmap[is]->inverse(_wfld);
// 				inj_rec->setStep(index[i]);
// 				inj_rec->adjoint(data,ph_wfld,1);

// 				// reflect->forward(_wfld,_wfld,1);
// 				// ph_wfld = cmap[is]->inverse(_wfld);
// 				// // //
// 				// up->setFreq(freq[i]);
// 				// up->forward(ph_wfld,ph_wfld,1);
// 				// // ph_wfld = cmap[is]->inverse(_wfld);
// 				// inj_rec->setStep(index[i]);
// 				// inj_rec->adjoint(data,ph_wfld,1);

// 			}
// 		}
// 	});
// 	});
// }
