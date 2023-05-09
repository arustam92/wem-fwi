
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <BornRefl.h>
#include <chrono>

using namespace std::chrono;

using namespace SEP;

void BornRefl::forward (std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	if(_verb_) std::cerr << "Forward complex born (refl) ... " << std::endl;
		auto start = high_resolution_clock::now();
	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is<rs.end(); ++is) {

		tbb::parallel_for(tbb::blocked_range<int> (fmin,fmax),
			[=] (const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ng",0)));
			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX()));
			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<Down> down (new Down(_slow,param,ref));
			std::shared_ptr<Up> up (new Up (_slow,param,ref));

			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));

			inj_rec->setShot(is);
			inj_src->setShot(is);

			for (int i=r.begin(); i<r.end(); i++) {

				// setting up
				inj_src->setStep(i);
				inj_rec->setStep(i);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(wave_f,_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// linearizetion of reflectivity
				born_refl->forward(model,data,1);

			}
		});
	}

	});

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	if(_verb_) std::cerr << duration.count()/1e6 << " s" << std::endl;

}

void BornRefl::c_adj_illum (std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) model->scale(0.);

	if(_verb_) std::cerr << "Adjoint complex born (refl) illumination... " << std::endl;
		auto start = high_resolution_clock::now();

		std::shared_ptr<Illumination> illum (new Illumination (_slow->getHyper()));
		std::shared_ptr<complex2DReg> temp = model->clone();

	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is<rs.end(); ++is) {

		// temp->scale(0);

		tbb::parallel_for(tbb::blocked_range<int> (fmin,fmax),
			[=] (const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ng",0)));
			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX()));
			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<Down> down (new Down(_slow,param,ref));
			std::shared_ptr<Up> up (new Up (_slow,param,ref));

			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));

			inj_rec->setShot(is);
			inj_src->setShot(is);

			for (int i=r.begin(); i<r.end(); i++) {

				// setting up
				inj_src->setStep(i);
				inj_rec->setStep(i);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(wave_f,_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				illum->accumulate(_bg_wfld);

				// linearizetion of reflectivity
				born_refl->adjoint(temp,data,1);

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				illum->accumulate(_bg_wfld);

			}
		});
		}

	});

	illum->compensate(temp);
	model -> scaleAdd(temp,1.,1.);

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	if(_verb_)	std::cerr << duration.count()/1e6 << " s" << std::endl;
}

void BornRefl::c_adj_reg (std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) model->scale(0.);

	if(_verb_) std::cerr << "Adjoint complex born (refl) ... " << std::endl;
		auto start = high_resolution_clock::now();
	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is<rs.end(); ++is) {

		tbb::parallel_for(tbb::blocked_range<int> (fmin,fmax),
			[=] (const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(getSrc()->getZ(is),getSrc()->getX(is),param->getInt("ng",0)));
			std::shared_ptr<Injection> inj_rec (new Injection(getRec()->getZ(),getRec()->getX()));
			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()->getAxis(1),_slow->getHyper()->getAxis(2)));
			std::shared_ptr<Down> down (new Down(_slow,param,ref));
			std::shared_ptr<Up> up (new Up (_slow,param,ref));

			std::unique_ptr<Born_refl> born_refl(new Born_refl(reflect,inj_rec,up,_bg_wfld,_sc_wfld));

			inj_rec->setShot(is);
			inj_src->setShot(is);

			for (int i=r.begin(); i<r.end(); i++) {

				// setting up
				inj_src->setStep(i);
				inj_rec->setStep(i);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(wave_f,_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// linearizetion of reflectivity
				born_refl->adjoint(model,data,1);

			}
		});
	}
	});

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	if(_verb_)	std::cerr << duration.count()/1e6 << " s" << std::endl;
}
