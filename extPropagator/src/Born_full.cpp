
// #include <complex2DReg.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <Born_full.h>
#include <chrono>
#include <boost/progress.hpp>

using namespace std::chrono;

using namespace SEP;

Born_full::Born_full (std::shared_ptr<float1DReg> wave, std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
std::shared_ptr<float2DReg> scoord, std::shared_ptr<float2DReg> rcoord) :
Acquisition(slow->getHyper(),par,scoord,rcoord) {

	param = par;
	_slow = slow;

	nref = par->getInt("nref",1);
	setBgSlow(slow);

	wave_f = std::make_shared<complex1DReg> (wave->getHyper());

	FFT1 fft1 (wave,wave_f,1,1);
	fft1.forward(wave,wave_f,0);

	construct();

}

// Born_full::Born_full (std::shared_ptr<float3DReg> wave, std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
// std::shared_ptr<float2DReg> scoord, std::shared_ptr<float2DReg> rcoord) :
// Acquisition(slow->getHyper(),par,scoord,rcoord) {
//
// 	param = par;
// 	_slow = slow;
//
// 	nref = par->getInt("nref",1);
// 	setBgSlow(slow);
//
// 	wave_f = std::make_shared<complex3DReg> (wave->getHyper());
//
// 	FFT1 fft1 (wave,wave_f,1,1);
// 	fft1.forward(wave,wave_f,0);
//
// 	construct();
//
// }

void Born_full::forward (std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));
	if(!add) data_f->scale(0.);

	std::cerr << "Forward born ... " << std::endl;
		auto start = high_resolution_clock::now();

	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is<rs.end(); ++is) {

		tbb::parallel_for(tbb::blocked_range<int> (fmin,fmax),
			[=] (const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(_src->getZ(is),_src->getX(is),param->getInt("ng",0)));
			std::shared_ptr<Injection> inj_rec (new Injection(_rec->getZ(),_rec->getX()));
			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()));
			std::shared_ptr<Down> down (new Down(_slow,param,ref));
			std::shared_ptr<Up> up (new Up (_slow,param,ref));

			std::unique_ptr<Born_refl> born_refl(new Born_refl(_slow,inj_rec,up,_bg_wfld,_sc_wfld));
			std::unique_ptr<Born_down> born_down(new Born_down (_slow,param,inj_rec,down,up,reflect,_bg_wfld,_sc_wfld));
			std::unique_ptr<Born_up> born_up (new Born_up(_slow,param,inj_rec,up,_bg_wfld,_sc_wfld));

			inj_rec->setShot(is);
			inj_src->setShot(is);

			for (int i=r.begin(); i<r.end(); i++) {

				// setting up
				inj_src->setStep(i);
				inj_rec->setStep(i);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				born_down->setFreq(freq[i]);
				born_up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(wave_f,_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// linearizetion of reflectivity
				born_refl->forward(model,data_f,1);

				// down-scattering
				born_down->forward(model,data_f,1);
				#ifdef SNAPSHOT
				snap_wfld(_sc_wfld,i);
				#endif

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				// up-scattering
				born_up->forward(model,data_f,1);
				// #ifdef SNAPSHOT
				// snap_wfld(_sc_wfld,i);
				// #endif
			}
		});
	}

	});

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

	FFT1 fft_data(data,data_f,1,1);
	fft_data.adjoint(data,data_f,0);
}

void Born_full::adj_illum (std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) model->scale(0.);
	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));
	FFT1 fft_data(data,data_f,1,1);
	fft_data.forward(data,data_f,0);

	std::cerr << "Adjoint born illumination ... " << std::endl;
		auto start = high_resolution_clock::now();

		std::shared_ptr<Illumination> illum (new Illumination (_slow->getHyper()));
		std::shared_ptr<float2DReg> temp = model->clone();

	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is<rs.end(); ++is) {

		// temp->scale(0);

		tbb::parallel_for(tbb::blocked_range<int> (fmin,fmax),
			[=] (const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(_src->getZ(is),_src->getX(is),param->getInt("ng",0)));
			std::shared_ptr<Injection> inj_rec (new Injection(_rec->getZ(),_rec->getX()));
			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()));
			std::shared_ptr<Down> down (new Down(_slow,param,ref));
			std::shared_ptr<Up> up (new Up (_slow,param,ref));

			std::unique_ptr<Born_refl> born_refl(new Born_refl(_slow,inj_rec,up,_bg_wfld,_sc_wfld));
			std::unique_ptr<Born_down> born_down(new Born_down (_slow,param,inj_rec,down,up,reflect,_bg_wfld,_sc_wfld));
			std::unique_ptr<Born_up> born_up (new Born_up(_slow,param,inj_rec,up,_bg_wfld,_sc_wfld));

			inj_rec->setShot(is);
			inj_src->setShot(is);

			for (int i=r.begin(); i<r.end(); i++) {

				// setting up
				inj_src->setStep(i);
				inj_rec->setStep(i);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				born_down->setFreq(freq[i]);
				born_up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(wave_f,_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				illum->accumulate(_bg_wfld);

				// linearizetion of reflectivity
				born_refl->adjoint(temp,data_f,1);
				#ifdef SNAPSHOT
				snap_wfld(_sc_wfld,i);
				#endif
				// down-scattering
				born_down->adjoint(temp,data_f,1);

				#ifdef SNAPSHOT
				snap_wfld(_sc_wfld,i);
				#endif

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				illum->accumulate(_bg_wfld);
				// up-scattering
				born_up->adjoint(temp,data_f,1);
				#ifdef SNAPSHOT
				snap_wfld(_sc_wfld,i);
				#endif
			}
		});
	}
	});

	illum->compensate(temp);
	model -> scaleAdd(temp,1.,1.);

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

		std::cerr << duration.count()/1e6 << " s" << std::endl;
}

void Born_full::adj_reg (std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) model->scale(0.);

	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));
	FFT1 fft_data(data,data_f,1,1);
	fft_data.forward(data,data_f,0);

	std::cerr << "Adjoint born ... " << std::endl;
		auto start = high_resolution_clock::now();
	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is<rs.end(); ++is) {

		tbb::parallel_for(tbb::blocked_range<int> (fmin,fmax),
			[=] (const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(_src->getZ(is),_src->getX(is),param->getInt("ng",0)));
			std::shared_ptr<Injection> inj_rec (new Injection(_rec->getZ(),_rec->getX()));
			std::shared_ptr<complex2DReg> _bg_wfld(new complex2DReg(_slow->getHyper()));
			std::shared_ptr<complex2DReg> _sc_wfld(new complex2DReg(_slow->getHyper()));
			std::shared_ptr<Down> down (new Down(_slow,param,ref));
			std::shared_ptr<Up> up (new Up (_slow,param,ref));

			std::unique_ptr<Born_refl> born_refl(new Born_refl(_slow,inj_rec,up,_bg_wfld,_sc_wfld));
			std::unique_ptr<Born_down> born_down(new Born_down (_slow,param,inj_rec,down,up,reflect,_bg_wfld,_sc_wfld));
			std::unique_ptr<Born_up> born_up (new Born_up(_slow,param,inj_rec,up,_bg_wfld,_sc_wfld));

			inj_rec->setShot(is);
			inj_src->setShot(is);

			for (int i=r.begin(); i<r.end(); i++) {

				// setting up
				inj_src->setStep(i);
				inj_rec->setStep(i);
				down->setFreq(freq[i]);
				up->setFreq(freq[i]);

				born_down->setFreq(freq[i]);
				born_up->setFreq(freq[i]);

				// down background wavefield
				inj_src->forward(wave_f,_bg_wfld,0);
				down->forward(_bg_wfld,_bg_wfld,1);

				// linearizetion of reflectivity
				born_refl->adjoint(model,data_f,1);
				#ifdef SNAPSHOT
				snap_wfld(_sc_wfld,i);
				#endif
				// down-scattering
				born_down->adjoint(model,data_f,1);

				#ifdef SNAPSHOT
				snap_wfld(_sc_wfld,i);
				#endif

				// up background wavefield
				reflect->forward(_bg_wfld,_bg_wfld,1);
				up->forward(_bg_wfld,_bg_wfld,1);

				// up-scattering
				born_up->adjoint(model,data_f,1);
				#ifdef SNAPSHOT
				snap_wfld(_sc_wfld,i);
				#endif
			}
		});
	}
	});

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

		std::cerr << duration.count()/1e6 << " s" << std::endl;
}
