#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include "WEM.h"
#include <chrono>

using namespace std::chrono;

using namespace SEP;


WEM::WEM (std::shared_ptr<float1DReg> wave, std::shared_ptr<float2DReg> modspace,
		 std::shared_ptr<paramObj> par,std::shared_ptr<float2DReg> scoord, std::shared_ptr<float2DReg> rcoord) :
		 Acquisition(modspace->getHyper(),par,scoord,rcoord) {

	param = par;

	_wave_f = std::make_shared<complex1DReg>(wave->getHyper());
	FFT1 fft_wave (wave,_wave_f, 1, 1);
	fft_wave.forward(wave,_wave_f,false);

	nref = par->getInt("nref",1);

	float df = 1./(_wave_f->getHyper()->getAxis(1).d*_wave_f->getHyper()->getAxis(1).n);
	fMIN = par->getFloat("fmin",1);
	fMAX = par->getFloat("fmax",0);
	if (fMAX == 0) {
		fmin = 1;
		fmax = par->getInt("nfreq",0);
		if (fmax==0) fmax = _wave_f->getHyper()->getAxis(1).n/2+1;
	}
	else {
		// calculate indexes
		fmin = fMIN/df;
		fmax = fMAX/df;
	}

	freq.resize(boost::extents[fmax+1]);
	for (int i=0; i<freq.size(); i++) {
		freq[i] = (i*df);
	}
	nshots = _src->getZ().size();

}

WEM::WEM (std::shared_ptr<complex1DReg> wave, std::shared_ptr<float2DReg> modspace,std::shared_ptr<paramObj> par,
		std::shared_ptr<float2DReg> scoord, std::shared_ptr<float2DReg> rcoord) :
		Acquisition(modspace->getHyper(),par,scoord,rcoord) {

	param = par;
	_wave_f = wave;

	nref = par->getInt("nref",1);

	float df = _wave_f->getHyper()->getAxis(1).d;

	fmin = 0;
	fmax = _wave_f->getHyper()->getAxis(1).n;

	freq.resize(boost::extents[fmax]);
	for (int i=0; i<freq.size(); i++) {
		freq[i] = (i*df) + _wave_f->getHyper()->getAxis(1).o;
	}
	nshots = _src->getZ().size();
}

void WEM::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) data->scale(0.);
	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));
	// data_f->scale(0.);

	ref = std::make_shared<RefSampler>(model,nref);
	reflect = std::make_shared<Reflect>(model);

	// loop over shots
	std::cerr << "Propagating shots ... " << std::endl;
			auto start = high_resolution_clock::now();
	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is < rs.end(); ++is) {

		tbb::parallel_for(tbb::blocked_range<int>(fmin,fmax),
			[=](const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(_src->getZ(is),
												_src->getX(is),param->getInt("ng",0)));
			std::unique_ptr<Injection> inj_rec (new Injection(_rec->getZ(),_rec->getX()));
			std::shared_ptr<complex2DReg> _wfld(new complex2DReg(model->getHyper()));
			std::unique_ptr<Down> down (new Down(model,param,ref));
			std::unique_ptr<Up> up (new Up (model,param,ref));

			inj_rec->setShot(is);

			for (int i=r.begin(); i < r.end(); i++) {

				inj_src->setStep(i);
				inj_src->forward(_wave_f,_wfld,0);
				//
				down->setFreq(freq[i]);
				down->forward(_wfld,_wfld,1);
				//
				// // if(wantwfld) inj_src->adjoint(_full_wfld_w,_wfld,1);
				//
				reflect->forward(_wfld,_wfld,1);

				up->setFreq(freq[i]);
				up->forward(_wfld,_wfld,1);

				inj_rec->setStep(i);
				inj_rec->adjoint(data_f,_wfld,1);

				// if(wantwfld) inj_rec[ithread]->adjoint(_full_wfld_w,_wfld[ithread],1);

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

void WEM::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	ref.reset(new RefSampler(model,nref));
	reflect.reset(new Reflect(model));

	// loop over shots
	std::cerr << "Propagating complex shots ... " << std::endl;
			auto start = high_resolution_clock::now();

	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is < rs.end(); ++is) {

			tbb::parallel_for(tbb::blocked_range<int>(fmin,fmax),
			[=](const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(_src->getZ(is),
												_src->getX(is),param->getInt("ng",0)));
			std::unique_ptr<Injection> inj_rec (new Injection(_rec->getZ(),_rec->getX()));
			std::shared_ptr<complex2DReg> _wfld(new complex2DReg(model->getHyper()));
			std::unique_ptr<Down> down (new Down(model,param,ref));
			std::unique_ptr<Up> up (new Up (model,param,ref));

			inj_rec->setShot(is);

			for (int i=r.begin(); i < r.end(); i++) {

				inj_src->setStep(i);
				inj_src->forward(_wave_f,_wfld,0);

				down->setFreq(freq[i]);
				down->forward(_wfld,_wfld,1);

				reflect->forward(_wfld,_wfld,1);

				up->setFreq(freq[i]);
				up->forward(_wfld,_wfld,1);

				inj_rec->setStep(i);
				inj_rec->adjoint(data,_wfld,1);


			}
		});
		}
	});

	auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

		std::cerr << duration.count()/1e6 << " s" << std::endl;

}
