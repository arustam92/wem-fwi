#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include "WEM.h"
#include <chrono>
#include <boost/progress.hpp>

#ifdef VTUNE
#include <thread>
#include <ittnotify.h>
#endif

using namespace std::chrono;

using namespace SEP;


WEM::WEM (std::shared_ptr<float3DReg> slow, std::shared_ptr<paramObj> par,
std::shared_ptr<float2DReg> scoord, std::shared_ptr<float2DReg> rcoord) :
Acquisition(slow->getHyper(),par,scoord,rcoord) {

	param = par;
	_slow = slow;

	nref = par->getInt("nref",1);
	ref = std::make_shared<RefSampler>(slow,nref);
	reflect = std::make_shared<Reflect>(slow);

	fMIN = par->getFloat("fmin",1);
	fMAX = par->getFloat("fmax",0);

	nshots = _src->getZ().size();
}

void WEM::forward(std::shared_ptr<float1DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) data->scale(0.);

	std::shared_ptr<complex1DReg> model_f (new complex1DReg(model->getHyper()));
	FFT1 fft1 (model,model_f,1,1);
	fft1.forward(model,model_f,0);
	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));
	// data_f->scale(0.);

	float df = 1./(model_f->getHyper()->getAxis(1).d*model_f->getHyper()->getAxis(1).n);
	if (fMAX == 0) {
		fmin = 1;
		fmax = param->getInt("nfreq",0);
		if (fmax==0) fmax = _wave_f->getHyper()->getAxis(1).n/2+1;
	}
	else {
		// calculate indexes
		fmin = fMIN/df;
		fmax = fMAX/df;
	}

	freq.resize(boost::extents[fmax]);
	for (int i=0; i<freq.size(); i++) {
		freq[i] = i*df;
	}

	// loop over shots
	std::cerr << "Propagating shots ... " << std::endl;
			auto start = high_resolution_clock::now();
	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {
	for (int is=rs.begin(); is < rs.end(); ++is) {

		// boost::progress_display show_progress(nfreq,std::cerr);

		tbb::parallel_for(tbb::blocked_range<int>(fmin,fmax),
			[=](const tbb::blocked_range<int> &r) {

			// for vtune
			#ifdef VTUNE
			std::string msg("Thread ");
			// // msg.append(std::string(std::this_thread::get_id()));
			__itt_event event = __itt_event_create(msg.c_str(),msg.size());
			__itt_event_start(event);
			#endif

			std::unique_ptr<Injection> inj_src (new Injection(_src->getZ(is),
												_src->getX(is),param->getInt("ng",0)));
			std::unique_ptr<Injection> inj_rec (new Injection(_rec->getZ(),_rec->getX()));
			std::shared_ptr<complex2DReg> _wfld(new complex2DReg(_slow->getHyper()));
			std::unique_ptr<Down> down (new Down(_slow,param,ref));
			std::unique_ptr<Up> up (new Up (_slow,param,ref));

			inj_rec->setShot(is);

			for (int i=r.begin(); i < r.end(); i++) {

				inj_src->setStep(i);
				inj_src->forward(model_f,_wfld,0);

				down->setFreq(freq[i]);
				down->forward(_wfld,_wfld,1);

				#ifdef SNAPSHOT
				snap_wfld(_wfld,i);
				#endif

				reflect->forward(_wfld,_wfld,1);

				up->setFreq(freq[i]);
				up->forward(_wfld,_wfld,1);

				#ifdef SNAPSHOT
				snap_wfld(_wfld,i);
				#endif

				inj_rec->setStep(i);
				inj_rec->adjoint(data_f,_wfld,1);

				// ++show_progress;

			}
			#ifdef VTUNE
		__itt_event_end(event);
			#endif

		});
		}
	});

	auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

		std::cerr << duration.count()/1e6 << " s" << std::endl;

 	FFT1 fft_data(data,data_f,1,1);
	fft_data.adjoint(data,data_f,0);

}

void WEM::adjoint(std::shared_ptr<float1DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) model->scale(0.);

	std::shared_ptr<complex1DReg> model_f (new complex1DReg(model->getHyper()));
	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));

	FFT1 fft_data(data,data_f,1,1);
	fft_data.forward(data,data_f,0);

	float df = 1./(model_f->getHyper()->getAxis(1).d*model_f->getHyper()->getAxis(1).n);
	if (fMAX == 0) {
		fmin = 1;
		fmax = param->getInt("nfreq",0);
		if (fmax==0) fmax = _wave_f->getHyper()->getAxis(1).n/2+1;
	}
	else {
		// calculate indexes
		fmin = fMIN/df;
		fmax = fMAX/df;
	}

	freq.resize(boost::extents[fmax]);
	for (int i=0; i<freq.size(); i++) {
		freq[i] = (i*df);
	}

	std::cerr << "Propagating adjoint shots ... " << std::endl;
	auto start = high_resolution_clock::now();
	tbb::parallel_for(tbb::blocked_range<int>(0,nshots),
		[=](const tbb::blocked_range<int> &rs) {

	for (int is=rs.begin(); is<rs.end(); ++is) {

		// std::cerr << is << std::endl;

		tbb::parallel_for(tbb::blocked_range<int>(fmin,fmax),
			[=](const tbb::blocked_range<int> &r) {

			std::unique_ptr<Injection> inj_src (new Injection(_src->getZ(is),_src->getX(is),param->getInt("ng",0)));
			std::unique_ptr<Injection> inj_rec (new Injection(_rec->getZ(),_rec->getX()));
			std::shared_ptr<complex2DReg> _wfld(new complex2DReg(_slow->getHyper()));
			std::unique_ptr<Down> down (new Down(_slow,param,ref));
			std::unique_ptr<Up> up (new Up (_slow,param,ref));

			inj_rec->setShot(is);

			for (int i=r.begin(); i<r.end(); i++) {

				inj_rec->setStep(i);
				inj_rec->forward(data_f,_wfld,0);

				up->setFreq(freq[i]);
				up->adjoint(_wfld,_wfld,1);

				#ifdef SNAPSHOT
				snap_wfld(_wfld,i);
				#endif

				reflect->adjoint(_wfld,_wfld,1);

				down->setFreq(freq[i]);
				down->adjoint(_wfld,_wfld,1);

				inj_src->setStep(i);
				inj_src->adjoint(model_f,_wfld,1);

				#ifdef SNAPSHOT
				snap_wfld(_wfld,i);
				#endif
			}
		});

	}

	});

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

		std::cerr << duration.count()/1e6 << " s" << std::endl;


	FFT1 fft1 (model,model_f,1,1);
	fft1.adjoint(model,model_f,0);

}
