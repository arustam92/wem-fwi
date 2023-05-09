
#include "WEM.h"
#include "Hilbert.h"
#include <chrono>

using namespace std::chrono;

using namespace SEP;


WEM::WEM (std::shared_ptr<float1DReg> wave, std::shared_ptr<hypercube> modspace,std::shared_ptr<paramObj> par) :
		 PropParam(wave,modspace,par) {

	// if (param->getBool("onepass",false)) {
	// 	// illum = std::make_shared<Illumination>(_slow->getHyper());
	// 	fwd_propagate = &WEM::onepass_fwd;
	// }
	// else {
		fwd_propagate = &WEM::full_fwd;
	// }
}

WEM::WEM (std::shared_ptr<complex1DReg> wave, std::shared_ptr<hypercube> modspace,std::shared_ptr<paramObj> par) :
		PropParam(wave,modspace,par) {

	_wave_f = wave;

	float df = _wave_f->getHyper()->getAxis(1).d;

	int fmin = 0;
	int fmax = _wave_f->getHyper()->getAxis(1).n;

	freq.resize(boost::extents[fmax]);
	for (int i=0; i<freq.size(); i++) {
		freq[i] = (i*df) + _wave_f->getHyper()->getAxis(1).o;
	}
}

void WEM::forward(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) data->scale(0.);
	int tpad = param->getInt("tpad",0);
	Pad pad_data(data->getHyper(),0,tpad,false);
	auto data_pad = std::make_shared<float3DReg>(pad_data.getPaddedHyper());

	std::shared_ptr<complex3DReg> data_f = std::make_shared<complex3DReg>(data_pad->getHyper());

	std::cerr << "Propagating shots ... " << std::endl;
			auto start = high_resolution_clock::now();

	(*this.*fwd_propagate)(model,data_f);

	auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

	// Hilbert hilbert(tpad);
	// auto data_h = data_f->clone();
	// hilbert.forward(data_f, data_h, 0);
 	FFT1 fft_data(data_pad,data_f,1,1);
	fft_data.adjoint(data_pad,data_f,0);
	pad_data.adjoint(data,data_pad,add);

}

void WEM::forward(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	// loop over shots
	std::cerr << "Propagating complex shots ... " << std::endl;
			auto start = high_resolution_clock::now();

	(*this.*fwd_propagate)(model,data);

	auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

}

// void WEM::wavefield(std::vector<std::shared_ptr<complex2DReg>> model, std::shared_ptr<float4DReg> data, bool add) {

// 	if(!add) data->scale(0.);
// 	std::shared_ptr<complex4DReg> data_f = std::make_shared<complex4DReg>(data->getHyper());
// 	// data_f->scale(0.);

// 	// ref = std::make_shared<RefSampler>(model,nref);
// 	// reflect = std::make_shared<Reflect>(model);

// 	// loop over shots
// 	std::cerr << "Recording wavefield ... " << std::endl;
// 			auto start = high_resolution_clock::now();

// 	full_fwd(model,data_f);

// 	auto stop = high_resolution_clock::now();
// 		auto duration = duration_cast<microseconds>(stop - start);

// 	std::cerr << duration.count()/1e6 << " s" << std::endl;

//  	FFT1 fft_data(data,data_f,1,1);
// 	fft_data.adjoint(data,data_f,0);

// }
