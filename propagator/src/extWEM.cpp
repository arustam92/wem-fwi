
#include "extWEM.h"
#include <chrono>

using namespace std::chrono;

using namespace SEP;

extWEM::extWEM (std::shared_ptr<float1DReg> wave, std::shared_ptr<hypercube> modspace,std::shared_ptr<paramObj> par) :
		 PropParam(wave,modspace,par) {

	if (modspace->getAxis(3).n != freq.size())
	throw SEPException(std::string("Model is of the wrong size!"));

	if (param->getBool("onepass",false)) {
		// illum = std::make_shared<Illumination>(_slow->getHyper());
		fwd_propagate = &extWEM::onepass_fwd;
	}
	else if (param->getBool("ext_partial",false)) {
		// illum = std::make_shared<Illumination>(_slow->getHyper());
		fwd_propagate = &extWEM::partial_fwd;
	}
	else {
		fwd_propagate = &extWEM::full_fwd;
	}
}

void extWEM::forward(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) data->scale(0.);
	int tpad = param->getInt("tpad",0);
	Pad pad_data(data->getHyper(),0,tpad,false);
	auto data_pad = std::make_shared<float3DReg>(pad_data.getPaddedHyper());

	std::shared_ptr<complex3DReg> data_f = std::make_shared<complex3DReg>(data_pad->getHyper());

	// loop over shots
	std::cerr << "Propagating shots ... " << std::endl;
	auto start = high_resolution_clock::now();

	(*this.*fwd_propagate)(model,data_f);

	auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

 	FFT1 fft_data(data_pad,data_f,1,1);
	fft_data.adjoint(data_pad,data_f,false);
	pad_data.adjoint(data,data_pad,add);

}

void extWEM::forward(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);


	// loop over shots
	std::cerr << "Propagating complex shots ... " << std::endl;
			auto start = high_resolution_clock::now();

	(*this.*fwd_propagate)(model,data);

	auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

}

void extWEM::wavefield(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<float4DReg> data, bool add) {

	if(!add) data->scale(0.);
	std::shared_ptr<complex4DReg> data_f = std::make_shared<complex4DReg>(data->getHyper());

	// loop over shots
	std::cerr << "Propagating shots ... " << std::endl;
			auto start = high_resolution_clock::now();

	full_fwd(model,data_f);

	auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

 	FFT1 fft_data(data,data_f,1,1);
	fft_data.adjoint(data,data_f,0);

}
