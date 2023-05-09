
// #include <complex2DReg.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <BornFull.h>
#include <chrono>
#include <boost/progress.hpp>

using namespace std::chrono;

using namespace SEP;

BornFull::BornFull (const std::shared_ptr<float1DReg>& wave, const std::vector<std::shared_ptr<complex2DReg>>& model, const std::shared_ptr<paramObj>& par) :
PropParam(wave,model[0]->getHyper(),par) {

	_model.push_back(std::make_shared<complex2DReg>(model[0]->getHyper()));
	_model.push_back(std::make_shared<complex2DReg>(model[0]->getHyper()));
	setBgSlow(model);

	// if (param->getBool("illum",false)) {
	// 	//
	// 	ptr_adjoint = &BornFull::adj_illum;
	// }
	// else {
	// 	ptr_adjoint = &BornFull::adj_reg;
	// }

	// if (param->getBool("onepass",false)) {
	// 	// illum = std::make_shared<Illumination>(_slow->getHyper());
	// 	fwd_propagate = &BornFull::onepass_fwd;
	// 	adj_propagate = &BornFull::onepass_adj;
	// }
	// else {
		fwd_propagate = &BornFull::full_fwd;
		adj_propagate = &BornFull::full_adj;
	// }
}

void BornFull::forward (const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) data->scale(0.);
	int tpad = param->getInt("tpad",0);
	Pad pad_data(data->getHyper(),0,tpad,false);
	auto data_pad = std::make_shared<float3DReg>(pad_data.getPaddedHyper());

	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data_pad->getHyper()));

	std::cerr << "Forward born ... " << std::endl;
		auto start = high_resolution_clock::now();

	(*this.*fwd_propagate)(model,data_f);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

	FFT1 fft_data(data_pad,data_f,1,1);
	fft_data.adjoint(data_pad,data_f,0);
	pad_data.adjoint(data,data_pad,add);

}

void BornFull::forward (const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	std::cerr << "Forward born cmplx... " << std::endl;
		auto start = high_resolution_clock::now();

	(*this.*fwd_propagate)(model,data);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

}

void BornFull::adjoint (std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) std::for_each(model.begin(), model.end(), [](std::shared_ptr<complex2DReg>& v) { v->scale(0); });

	int tpad = param->getInt("tpad",0);
	Pad pad_data(data->getHyper(),0,tpad,false);
	auto data_pad = std::make_shared<float3DReg>(pad_data.getPaddedHyper());
	pad_data.forward(data,data_pad,0);

	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data_pad->getHyper()));
	FFT1 fft_data(data_pad,data_f,1,1);
	fft_data.forward(data_pad,data_f,0);

	std::cerr << "Adjoint born ... " << std::endl;
		auto start = high_resolution_clock::now();

		(*this.*adj_propagate)(model,data_f);

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;
}

void BornFull::adjoint (std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) std::for_each(model.begin(), model.end(), [](std::shared_ptr<complex2DReg>& v) { v->scale(0); });

	std::cerr << "Adjoint born cmplx... " << std::endl;
		auto start = high_resolution_clock::now();

		(*this.*adj_propagate)(model,data);

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;
}
