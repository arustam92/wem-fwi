
// #include <complex2DReg.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <BornTomo.h>
#include <chrono>
#include <boost/progress.hpp>

using namespace std::chrono;

using namespace SEP;

void BornTomo::forward (std::shared_ptr<complex2DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) data->scale(0.);
	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));

	std::cerr << "Forward BornTomo ... " << std::endl;
		auto start = high_resolution_clock::now();

	fwd_propagate(model,data_f);

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;

	FFT1 fft_data(data,data_f,1,1);
	fft_data.adjoint(data,data_f,add);
}

void BornTomo::adjoint (std::shared_ptr<complex2DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	if(!add) model->scale(0.);
	std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));
	FFT1 fft_data(data,data_f,1,1);
	fft_data.forward(data,data_f,0);

	std::cerr << "Adjoint BornTomo illumination ... " << std::endl;
		auto start = high_resolution_clock::now();

	adj_propagate(model, data_f);

	// illum->compensate(model);

		auto stop = high_resolution_clock::now();
		auto duration = duration_cast<microseconds>(stop - start);

	std::cerr << duration.count()/1e6 << " s" << std::endl;
}

void BornTomo::adj_reg (std::shared_ptr<complex2DReg> model, std::shared_ptr<float3DReg> data, bool add) {

	// if(!add) model->scale(0.);
	//
	// std::shared_ptr<complex3DReg> data_f (new complex3DReg(data->getHyper()));
	// FFT1 fft_data(data,data_f,1,1);
	// fft_data.forward(data,data_f,0);
	//
	// std::cerr << "Adjoint born ... " << std::endl;
	// 	auto start = high_resolution_clock::now();
	//
	// 	(*this.*adj_propagate)(model,data_f);
	//
	// 	auto stop = high_resolution_clock::now();
	// 	auto duration = duration_cast<microseconds>(stop - start);
	//
	// std::cerr << duration.count()/1e6 << " s" << std::endl;
}
