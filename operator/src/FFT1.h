#pragma once

#include <tbb/mutex.h>
#include <complexHyper.h>
#include <complex1DReg.h>
#include <float1DReg.h>
#include <floatHyper.h>
#include <Operator.h>
#include <fftw3.h>
#include <axis.h>


namespace SEP{

static tbb::mutex mutex;

class fft_data
{
public:
	fft_data(int size) {
		in = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*size);
		// out = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*size);
	};
	fftwf_complex* getModel() {return in;}
	fftwf_complex* getData() {return out;}
	inline void setData(std::complex<float> *ptr) {
		out = reinterpret_cast<fftwf_complex*>(ptr);
	}
	~fft_data() {
		fftwf_free(in);
		// fftwf_free(out);
	};

	fftwf_complex* in;
	fftwf_complex* out;

private:

};

class FFT1 :
public Operator<complexHyper,complexHyper> {

public:
	FFT1(std::shared_ptr<complexHyper> model, std::string mode, int rank, int a);
	void forward(std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add);
	void adjoint(std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add);

	~FFT1() {
		fftwf_destroy_plan(_fwd_plan);
		fftwf_destroy_plan(_adj_plan);
		// fftwf_cleanup();
	}

	FFT1(std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data, int rank, int a);
	void forward (std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data, bool add);
	void adjoint (std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data, bool add);

private:
	int _axis, _rank;
	int idist, odist, istride, ostride, howmany, size;
	axis new_axis, old_axis;

	fftwf_plan _fwd_plan, _adj_plan;
	fftwf_complex *_in, *_out;
	std::shared_ptr<fft_data> dataFFT;
	complex1D shift;
	float1D k;

};



// class PlannerFFT
// {
// public:
// 	PlannerFFT(std::shared_ptr<floatHyper> model, int rank, int a) {

// 	int ndim = model->getHyper()->getNdimG1();
// 	int n1 = model->getHyper()->getAxis(a).n;

// 	int howmany = 1;
// 	for (int i=ndim;i>rank;i--) howmany *= model->getHyper()->getAxis(i).n;

// 	/*
// 	j-th element of k-th transform is
// 	out = in + j*stride + k*idist
// 	*/
//  	int idist = n1;
//  	int odist = n1;
// 	int istride = 1;
// 	int ostride = 1;

// 	_fwd_plan = fftwf_plan_many_dft(_rank, &n1, howmany,
//                              _in, NULL, istride, idist,
//                              _out, NULL, ostride, odist,
//                              FFTW_FORWARD,FFTW_ESTIMATE);
//  	_adj_plan = fftwf_plan_many_dft(_rank, &n1, howmany,
//                              _out, NULL, istride, idist,
//                              _in, NULL, ostride, odist,
//                              FFTW_BACKWARD,FFTW_ESTIMATE);
// 	};
// 	~PlannerFFT() {
// 		fftwf_destroy_plan(_fwd_plan);
// 		fftwf_destroy_plan(_adj_plan);
// 		fftwf_cleanup();
// 	};

// 	fftwf_plan getFwd() {return _fwd_plan;}
// 	fftwf_plan getAdj() {return _adj_plan;}

// private:
// 	fftwf_plan _fwd_plan, _adj_plan;


// };

}
