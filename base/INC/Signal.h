#include "float2DReg.h"
#include "float1DReg.h"
#include "complex2DReg.h"
#include "FFT1.h"
#include <omp.h>

namespace SEP { 
class Signal : public float2DReg {
	
public:
	Signal(std::shared_ptr<hypercube> hyper) : float2DReg(hyper) { }


	std::shared_ptr<Signal> xcorr(const std::shared_ptr <float2DReg> x); 	// compute trace by trace
	std::shared_ptr<Signal> xcorr(const std::shared_ptr <float1DReg> x);	// correlate each trace with 1d signal

	std::shared_ptr<Signal> conv(const std::shared_ptr <float2DReg> x);

	void levinson(const std::shared_ptr <float2DReg> acorr, const std::shared_ptr <float2DReg> xcorr, float eps);

	std::shared_ptr<float2DReg> kolmogorov(float eps, bool wantwave);

	
};
}