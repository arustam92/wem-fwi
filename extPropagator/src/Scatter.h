#pragma once
#include <complex2DReg.h>
#include <float1DReg.h>
#include <float2DReg.h>
#include <Operator.h>
#include <cmath>
#include <FFT1.h>
#include <hypercube.h>
#include <Phshift.h>


namespace SEP {

class Scatter : public Operator<complex1DReg,complex1DReg>
{

public:

	Scatter(std::shared_ptr<float2DReg> slow, int ntaylor);

	inline void setFreq(float freq){w = 2*freq*pi;}
	inline void setDepth(int iz) {_iz = iz;}

	void scatter(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, int it);
	void slow_scale(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, int iz, int it);
	void set_next_wfld(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, int sign);

	void forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add);

private:

	std::vector<float> coef = {1., 1./2. , 3./8. , 5./16. , 35./128.};
	// std::vector<float> coef = {1./2., 1./4.};

	float _dz;
	int nx, _iz, _ntaylor;
	float w;
	float1D k;

	std::shared_ptr<float2DReg> _slow;
	std::shared_ptr<complex1DReg> _wfld_ref, _wfld_kx;
	std::shared_ptr<FFT1> fft_in, fft_out;


};

}
