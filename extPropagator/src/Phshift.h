#pragma once

#include <Operator.h>
#include <complex1DReg.h>
#include <float1DReg.h>
#include <cmath>
#include <hypercube.h>

namespace SEP {

class Scatter;

const double pi = 4*std::atan(1.);

class Phshift : public Operator<complex1DReg,complex1DReg>
{
public:
	Phshift (float dz, std::shared_ptr<hypercube> hyp);

	void forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);

	inline void setFreq(float freq) {w = 2*pi*freq;}
	inline void setSlow(float slow1d) {
		_s = slow1d;
	}

	inline void lookupTable();

protected:
	float w;
	float _dz;
	float1D k;
	std::shared_ptr<complex1D> _sqrt;
	std::shared_ptr<complex1D> kz_prev, kz;
	std::shared_ptr<float1D> _amp;
	int ntable;
	float dtable, dk;

	float min,max,std_dev;

	float _s;
};


}
