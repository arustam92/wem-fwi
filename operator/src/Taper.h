#pragma once
#include <Operator.h>
#include "float1DReg.h"
#include <complex1DReg.h>

namespace SEP {

class Taper : public Operator<complex1DReg,complex1DReg>
{
public:
	Taper(int n, int tap) {

		_tap = tap;
		filter.resize(boost::extents[n]);
		std::fill(filter.data(),filter.data()+n,1);
		double pi = 4. * std::atan(1.);
		// fill in values
		for (int i = 0; i < tap; i++) {
			filter[i] = .5f * (1 - std::cos(pi*i/(tap-1)));
			filter[n - 1 - i] = filter[i];
		}

	};

	void forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data,bool add);
	void adjoint(std::shared_ptr<complex1DReg> model,std::shared_ptr<complex1DReg> data,bool add);

private:
	int _tap;
	float1D filter;
};

}
