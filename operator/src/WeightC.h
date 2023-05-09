#pragma once 
#include <float2DReg.h>
#include <complex2DReg.h>
#include <Operator.h>

namespace SEP {

class WeightC : public Operator<complex2DReg,complex2DReg>
{


public:
	WeightC() {}
	WeightC(std::shared_ptr<complex2DReg> weight) {_w = weight;}

	// void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<float2DReg> data, bool add);
	// void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<float2DReg> data, bool add);

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);

	void setWeight(std::shared_ptr<complex2DReg> w) {_w = w;};

private:
	std::shared_ptr<complex2DReg> _w;

};

}