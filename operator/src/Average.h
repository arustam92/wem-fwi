#pragma once

#include <float1DReg.h>
#include <float2DReg.h>
#include <Operator.h>

namespace SEP {

class Average : public Operator<float2DReg,float1DReg> {
public:
	Average() { }
	
	void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<float1DReg> data, bool add);
	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<float1DReg> data, bool add);

};

}
