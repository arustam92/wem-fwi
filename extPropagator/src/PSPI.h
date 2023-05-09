#pragma once
#include <OneStep.h>

namespace SEP{

class PSPI : public OneStep

{
public:
	PSPI (std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par, std::shared_ptr<RefSampler> _ref) : 
	OneStep(slow,par,_ref) {}

	inline void setDepth(int iz) {_iz = iz;}
	inline void setFreq(float freq) {ps->setFreq(freq);}

	void forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add);

};
}