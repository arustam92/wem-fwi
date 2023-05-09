#pragma once
#include <OneStep.h>
#include <SplitStep.h>

namespace SEP{

class SSF : public OneStep

{
public:
	SSF (std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par, std::shared_ptr<RefSampler> _ref) : 
	OneStep(slow,par,_ref) {
		ss = std::make_shared<SplitStep>(slow);
	}

	inline void setDepth(int iz) {_iz = iz;ss->setDepth(iz);}
	inline void setFreq(float freq) {ps->setFreq(freq);ss->setFreq(freq);}

	void forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add);

private:
	std::shared_ptr<SplitStep> ss;

};
}