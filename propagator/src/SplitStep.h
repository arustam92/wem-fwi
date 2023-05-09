#pragma once
#include <complex1DReg.h>
#include <float2DReg.h>
#include <Operator.h>

namespace SEP {

class SplitStep : public Operator<complex1DReg,complex1DReg>
{
public:
	SplitStep(std::shared_ptr<complex2DReg> slow) : _slow(slow) {
		_dz = slow->getHyper()->getAxis(2).d;
	};

	inline void setFreq(float& freq) {w = 2*pi*freq;}
	inline void setLocation(std::vector<int>& loc) {_loc = loc;}
	inline void setSlow(std::complex<float>& s) {_sref = s;}
	inline void setDepth(int iz) {_iz = iz;};

	void forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);

private:
	constexpr static double pi = 4*std::atan(1.);
	float w;
	int _iz;
	float _dz;
	std::complex<float> _sref;
	std::vector<int> _loc;
	std::shared_ptr<complex2DReg> _slow;

};

}
