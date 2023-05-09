#pragma once
#include <Operator.h>
#include <paramObj.h>
#include <float2DReg.h>
#include <complex1DReg.h>
#include <Phshift.h>
#include <Taper.h>
#include <Selector.h>
#include <FFT1.h>
#include <RefSampler.h>

namespace SEP{

class OneStep : public Operator<complex1DReg,complex1DReg>

{
public:
	OneStep (std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par, std::shared_ptr<RefSampler> _ref) {

		_dz = slow->getHyper()->getAxis(2).d;

		nref = par->getInt("nref",1);
		ref = _ref;

		_wfld_ref = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
		model_kx = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
		fft_in = std::make_shared<FFT1>(_wfld_ref,"in",1,1);
		fft_out = std::make_shared<FFT1>(_wfld_ref,"out",1,1);
		taper = std::make_shared<Taper>(par->getInt("tap",50));
		ps = std::make_shared<Phshift>(_dz,slow->getHyper());
		select = std::make_shared<Selector>();
	}

	virtual void setFreq(float freq) {};
	virtual void setDepth(int iz) {};

	virtual void forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {};
	virtual void adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add) {};

protected:

	int nref, _iz;
	float _dz;

	std::shared_ptr<complex1DReg> _wfld_ref , model_kx;

	std::shared_ptr<FFT1> fft_in ,fft_out;
	std::shared_ptr<Phshift> ps;
	std::shared_ptr<Taper> taper;
	std::shared_ptr<RefSampler> ref;
	std::shared_ptr<Selector> select;

};
}
