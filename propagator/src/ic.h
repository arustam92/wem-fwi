
#pragma once
#include <Operator.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <functional>

namespace SEP {

class IC : public Operator<complex2DReg,complex2DReg>
{
public:
	IC(const std::shared_ptr<complex2DReg> bg_wfld) : _bg_wfld(bg_wfld) {

		_nx = bg_wfld->getHyper()->getAxis(1).n;
		_N123 = bg_wfld->getHyper()->getN123();
	};

	inline void setDepth(int iz) {_iz = iz;}
	// inline void setBgWfld(std::shared_ptr<complex2DReg> bg_wfld) {_bg_wfld = bg_wfld->clone();};

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);

	void forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);

protected:
	const std::shared_ptr<complex2DReg> _bg_wfld;
	int _N123, _iz, _nx;
};

}

// // simple cross correlation
// void ic(std::shared_ptr<float2DReg> image, std::shared_ptr<complex2DReg> wfld1, std::shared_ptr<complex2DReg> wfld2);
