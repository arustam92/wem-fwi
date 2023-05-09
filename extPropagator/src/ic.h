
#pragma once
#include <Operator.h>
#include <float2DReg.h>
#include <complex2DReg.h>
#include <functional>

namespace SEP {

class IC : public Operator<float2DReg,complex2DReg>
{
public:
	IC(std::shared_ptr<complex2DReg> bg_wfld): _bg_wfld(bg_wfld) {};

	void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add);

protected:
	std::shared_ptr<complex2DReg> _bg_wfld;
};

}

// // simple cross correlation
// void ic(std::shared_ptr<float2DReg> image, std::shared_ptr<complex2DReg> wfld1, std::shared_ptr<complex2DReg> wfld2);
