#pragma once

// #include <Born_base.h>
#include <Scatter.h>
#include <Weight.h>
#include <LinOneWay.h>
#include <OneWay.h>
#include <Injection.h>
#include <Reflect.h>

namespace SEP {

class Born_down : public Operator<complex2DReg,complex3DReg>
{
public:
	Born_down(const std::shared_ptr<complex2DReg>& slow, std::shared_ptr<paramObj>& par,
			std::shared_ptr<Injection>& inj_rec,
			std::shared_ptr<Down>& down, std::shared_ptr<Up>& up, std::shared_ptr<Operator<complex2DReg,complex2DReg>>& reflect,
			std::shared_ptr<complex2DReg>& bg_wfld,std::shared_ptr<complex2DReg>& sc_wfld);

	inline void setFreq(float freq) {lin_down->setFreq(freq);}

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);

private:
	std::shared_ptr<Down>& _down;
	std::shared_ptr<Up>& _up;
	std::shared_ptr<Operator<complex2DReg,complex2DReg>>& _reflect;
	std::shared_ptr<Injection>& _inj_rec;
	std::shared_ptr<LinDown> lin_down;

	std::shared_ptr<complex2DReg>& _sc_wfld;

};

}
