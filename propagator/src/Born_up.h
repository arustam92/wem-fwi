#pragma once

#include <LinOneWay.h>
#include <Injection.h>
#include <OneWay.h>

namespace SEP {

class Born_up : public Operator<complex2DReg,complex3DReg>
{
public:
	Born_up(std::shared_ptr<complex2DReg> slow, std::shared_ptr<paramObj> par,
		std::shared_ptr<Injection> inj_rec, std::shared_ptr<Up> up,
		std::shared_ptr<complex2DReg> bg_wfld, std::shared_ptr<complex2DReg> sc_wfld);

	inline void setFreq(float freq) {lin_up->setFreq(freq);}

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);

private:
	std::shared_ptr<Up> _up;
	std::shared_ptr<Injection> _inj_rec;
	std::shared_ptr<LinUp> lin_up;
	std::shared_ptr<complex2DReg> _sc_wfld;

};

}
