#pragma once

#include <float2DReg.h>
#include <complex2DReg.h>
#include <dReflect.h>
#include <Weight.h>
#include <OneWay.h>
#include <Injection.h>
#include <ic.h>


namespace SEP {


class Born_refl : public Operator<float2DReg,complex3DReg>
{
public:
	Born_refl(std::shared_ptr<float2DReg> slow, std::shared_ptr<Injection> inj_rec,
	 std::shared_ptr<Up> up, std::shared_ptr<complex2DReg> bg_wfld,std::shared_ptr<complex2DReg> sc_wfld);

	void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);

	void setWfld(std::shared_ptr<complex2DReg> wfld) {_sc_wfld = wfld;};

	std::shared_ptr<dReflect> drefl;

private:
	std::shared_ptr<Weight> w;
	std::shared_ptr<float2DReg> dr;
	std::shared_ptr<complex2DReg> _bg_wfld, _sc_wfld;
	std::shared_ptr<Up> _up;
	std::shared_ptr<Injection> _inj_rec;
	std::shared_ptr<IC> ic;

	// void (Born_refl::*ptr_adjoint)(std::shared_ptr<float2DReg>, std::shared_ptr<complex3DReg>, bool);
	// void adj_reg(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	// void adj_mixing(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);

};

}
