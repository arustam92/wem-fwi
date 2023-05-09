#pragma once

#include <float2DReg.h>
#include <complex2DReg.h>
#include <dReflect.h>
#include <Weight.h>
#include <OneWay.h>
#include <Injection.h>
#include <ic.h>


namespace SEP {


class Born_trans : public Operator<complex2DReg,complex3DReg>
{
public:
	Born_trans(std::shared_ptr<complex2DReg>& slow, std::shared_ptr<Injection>& inj_rec,
	 std::shared_ptr<Down>& down, std::shared_ptr<complex2DReg>& bg_wfld,std::shared_ptr<complex2DReg>& sc_wfld);

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);

	std::shared_ptr<dTransmission> dtrans;

private:
	std::shared_ptr<complex2DReg> dt;
	std::shared_ptr<complex2DReg>& _bg_wfld, _sc_wfld;
	std::shared_ptr<Down>& _down;
	std::shared_ptr<Injection>& _inj_rec;
	std::shared_ptr<IC> ic;

	// void (Born_refl::*ptr_adjoint)(std::shared_ptr<float2DReg>, std::shared_ptr<complex3DReg>, bool);
	// void adj_reg(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	// void adj_mixing(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);

};

}
