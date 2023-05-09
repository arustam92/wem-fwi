#pragma once
#include <PropParam.h>
#include <Born_refl.h>
#include <Born_up.h>
#include <Born_down.h>
#include <illumination.h>
namespace SEP {

class BornRefl : public PropParam
{
public:

	BornRefl(std::shared_ptr<float1DReg> wave, std::shared_ptr<complex2DReg> slow, std::shared_ptr<paramObj> par) :
					PropParam(wave,slow->getHyper(),par) {
						setBgSlow(slow);
						nx = slow->getHyper()->getAxis(1).n;
						nz = slow->getHyper()->getAxis(2).n;
					};

	void forward(std::shared_ptr<complex3DReg> model, std::shared_ptr<float3DReg> data, bool add);
	void adjoint(std::shared_ptr<complex3DReg> model, std::shared_ptr<float3DReg> data, bool add);

	// BornRefl(std::shared_ptr<complex1DReg> wave, std::shared_ptr<complex2DReg> slow, std::shared_ptr<paramObj> par) :
	// 				PropParam(wave,slow->getHyper(),par) {};
	// void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	// void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {
	// 	(*this.*c_ptr_adjoint)(model,data,add);
	// };

	void setBgSlow(std::shared_ptr<complex2DReg> slow) {
		_slow = slow->clone();
		// reflect.reset(new Reflect(slow));
		ref.reset(new RefSampler(slow,nref));
	}



private:

	void (BornRefl::*ptr_adjoint)(std::shared_ptr<complex3DReg>, std::shared_ptr<float3DReg>, bool);
	void adj_reg(std::shared_ptr<complex3DReg> model, std::shared_ptr<float3DReg> data, bool add);
	void adj_illum(std::shared_ptr<complex3DReg> model, std::shared_ptr<float3DReg> data, bool add);

	void (BornRefl::*c_ptr_adjoint)(std::shared_ptr<complex2DReg>, std::shared_ptr<complex3DReg>, bool);
	void c_adj_reg(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	void c_adj_illum(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add);

	void fwd_propagate(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data);
	void adj_propagate(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data);

	void sliceModel(std::shared_ptr<complex3DReg> input,std::shared_ptr<complex2DReg> output, int ifreq);
	void sliceModel(std::shared_ptr<complex2DReg> input,std::shared_ptr<complex3DReg> output, int ifreq);

	int nz, nx;

	std::shared_ptr<complex2DReg> _slow;
	std::shared_ptr<RefSampler> ref;
};

}
