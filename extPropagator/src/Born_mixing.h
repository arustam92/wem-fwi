#pragma once
#include <WEM.h>
#include <Born_full.h>
#include <illumination.h>

namespace SEP {

class Born_mixing : public Born_full
{
public:

	Born_mixing(std::shared_ptr<float1DReg> wave, std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
							std::vector<float> eps) {
		_eps = eps;

		if (par->getBool("illum",false)) {
			illum = std::make_shared<Illumination>(_slow->getHyper());
			ptr_adjoint = &Born_mixing::adj_illum;
		}
		else {
			ptr_adjoint = &Born_mixing::adj_reg;
		}
	};

	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add) {
		iter++;
		(*this.*ptr_adjoint)(model,data,add);
	};

private:

	std::vector<float> _eps;

	void (Born_mixing::*ptr_adjoint)(std::shared_ptr<float2DReg>, std::shared_ptr<float3DReg>, bool);
	void adj_reg(std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add);
	void adj_illum(std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add);
};

}
