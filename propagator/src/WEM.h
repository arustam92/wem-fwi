#pragma once
#include <Operator.h>
#include <float1DReg.h>
#include <float2DReg.h>
#include <float3DReg.h>
#include <Reflect.h>
#include <Injection.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <complex3DReg.h>
#include <PropParam.h>
#include <OneWay.h>

namespace SEP {

class WEM : public PropParam
{
public:
	// WEM(std::shared_ptr<complex2DReg> slow, std::shared_ptr<paramObj> par);

	WEM(std::shared_ptr<float1DReg> wave, std::shared_ptr<hypercube> modspace, std::shared_ptr<paramObj> par);
	WEM(std::shared_ptr<complex1DReg> wave, std::shared_ptr<hypercube> modspace, std::shared_ptr<paramObj> par);

	// void forward(std::shared_ptr<float1DReg> model, std::shared_ptr<float3DReg> data, bool add);
	// void adjoint(std::shared_ptr<float1DReg> model, std::shared_ptr<float3DReg> data, bool add);

	void forward(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<float3DReg> data, bool add);
	void forward(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data, bool add);
	// void wavefield(const std::vector<std::shared_ptr<complex2DReg>> model, std::shared_ptr<float4DReg> data, bool add);

private:

	void (WEM::*fwd_propagate)(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data);
	void full_fwd(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data);
	// void full_fwd(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex4DReg> data);
	// void onepass_fwd(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data);

	std::shared_ptr<complex3DReg> _full_wfld_w;
	std::shared_ptr<float3DReg> _full_wfld_t;
	std::shared_ptr<complex2DReg> _slow; // for linear op

};

}
