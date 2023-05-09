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

class extWEM : public PropParam
{
public:

	extWEM(std::shared_ptr<float1DReg> wave, std::shared_ptr<hypercube> modspace, std::shared_ptr<paramObj> par);

	// extWEM(std::shared_ptr<complex1DReg> wave, std::shared_ptr<hypercube> modspace, std::shared_ptr<paramObj> par) :
	// 			PropParam(wave,modspace,par) {;	};

	void forward(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<float3DReg> data, bool add);
	void forward(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data, bool add);

	void wavefield(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<float4DReg> data, bool add);

private:
	void (extWEM::*fwd_propagate)(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data);
	void full_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data);
	void full_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex4DReg> data);
	void onepass_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data);
	void partial_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data);


};

}
