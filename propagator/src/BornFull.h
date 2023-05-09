#pragma once
#include <PropParam.h>
#include <Born_refl.h>
#include <Born_up.h>
#include <Born_down.h>
#include <illumination.h>

namespace SEP {

class BornFull : public Operator<complex2DReg,complex3DReg>, public PropParam
{
public:
	// BornFull();

	BornFull(const std::shared_ptr<float1DReg>& wave, const std::vector<std::shared_ptr<complex2DReg>>& slow, const std::shared_ptr<paramObj>& par);
	// BornFull(std::shared_ptr<float3DReg> wave, std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
	// 				std::shared_ptr<float2DReg> scoord = nullptr, std::shared_ptr<float2DReg> rcoord = nullptr);
	void forward(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<float3DReg> data, bool add);
	void adjoint(std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<float3DReg> data, bool add);

	void forward(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data, bool add);
	void adjoint(std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data, bool add);

	void setBgSlow(const std::vector<std::shared_ptr<complex2DReg>>& model) {
		_model[0] = model[0]->clone();
		_model[1] = model[1]->clone();
		// ref.reset(new RefSampler(slow,nref));
		// reflect = std::make_shared<Reflect>(slow);
	}

	inline std::shared_ptr<complex2DReg> getBgSlow() {return _model[0];};
	inline std::shared_ptr<complex2DReg> getBgDensity() {return _model[1];};

private:
	void (BornFull::*fwd_propagate)(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data);
	void (BornFull::*adj_propagate)(std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data);

	// void onepass_fwd(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data);
	// void onepass_adj(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data);

	void full_fwd(const std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data);
	void full_adj(std::vector<std::shared_ptr<complex2DReg>>& model, std::shared_ptr<complex3DReg> data);

protected:

	std::vector<std::shared_ptr<complex2DReg>> _model;



};

}
