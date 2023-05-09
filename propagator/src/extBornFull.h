#pragma once
#include "WEM.h"
#include "BornFull.h"

namespace SEP {

class extBornFull : public Operator<complex3DReg,complex3DReg>, PropParam
{
public:
	// BornFull();

	extBornFull(std::shared_ptr<float1DReg> wave, std::vector<std::shared_ptr<complex3DReg>> model, std::shared_ptr<paramObj> par);
	void forward(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<float3DReg> data, bool add);
	void adjoint(std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<float3DReg> data, bool add);

	void forward(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data, bool add);
	void adjoint(std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data, bool add);

	// BornFull(std::shared_ptr<complex1DReg> wave, std::shared_ptr<complex3DReg> slow, std::shared_ptr<paramObj> par, bool verb = true,
	// 				std::shared_ptr<float2DReg> scoord = nullptr, std::shared_ptr<float2DReg> rcoord = nullptr);
	// void forward(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	// void adjoint(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, bool add) {
	// 	(*this.*c_ptr_adjoint)(model,data,add);
	// };

	void setBgSlow(const std::vector<std::shared_ptr<complex3DReg>>& bg) {
		pad_slow->forward(bg[0],_bg[0],0);
		pad_slow->forward(bg[1],_bg[1],0);
		// stack->forward(_slow,_slow2d,0);
		// ref.reset(new RefSampler(slow,nref));
		// reflect.reset(new Reflect(slow));
	}
	inline std::shared_ptr<complex3DReg> getBgSlow() {return _bg[0];};
	inline std::shared_ptr<complex3DReg> getBgDensity() {return _bg[1];};

private:
	void (extBornFull::*fwd_propagate)(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data);
	void (extBornFull::*adj_propagate)(std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data);

	// void onepass_fwd(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data);
	// void onepass_adj(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data);

	void full_fwd(const std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data);
	void full_adj(std::vector<std::shared_ptr<complex3DReg>>& model, std::shared_ptr<complex3DReg> data);

	// void partial_fwd(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data);
	// void partial_adj(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data);

	// void dru_fwd(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data);
	// void dru_adj(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data);

	// void (extBornFull::*ptr_adjoint)(std::shared_ptr<complex3DReg>, std::shared_ptr<float3DReg>, bool);
	// void adj_reg(std::shared_ptr<complex3DReg> model, std::shared_ptr<float3DReg> data, bool add);
	// void adj_illum(std::shared_ptr<complex3DReg> model, std::shared_ptr<float3DReg> data, bool add);

	// void (extBornFull::*c_ptr_adjoint)(std::shared_ptr<complex3DReg>, std::shared_ptr<complex3DReg>, bool);
	// void c_adj_reg(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	// void c_adj_illum(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, bool add);

protected:

	std::vector<std::shared_ptr<complex3DReg>> _bg;
	std::shared_ptr<complex2DReg> _slow2d;

	// typedef complex4D::array_view<3>::type c4Dview;
	// typedef boost::multi_array_types::index_range range;
	// complex4D::index_gen indices;

	// tbb::cache_aligned_allocator<Injection> inj_alloc;
	// tbb::cache_aligned_allocator<complex2DReg> c2d_alloc;
	// tbb::cache_aligned_allocator<complex3DReg> c3d_alloc;
	// tbb::cache_aligned_allocator<Down> ow_alloc;
	// tbb::cache_aligned_allocator<RefSampler> ref_alloc;
	// tbb::cache_aligned_allocator<Reflect> refl_alloc;

};

}
