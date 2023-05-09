#pragma once
#include <OneStep.h>
#include <SplitStep.h>

namespace SEP{

class StatPhase : public OneStep

{
public:
	StatPhase (std::shared_ptr<complex2DReg> slow, std::shared_ptr<paramObj> par, std::shared_ptr<RefSampler> _ref) :
	OneStep(slow,par,_ref) {
		_slow_ = slow;
		computeLookup();

		nx = _slow_->getHyper()->getAxis(1).n;
		dk = 2*pi/(nx*_slow_->getHyper()->getAxis(1).d);
		k.resize(boost::extents[nx]);
		for (int ik=1 - k.size()%2; ik<k.size()/2; ik++) {
			k[ik] = ik*dk;
			k[nx-ik - k.size()%2] = -k[ik];
		}
		k[nx/2] = (nx/2)*dk;

		_slow_slice.resize(boost::extents[nx]);
	}

	inline void setDepth(int iz) {_iz = iz; slice_model();}
	inline void setFreq(float freq) {w = 2*pi*freq;}

	void forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add);

private:
	constexpr static float pi = 4*std::atan(1.f);
	void computeLookup() {
		axis ax1 = _slow_->getHyper()->getAxis(1);
		scale1.resize(boost::extents[ax1.n]);
		scale21.resize(boost::extents[ax1.n]);
		scale22.resize(boost::extents[ax1.n]);
		scale3.resize(boost::extents[ax1.n]);

		for (int i=0; i < ax1.n; ++i) {
			float x = ax1.o + i*ax1.d;
			float t = x / std::sqrt(x*x + _dz*_dz);
			float tt = 1 - t*t;
			scale1[i] = t;
			scale21[i] = _dz * std::sqrt(tt) - t*x;
			scale22[i] = _dz * std::sqrt(tt) + t*x;
			scale3[i] = -2*pi*std::sqrt(tt*tt*tt) / _dz;
		}

	};

	void slice_model() {
		std::copy_n(_slow_->_mat->data()+_iz*nx,nx,_slow_slice.data());
	};

	float1D k;
	float1D scale1, scale21,scale22, scale3;
	std::shared_ptr<complex2DReg> _slow_;
	complex1D _slow_slice;
	float k1, k2, k0;	//stationary points
	int nx;
	float w, dk;
};
}
