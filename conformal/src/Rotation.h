#pragma once

#include "float1DReg.h"
#include "complex1DReg.h"
#include <complex2DReg.h>
#include "ConformalMap.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace SEP {

class Rotation : public ConformalMap {
public:
	// f(z) = a(z - b)
	Rotation(std::shared_ptr<hypercube> inHyper, float rad) : ConformalMap(inHyper) {
		v.resize(4); u.resize(4);
		std::complex<float> I = {0,1};
		_rad_ = rad;
		_rot_ = std::exp(-I*rad);
		_irot_ = std::exp(I*rad);
		allocateOutSpace();
		complexify();
		mapCoord();
		mapJacob();
	}

	~Rotation() {};

private:

	void allocateOutSpace() {
		axis ax1 = getInHyper()->getAxis(1);
		axis ax2 = getInHyper()->getAxis(2);
		// vertices
		v[0] = {ax2.o, ax1.o};
		v[1] = {ax2.o, ax1.o + (ax1.n-1)*ax1.d};
		v[2] = {ax2.o + (ax2.n-1)*ax2.d, ax1.o};
		v[3] = {ax2.o + (ax2.n-1)*ax2.d, ax1.o + (ax1.n-1)*ax1.d};
		for (int i=0; i < 4; ++i) {
			u[i] = v[i]*_rot_;
			o1 = std::min(o1, u[i].imag());
			o2 = std::min(o2, u[i].real());
			max1 = std::max(max1, u[i].imag());
			max2 = std::max(max2, u[i].real());
		}
		_outSpace_ = std::make_shared<complex2DReg>(axis(ax1.n, o1, (max1-o1)/(ax1.n-1)), axis(ax2.n, o2, (max2-o2)/(ax2.n-1)));
	}

	// construct a map between coordinates of the 2 domains as z = inv(f)(u)
	void mapCoord(){
		std::complex<float> z;
		int i1, i2;
		axis ax1 = getInHyper()->getAxis(1);
		axis ax2 = getInHyper()->getAxis(2);
		float m1 = (ax1.n - 1)*ax1.d + ax1.o;
		float m2 = (ax2.n - 1)*ax2.d + ax2.o;
		float x;
		for (int i=0; i < zcoord.size(); ++i){
			z = getOutCoord()[i] * _irot_;
			i1 = (z.imag() - ax1.o) / ax1.d;
			i2 = (z.real() - ax2.o) / ax2.d;
			if (i1 < 0) {
				x = -(z.imag()-ax1.o) / std::tan(_rad_);
				x = std::min(z.real() + x, m2);
				z = {x, ax1.o};
			}
			if (i2 < 0) {
				x = -(z.real()-ax2.o) * std::tan(_rad_);
				x = std::min(z.imag() + x, m1);
				z = {ax2.o, x};
			}
			if (i1 >= ax1.n) {
				x = (z.imag()-m1) / std::tan(_rad_);
				x = std::max(z.real() - x, ax1.o);
				z = {x, m1};
			}
			if (i2 >= ax2.n) {
				x = (z.real()-m2) * std::tan(_rad_);
				x = std::min(z.imag() - x, ax1.o);
				z = {m2, x};
			}
			zmap.data()[i] = z;
			umap.data()[i] = getInCoord()[i] * _rot_;
		}
	}

	void mapJacob() {
		for (int i=0; i < zcoord.size(); ++i) {
			jacob.data()[i] = 1.f;
		}
	};


	float _rad_;
	float max1 = 0; float max2 = 0;
	float o1 = 0; float o2 = 0;
	std::complex<float> _rot_, _irot_;
	std::vector<std::complex<float>> v;
	std::vector<std::complex<float>> u;
};

}
