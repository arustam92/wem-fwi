#pragma once

#include "float1DReg.h"
#include "complex1DReg.h"
#include <complex2DReg.h>
#include "ConformalMap.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace SEP {

class ShiftAndScale : public ConformalMap {
public:
	// f(z) = a(z - b)
	ShiftAndScale(std::shared_ptr<hypercube> inHyper, std::complex<float>& a, std::complex<float>& b) : ConformalMap(inHyper) {
		_a_ = a; _b_ = b;
		allocateOutSpace();
		complexify();
		mapCoord();
		mapJacob();
	}

	~ShiftAndScale() {};

private:
	std::complex<float>& get_a() {return _a_;}
	std::complex<float>& get_b() {return _b_;}

	void allocateOutSpace() {
		axis ax1 = getInHyper()->getAxis(1);
		axis ax2 = getInHyper()->getAxis(2);
 		std::complex<float> min = {ax2.o,ax1.o};
		std::complex<float> max = {ax2.o + (ax2.n-1)*ax2.d, ax1.o + (ax1.n-1)*ax1.d};
		min = get_a() * (min - get_b());
		max = get_a() * (max - get_b());
		_outSpace_ = std::make_shared<complex2DReg>(axis(ax1.n, min.imag(), (max-min).imag()/(ax1.n-1)),
																							axis(ax2.n, min.real(), (max-min).real()/(ax2.n-1)));
	}

	// construct a map between coordinates of the 2 domains as z = inv(f)(u)
	void mapCoord(){
		for (int i=0; i < zcoord.size(); ++i){
			zmap.data()[i] = getOutCoord()[i] / get_a() + get_b();	// pointer to function
			umap.data()[i] = get_a() * (getInCoord()[i] - get_b());
		}
	}

	void mapJacob() {
		for (int i=0; i < zcoord.size(); ++i) {
			jacob.data()[i] = std::abs(1.f/get_a());
		}
	};


	std::complex<float> _a_, _b_;
};

}
