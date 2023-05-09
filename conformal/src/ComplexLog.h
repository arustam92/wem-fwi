#pragma once

#include "float1DReg.h"
#include "complex1DReg.h"
#include <complex2DReg.h>
#include "ConformalMap.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace SEP {

class ComplexLog : public ConformalMap {
public:
	// eps = e^(-eps) - circle around zero
	ComplexLog(std::shared_ptr<hypercube> inHyper, float eps=0.01) : ConformalMap(inHyper) {
		_eps_ = std::log(eps);
		allocateOutSpace();
		complexify();
		mapCoord();
		mapJacob();
	}

	~ComplexLog() {};

	std::shared_ptr<complex2DReg> forward(std::shared_ptr<complex2DReg> inSpace);
	std::shared_ptr<complex2DReg> inverse(std::shared_ptr<complex2DReg> outSpace);
	void forward(const std::shared_ptr<complex2DReg>& inSpace, std::shared_ptr<complex2DReg>& outSpace);
	void inverse(std::shared_ptr<complex2DReg>& inSpace, const std::shared_ptr<complex2DReg>& outSpace);
	std::shared_ptr<Point> forward(std::vector<int>& re, std::vector<int>& im);
	std::shared_ptr<Point> forward(int& re, int& im);

private:
	/*
		Constructs domain for action of conformal mapping.
		Assumes "canonical domain" as a map of unit circle.
		Have to prepare the input domain accordingly.
	*/
	void allocateOutSpace() {
		// fast = radians of the ray
		// slow = -eps - 0
		_outSpace_ = std::make_shared<complex2DReg>(axis(getInHyper()->getAxis(1).n, 0, 2*PI/(getInHyper()->getAxis(1).n-1)),
																							axis(getInHyper()->getAxis(2).n, _eps_, -_eps_/(getInHyper()->getAxis(2).n-1)));
	}

	// construct a map between coordinates of the 2 domains as z = inv(f)(u)
	void mapCoord(){
		for (int i=0; i < zcoord.size(); ++i){
			zmap.data()[i] = std::exp(getOutCoord()[i]);	// pointer to function
			constexpr std::complex<float> ip = {0,PI};
			std::complex<float> zz = getInCoord()[i];
			if (std::abs(zz) < std::exp(_eps_)) {
				umap.data()[i] = {_eps_,std::arg(zz)};
			}
			else umap.data()[i] = std::log(zz);

			if (umap.data()[i].imag() < 0) umap.data()[i] += 2.f*ip;

		}
}

		void mapJacob() {
			for (int i=0; i < zcoord.size(); ++i) {
			jacob.data()[i] = std::exp(getOutCoord()[i].real());
		}
	};


	float _eps_;
};

}
