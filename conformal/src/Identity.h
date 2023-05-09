#pragma once

#include "float1DReg.h"
#include "complex1DReg.h"
#include <complex2DReg.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include "ConformalMap.h"

namespace SEP {

/*
	Map g: Z -> U
*/
class Identity : public ConformalMap{
public:
	Identity() {}
	~Identity() {}

	inline std::shared_ptr<complex2DReg> forward(std::shared_ptr<complex2DReg> inSpace) {
		return inSpace;};
	inline std::shared_ptr<complex2DReg> inverse(std::shared_ptr<complex2DReg> outSpace) {return outSpace;};

	inline void forward(const std::shared_ptr<complex2DReg>& inSpace, std::shared_ptr<complex2DReg>& outSpace) {outSpace = inSpace->clone();};
	inline void inverse(std::shared_ptr<complex2DReg>& inSpace, const std::shared_ptr<complex2DReg>& outSpace) {inSpace = outSpace->clone();};

	inline std::shared_ptr<ConformalMap::Point> forward(int& re, int& im) {
		std::shared_ptr<Point> point = std::make_shared<Point>();
		point->real.push_back(re);
		point->imag.push_back(im);
		return point;
	};
	inline std::shared_ptr<ConformalMap::Point> forward(std::vector<int>& re, std::vector<int>& im) {
		std::shared_ptr<Point> point = std::make_shared<Point>();
		point->real = re;
		point->imag = im;
		return point;
	};

	inline void scaleByJacobian(std::shared_ptr<complex2DReg> inSpace) {return;};


};


}
