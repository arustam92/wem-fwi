#pragma once

#include "float1DReg.h"
#include "complex1DReg.h"
#include <complex2DReg.h>
#include "ConformalMap.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace SEP {

class MapChain : public ConformalMap {
public:
	// Chain of conformal maps
	MapChain(std::shared_ptr<ConformalMap> map1, std::shared_ptr<ConformalMap> map2) : _map1_(map1), _map2_(map2) {}
	~MapChain() {};

	std::shared_ptr<complex2DReg> forward(std::shared_ptr<complex2DReg> inSpace) {
		std::shared_ptr<complex2DReg> out = _map1_->forward(inSpace);
		out = _map2_->forward(out);
		return out;
	}

	std::shared_ptr<complex2DReg> inverse(std::shared_ptr<complex2DReg> outSpace) {
		std::shared_ptr<complex2DReg> in = _map2_->inverse(outSpace);
		in = _map1_->inverse(in);
		return in;
	}

	void forward(const std::shared_ptr<complex2DReg>& inSpace, std::shared_ptr<complex2DReg>& outSpace) {
		auto temp = std::make_shared<complex2DReg>(_map1_->getOutHyper());
		_map1_->forward(inSpace, temp);
		_map2_->forward(temp, outSpace);
	}
	void inverse(std::shared_ptr<complex2DReg>& inSpace, const std::shared_ptr<complex2DReg>& outSpace) {
		auto temp = std::make_shared<complex2DReg>(_map2_->getInHyper());
		_map2_->inverse(temp, outSpace);
		_map1_->inverse(inSpace, temp);
	}

	std::shared_ptr<Point> forward(std::vector<int>& re, std::vector<int>& im) {
		std::shared_ptr<Point> point = _map1_->forward(re,im);
		point = _map2_->forward(point->real,point->imag);
		return point;
	};

	std::shared_ptr<Point> forward(int& re, int& im) {
		std::shared_ptr<Point> point = _map1_->forward(re,im);
		point = _map2_->forward(point->real,point->imag);
		return point;
	};


	virtual void scaleByJacobian(std::shared_ptr<complex2DReg> inSpace){
		_map1_->scaleByJacobian(inSpace);
		_map2_->scaleByJacobian(inSpace);
	};

private:
	std::shared_ptr<ConformalMap> _map1_, _map2_;
};

}
