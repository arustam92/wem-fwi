#pragma once

#include "float1DReg.h"
#include "complex1DReg.h"
#include <complex2DReg.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace SEP {

/*
	Map g: Z -> U
*/
class ConformalMap {
public:
	struct Point {
		std::vector<int> real;
		std::vector<int> imag;
	};
	ConformalMap() {}
	ConformalMap(std::shared_ptr<hypercube> inHyper) {
		_inSpace_ = std::make_shared<complex2DReg>(inHyper);
	}
	virtual ~ConformalMap() {}

	virtual std::shared_ptr<complex2DReg> forward(std::shared_ptr<complex2DReg> inSpace);
	virtual std::shared_ptr<complex2DReg> inverse(std::shared_ptr<complex2DReg> outSpace);

	virtual void forward(const std::shared_ptr<complex2DReg>& inSpace, std::shared_ptr<complex2DReg>& outSpace);
	virtual void inverse(std::shared_ptr<complex2DReg>& inSpace, const std::shared_ptr<complex2DReg>& outSpace);

	virtual std::shared_ptr<Point> forward(int& re, int& im) {

		std::shared_ptr<Point> point = std::make_shared<Point>();
		int n1 = getInHyper()->getAxis(1).n;
		axis outax1 = getOutHyper()->getAxis(1);
		axis outax2 = getOutHyper()->getAxis(2);
		std::complex<float> c;
		int ind = im + re*n1;
		c = getOutMap()[ind];
		int iz = (c.real() - outax2.o) / outax2.d;
		int ix = (c.imag() - outax1.o) / outax1.d;

		point->real.push_back(iz);
		point->imag.push_back(ix);
		return point;
	};

	virtual std::shared_ptr<Point> forward(std::vector<int>& re, std::vector<int>& im) {

		std::shared_ptr<Point> point = std::make_shared<Point>();
		int n1 = getInHyper()->getAxis(1).n;
		axis outax1 = getOutHyper()->getAxis(1);
		axis outax2 = getOutHyper()->getAxis(2);
		std::complex<float> c;
		for (int i=0; i < re.size(); ++i) {
			int ind = im[i] + re[i]*n1;
			c = getOutMap()[ind];
			int iz = (c.real() - outax2.o) / outax2.d;
			int ix = (c.imag() - outax1.o) / outax1.d;
			point->real.push_back(iz);
			point->imag.push_back(ix);
		}

		return point;
	};


		virtual void scaleByJacobian(std::shared_ptr<complex2DReg> inSpace){
			for (int i=0; i < jacob.size(); ++i) {
				float j = getJacobian()[i];
				inSpace->getVals()[i] *= j;
			}
		};

		std::shared_ptr<hypercube> getInHyper() {return _inSpace_->getHyper();}
		std::shared_ptr<hypercube> getOutHyper() {return _outSpace_->getHyper();}

protected:

	const std::complex<float>* getInMap() {return zmap.data();}
	const std::complex<float>* getOutMap() {return umap.data();}

	const std::complex<float>* getInCoord() {return zcoord.data();}
	const std::complex<float>* getOutCoord() {return ucoord.data();}
	const float* getJacobian() {return jacob.data();}

	void complexify() {
		zcoord.resize(boost::extents[getInHyper()->getAxis(1).n*getInHyper()->getAxis(2).n]);
		zmap.resize(boost::extents[zcoord.size()]);

		ucoord.resize(boost::extents[getInHyper()->getAxis(1).n*getInHyper()->getAxis(2).n]);
		umap.resize(boost::extents[ucoord.size()]);
		jacob.resize(boost::extents[umap.size()]);

		regularCoord(getInHyper(),zcoord.data());
		regularCoord(getOutHyper(),ucoord.data());
	}

	void regularCoord(std::shared_ptr<hypercube> hyper, std::complex<float>* coord) {
		axis ax2 = hyper->getAxis(2);
		axis ax1 = hyper->getAxis(1);

		float x,z;
		for (int i2=0; i2 < ax2.n; ++i2) {
			z = ax2.o + i2*ax2.d;
			for (int i1=0; i1 < ax1.n; ++i1) {
				x = ax1.o + i1*ax1.d;
				coord[i1 + i2*ax1.n] = {z,x};
			}
		}
	}

	std::shared_ptr<complex2DReg> _inSpace_,_outSpace_;
	complex1D zmap, umap;
	complex1D ucoord, zcoord;
	float1D jacob;

	constexpr static float PI = 4*std::atan(1.f);

};

}
