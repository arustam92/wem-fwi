
#pragma once
#include <hypercube.h>
#include <float2DReg.h>
#include <complex2DReg.h>
#include <algorithm>

namespace SEP {

class Illumination
{
public:
	Illumination(std::shared_ptr<hypercube> hyper) {
		nz = hyper->getAxis(2).n;
		nx = hyper->getAxis(1).n;
		_illum = std::make_shared<float2DReg>(hyper);
		maxAmp = 0;
		eps = 0.1;
	};

	void accumulate(std::shared_ptr<complex2DReg> bg_wfld) {
		for (size_t iz=0; iz < nz; ++iz) {
			for (size_t ix = 0; ix < nx; ++ix) {
				(*_illum->_mat)[iz][ix] += std::abs((*bg_wfld->_mat)[iz][ix]);
			}
		}
	}

	void compensate(std::shared_ptr<float2DReg> in) {

		maxAmp = *std::max_element(_illum->_mat->data(), _illum->_mat->data()+nz*nx);

		for (size_t iz=0; iz < nz; ++iz) {
			for (size_t ix = 0; ix < nx; ++ix) {
				(*in->_mat)[iz][ix] /= ((*_illum->_mat)[iz][ix] + eps*maxAmp);
			}
		}
		_illum->scale(0);
	}

private:
	int nz,nx;
	std::shared_ptr<float2DReg> _illum;
	float maxAmp, eps;
};

}

// // simple cross correlation
// void ic(std::shared_ptr<float2DReg> image, std::shared_ptr<complex2DReg> wfld1, std::shared_ptr<complex2DReg> wfld2);
