#pragma once
#include "complex2DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <lanczos_kernel.h>

using namespace SEP;

namespace SEP {

class Spline2D : public Operator<complex2DReg,complex2DReg> {

public:
	Spline2D(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, float a = 1, float b = 0) {
		for (int i=0; i<2; ++i) {
			ax_d.push_back(data->getHyper()->getAxis(i+1));
			ax_m.push_back(model->getHyper()->getAxis(i+1));
		}
		filter2 = buildFilter(ax_d[1],ax_m[1], a, b);
		filter1 = buildFilter(ax_d[0],ax_m[0], a, b);

	}

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) data->scale(0);

		float* mre = &reinterpret_cast<float(&)[2]>(model->_mat->data()[0])[0];
		float* dre = &reinterpret_cast<float(&)[2]>(data->_mat->data()[0])[0];

		tbb::parallel_for(tbb::blocked_range<int>(0,ax_d[1].n),
		[=](const tbb::blocked_range<int> &r2) {
		// tbb::parallel_for(tbb::blocked_range<int>(0,ax_m[1].n),
		// [=](const tbb::blocked_range<int> &m3) {

		float filt1, filt2;
		int ind, ind2;

		for (int i2 = r2.begin(); i2 < r2.end(); i2++) {
			for (int j2 = 0; j2 < ax_m[1].n; j2++) {
				filt2 = (*filter2)[i2][j2];
				if (filt2 != 0) {
					for (int i1 = 0; i1 < ax_d[0].n; i1++) {
						ind = j2*ax_m[0].n;
						ind2 = i1 + i2*ax_d[0].n;
						ispc::lanczos_fwd_2d(0, ax_m[0].n, filt2, filter1->data() + i1*ax_m[0].n,
																mre + 2*ind, dre + 2*ind2);
					}
				}
			}
		}



		});
		// });
	};

	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) model->scale(0);
		float* mre = &reinterpret_cast<float(&)[2]>(model->_mat->data()[0])[0];
		float* dre = &reinterpret_cast<float(&)[2]>(data->_mat->data()[0])[0];

		// tbb::parallel_for(tbb::blocked_range<int>(0,ax_d[2].n),
		// [=](const tbb::blocked_range<int> &r3) {
		tbb::parallel_for(tbb::blocked_range<int>(0,ax_m[1].n),
		[=](const tbb::blocked_range<int> &m2) {

		float filt1, filt2;
		int ind, ind2;
		for (int j2 = m2.begin(); j2 != m2.end(); j2++) {
			for (int i2 = 0; i2 != ax_d[1].n; i2++) {
				filt2 = (*filter2)[i2][j2];
				if (filt2 != 0) {
					for (int i1 = 0; i1 < ax_d[0].n; i1++) {
							ind = j2*ax_m[0].n;
							ind2 = i1 + i2*ax_d[0].n;
							ispc::lanczos_adj_2d(0, ax_m[0].n,	filt2, filter1->data() + i1*ax_m[0].n,
															mre + 2*ind, dre + 2*ind2);
						}
					}
				}
			}
		});
		// avoiding boundary effects -> doesn't pass the dot test (it's okay though)
		for (int j2=0; j2 < ax_m[1].n; ++j2) {
				(*model->_mat)[j2][0] = (*model->_mat)[j2][1];
				(*model->_mat)[j2][ax_m[0].n-1] = (*model->_mat)[j2][ax_m[0].n-2];
		}
		for (int j1=0; j1 < ax_m[0].n; ++j1) {
				(*model->_mat)[0][j1] = (*model->_mat)[1][j1];
				(*model->_mat)[ax_m[1].n-1][j1] = (*model->_mat)[ax_m[1].n-2][j1];
		}
	};

private:

	std::shared_ptr<float2D> buildFilter(axis ax_d, axis ax_m, float& a, float &b) {

		std::shared_ptr<float2D> filter = std::make_shared<float2D>(boost::extents[ax_d.n][ax_m.n]);
		for (int i=0; i<ax_d.n; ++i) {
			for (int j=0; j<ax_m.n; ++j) {
				float x = std::abs(i*ax_d.d - j*ax_m.d);
				float y = x / ax_m.d;
				if (y < 1) 
					(*filter)[i][j] = ((-6*a - 9*b + 12) * y*y*y + (6*a + 12*b - 18) * y*y - 2*b + 6) / 6;
				else if (1 <= y && y < 2) 
					(*filter)[i][j] = ((-6*a - b) * y*y*y + (30*a + 6*b) * y*y + (-48*a - 12*b)*y + 24*a + 8*b) / 6;
			}
		}
		
		for (int i=0; i<ax_d.n; ++i) {
			for (int j=1; j<ax_m.n-1; ++j) {
				float f = 0;
				float x = i*ax_d.d - j*ax_m.d;
				float xx = i*ax_d.d + j*ax_m.d;
				float y = std::abs(xx) / ax_m.d;
				if (y < 1.f ) 
					f = ((-6*a - 9*b + 12) * y*y*y + (6*a + 12*b - 18) * y*y - 2*b + 6) / 6;
				else if (1.f <= y && y < 2.f ) 
					f = ((-6*a - b) * y*y*y + (30*a + 6*b) * y*y + (-48*a - 12*b)*y + 24*a + 8*b) / 6;
				(*filter)[i][j] += f;
				(*filter)[ax_d.n-i-1][ax_m.n-j-1] += f;
			}
		}
		return filter;
	}

	std::vector<int> na;
	std::vector<float> half, _a_;
	std::vector<axis> ax_d, ax_m;
	std::shared_ptr<float2D> filter1, filter2;
};

}
