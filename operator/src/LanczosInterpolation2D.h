#pragma once
#include "complex2DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <lanczos_kernel.h>

using namespace SEP;

namespace SEP {

class LanczosInterpolation2D : public Operator<complex2DReg,complex2DReg> {

public:
	LanczosInterpolation2D(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, std::vector<float> a) {
		for (int i=0; i<2; ++i) {
			_a_.push_back(a[i] * model->getHyper()->getAxis(i+1).d);
			int nf = 2*_a_[i] / data->getHyper()->getAxis(i+1).d + 1;
			na.push_back((nf-1)/2);
			half.push_back(na[i]*data->getHyper()->getAxis(i+1).d);
			ax_d.push_back(data->getHyper()->getAxis(i+1));
			ax_m.push_back(model->getHyper()->getAxis(i+1));
		}
		filter2 = buildFilter(ax_d[1],ax_m[1],_a_[1]);
		filter1 = buildFilter(ax_d[0],ax_m[0],_a_[0]);

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
		// });
	};

private:

	std::shared_ptr<float2D> buildFilter(axis ax_d, axis ax_m, float& a) {
		int nf = 2*a / ax_d.d + 1;
		int	na = (nf-1)/2;
		float half = na*ax_d.d;
		double PI = 4.*std::atan(1);

		std::shared_ptr<float2D> filter = std::make_shared<float2D>(boost::extents[ax_d.n][ax_m.n]);
		for (int i=0; i<ax_d.n; ++i) {
			for (int j=0; j<ax_m.n; ++j) {
				float x = i*ax_d.d - j*ax_m.d;
				float xx = i*ax_d.d + j*ax_m.d;
				if (std::abs(x) <= half) {
					if (x==0) (*filter)[i][j] = 1;
					else (*filter)[i][j] = ax_m.d * a * std::sin(PI*x/ax_m.d) * std::sin(PI*x/a) / (PI*PI*x*x);
				}
			}
		}
		for (int i=0; i<na; ++i) {
			for (int j=0; j<ax_m.n; ++j) {
				float x = i*ax_d.d - j*ax_m.d;
				float xx = i*ax_d.d + j*ax_m.d;
				if (std::abs(xx) < half && x != xx) {
					(*filter)[i][j] += ax_m.d * a * std::sin(PI*xx/ax_m.d) * std::sin(PI*xx/a) / (PI*PI*xx*xx);
					(*filter)[ax_d.n-i-1][ax_m.n-j-1] += ax_m.d * a * std::sin(PI*xx/ax_m.d) * std::sin(PI*xx/a) / (PI*PI*xx*xx);
				}
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
