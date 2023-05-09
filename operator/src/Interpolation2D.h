#pragma once
#include "complex2DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <lanczos_kernel.h>

using namespace SEP;

namespace SEP {

class Interpolation2D : public Operator<complex2DReg,complex2DReg> {

public:
	Interpolation2D() {;}

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

	std::vector<axis> ax_d, ax_m;
	std::shared_ptr<float2D> filter1, filter2;
};

}
