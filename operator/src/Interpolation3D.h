#pragma once
#include "complex3DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <lanczos_kernel.h>

namespace SEP {

double PI = 4.*std::atan(1);

class Interpolation3D : public Operator<complex3DReg,complex3DReg> {

public:
	Interpolation3D(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data) {
		for (int i=0; i<3; ++i) {
			ax_d.push_back(data->getHyper()->getAxis(i+1));
			ax_m.push_back(model->getHyper()->getAxis(i+1));
		}
	}

	void forward(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, bool add) {
		if (!add) data->scale(0);

		float* mre = &reinterpret_cast<float(&)[2]>(model->_mat->data()[0])[0];
		float* dre = &reinterpret_cast<float(&)[2]>(data->_mat->data()[0])[0];

		tbb::parallel_for(tbb::blocked_range<int>(0,ax_d[2].n),
		[=](const tbb::blocked_range<int> &r3) {
		// tbb::parallel_for(tbb::blocked_range<int>(0,ax_m[1].n),
		// [=](const tbb::blocked_range<int> &m3) {

		float filt1, filt2, filt3;
		int ind, ind2;

		for (int i3 = r3.begin(); i3 != r3.end(); i3++) {
			for (int j3 = 0; j3 != ax_m[2].n; j3++) {
				filt3 = (*filter3)[i3][j3];
				if (filt3 != 0) {
					for (int i2 = 0; i2 < ax_d[1].n; i2++) {
						for (int j2 = 0; j2 < ax_m[1].n; j2++) {
							filt2 = (*filter2)[i2][j2];
							if (filt2 != 0) {
								for (int i1 = 0; i1 < ax_d[0].n; i1++) {
									ind = j2*ax_m[0].n + j3*ax_m[1].n*ax_m[0].n;
									ind2 = i1 + i2*ax_d[0].n + i3*ax_d[0].n*ax_d[1].n;
									ispc::lanczos_fwd(0, ax_m[0].n, filt3, filt2, filter1->data() + i1*ax_m[0].n,
																			mre + 2*ind, dre + 2*ind2);
								}
							}
						}
					}
				}
			}
		}
		});
		// });
	};

	void adjoint(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, bool add) {
		if (!add) model->scale(0);
		float* mre = &reinterpret_cast<float(&)[2]>(model->_mat->data()[0])[0];
		float* dre = &reinterpret_cast<float(&)[2]>(data->_mat->data()[0])[0];

		// tbb::parallel_for(tbb::blocked_range<int>(0,ax_d[2].n),
		// [=](const tbb::blocked_range<int> &r3) {
		tbb::parallel_for(tbb::blocked_range<int>(0,ax_m[2].n),
		[=](const tbb::blocked_range<int> &m3) {

		float filt1, filt2, filt3;
		int ind, ind2;
		for (int j3 = m3.begin(); j3 != m3.end(); j3++) {
		for (int i3 = 0; i3 != ax_d[2].n; i3++) {
				filt3 = (*filter3)[i3][j3];
				if (filt3 != 0) {
					for (int i2 = 0; i2 < ax_d[1].n; i2++) {
						for (int j2 = 0; j2 < ax_m[1].n; j2++) {
							filt2 = (*filter2)[i2][j2];
							if (filt2 != 0) {
								for (int i1 = 0; i1 < ax_d[0].n; i1++) {
										ind = j2*ax_m[0].n + j3*ax_m[1].n*ax_m[0].n;
										ind2 = i1 + i2*ax_d[0].n + i3*ax_d[0].n*ax_d[1].n;
										ispc::lanczos_adj(0, ax_m[0].n,	filt3, filt2, filter1->data() + i1*ax_m[0].n,
																		mre + 2*ind, dre + 2*ind2);
									}
								}
							}
						}
					}
				}
			}
		});
		// });
	};

protected:

	std::vector<axis> ax_d, ax_m;
	std::shared_ptr<float2D> filter1, filter2, filter3;
};

}
