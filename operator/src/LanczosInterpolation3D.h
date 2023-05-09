#pragma once
#include "complex3DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <lanczos_kernel.h>
#include "Interpolation3D.h"

namespace SEP {

class LanczosInterpolation3D : public Interpolation3D {

public:
	LanczosInterpolation3D(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data,
		std::vector<float> a, std::vector<float> taper_perc) :
	Interpolation3D(model,data) {
		for (int i=0; i<3; ++i) {
			_a_.push_back(a[i] * model->getHyper()->getAxis(i+1).d);
			int nf = 2*_a_[i] / data->getHyper()->getAxis(i+1).d + 1;
			na.push_back((nf-1)/2);
			half.push_back(na[i]*data->getHyper()->getAxis(i+1).d);
		}
		filter3 = buildFilter(ax_d[2],ax_m[2],_a_[2], taper_perc[2]);
		filter2 = buildFilter(ax_d[1],ax_m[1],_a_[1], taper_perc[1]);
		filter1 = buildFilter(ax_d[0],ax_m[0],_a_[0], taper_perc[0]);
	}


private:

	std::shared_ptr<float2D> buildFilter(axis ax_d, axis ax_m, float& a, float tap_perc) {
		int nf = 2*a / ax_d.d + 1;
		int	na = (nf-1)/2;
		float half = na*ax_d.d;

		float new_d = (1.f + tap_perc) * ax_m.d;
		std::shared_ptr<float2D> filter = std::make_shared<float2D>(boost::extents[ax_d.n][ax_m.n]);
		for (int i=0; i<ax_d.n; ++i) {
			for (int j=0; j<ax_m.n; ++j) {
				float x = i*ax_d.d-j*ax_m.d;
				float xx = i*ax_d.d+j*ax_m.d;
				if (std::abs(x) <= half) {
					if (x==0) (*filter)[i][j] = 1;
					else (*filter)[i][j] = new_d * a * std::sin(PI*x/new_d) * std::sin(PI*x/a) / (PI*PI*x*x);
				}
			}
		}
		for (int i=0; i<na; ++i) {
			for (int j=0; j<ax_m.n; ++j) {
				float x = i*ax_d.d-j*ax_m.d;
				float xx = i*ax_d.d+j*ax_m.d;
				if (std::abs(xx) < half && x != xx) {
					(*filter)[i][j] += new_d * a * std::sin(PI*xx/new_d) * std::sin(PI*xx/a) / (PI*PI*xx*xx);
					(*filter)[ax_d.n-i-1][ax_m.n-j-1] += new_d * a * std::sin(PI*xx/new_d) * std::sin(PI*xx/a) / (PI*PI*xx*xx);
				}
			}
		}
		return filter;
	}

	std::vector<int> na;
	std::vector<float> half, _a_;
};

}
