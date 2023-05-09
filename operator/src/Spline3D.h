#pragma once
#include "complex2DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <lanczos_kernel.h>
#include "Interpolation3D.h"

using namespace SEP;

namespace SEP {

class Spline3D : public Interpolation3D{

public:
	Spline3D(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, float a, float b, std::vector<float> taper_perc) :
	Interpolation3D(model,data) {

		filter3 = buildFilter(ax_d[2],ax_m[2], a, b, taper_perc[2]);
		filter2 = buildFilter(ax_d[1],ax_m[1], a, b, taper_perc[1]);
		filter1 = buildFilter(ax_d[0],ax_m[0], a, b, taper_perc[0]);

	}


private:

	std::shared_ptr<float2D> buildFilter(axis ax_d, axis ax_m, float& a, float &b, float tap_perc) {

		float new_d = (1.f + tap_perc) * ax_m.d;
		std::shared_ptr<float2D> filter = std::make_shared<float2D>(boost::extents[ax_d.n][ax_m.n]);
		for (int i=0; i<ax_d.n; ++i) {
			for (int j=0; j<ax_m.n; ++j) {
				float x = std::abs(i*ax_d.d - j*ax_m.d);
				float y = x / new_d;
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
				float y = std::abs(xx) / new_d;
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

};

}
