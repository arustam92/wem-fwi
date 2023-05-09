#pragma once
#include "complex4DReg.h"
#include "complex3DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace SEP {

class LinearInterpolation3D : public Operator<complex4DReg,complex4DReg> {
public:
	LinearInterpolation3D(std::shared_ptr<complex4DReg> data, std::vector<std::vector<float>> pts) {
		for (int i=0; i < data->getHyper()->getNdim();++i) ax.push_back(data->getHyper()->getAxis(i+1));
		_pts_ = pts;
	}

	void forward(std::shared_ptr<complex4DReg> model, std::shared_ptr<complex4DReg> data, bool add) {
		if (!add) data->scale(0);

		tbb::parallel_for(tbb::blocked_range<int>(0,_pts_[0].size()-1),
			[=](const tbb::blocked_range<int> &r4) {
		tbb::parallel_for(tbb::blocked_range<int>(0,_pts_[1].size()-1),
				[=](const tbb::blocked_range<int> &r2) {
		tbb::parallel_for(tbb::blocked_range<int>(0,_pts_[2].size()-1),
			[=](const tbb::blocked_range<int> &r1) {

		for (int ii4 = r4.begin(); ii4 < r4.end(); ++ii4) {
			int i4cur = (_pts_[2][ii4] - ax[3].o)/ax[3].d + 0.5;
			int i4next = (_pts_[2][ii4+1] - ax[3].o)/ax[3].d + 0.5;
		for (int ii3 = r4.begin(); ii3 < r4.end(); ++ii3) {
			int i3cur = (_pts_[1][ii3] - ax[3].o)/ax[3].d + 0.5;
			int i3next = (_pts_[1][ii3+1] - ax[3].o)/ax[3].d + 0.5;

		for (int i4 = i4cur; i4 < i4next; ++i4) {
			float alpha = std::max(0.f,(_pts_[ii4+1] - ax[3].o - i4*ax[3].d) / (_pts_[ii4+1] - _pts_[ii4] ));
			for (int i3 = 0; i3 < ax[2].n; ++i3) {
				for (int i2 = 0; i2 < ax[1].n; ++i2) {
					for (int i1 = 0; i1 < ax[0].n; ++i1) {
						(*data->_mat)[i4][i3][i2][i1] += alpha * (*model->_mat)[ii4][i3][i2][i1] + (1-alpha) * (*model->_mat)[ii4+1][i3][i2][i1];
					}
				}
			}
		}
	}});
	};
	void adjoint(std::shared_ptr<complex4DReg> model, std::shared_ptr<complex4DReg> data, bool add) {
		if (!add) model->scale(0);

		tbb::parallel_for(tbb::blocked_range<int>(0,_pts_.size()-1),
			[=](const tbb::blocked_range<int> &r4) {
		for (int j4 = r4.begin(); j4 < r4.end(); ++j4) {
			int icur = (_pts_[j4] - ax[3].o)/ax[3].d + 0.5;
			int inext = (_pts_[j4+1] - ax[3].o)/ax[3].d + 0.5;
		for (int i4 = icur; i4 < inext; ++i4) {
			float alpha = std::max(0.f,(_pts_[j4+1] - ax[3].o - i4*ax[3].d) / (_pts_[j4+1] - _pts_[j4] ));
			for (int i3 = 0; i3 < ax[2].n; ++i3) {
				for (int i2 = 0; i2 < ax[1].n; ++i2) {
					for (int i1 = 0; i1 < ax[0].n; ++i1) {
						(*model->_mat)[j4][i3][i2][i1] += alpha * (*data->_mat)[i4][i3][i2][i1];
						(*model->_mat)[j4+1][i3][i2][i1] += (1-alpha) * (*data->_mat)[i4][i3][i2][i1];
					}
				}
			}
		}
	}});
	};

private:
	std::vector<axis> ax;
	std::vector<std::vector<float>> _pts_;

};

}
