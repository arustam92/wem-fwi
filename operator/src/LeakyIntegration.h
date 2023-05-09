#pragma once
#include "complex4DReg.h"
#include "complex3DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

namespace SEP {

class LeakyIntegration : public Operator<complex4DReg,complex4DReg> {
public:
	LeakyIntegration(std::shared_ptr<complex4DReg> model, float alpha = 0.9) {
		_alpha_ = alpha;
	}

	void forward(std::shared_ptr<complex4DReg> model, std::shared_ptr<complex4DReg> data, bool add) {
		if (!add) data->scale(0);
		int n4 = model->getHyper()->getAxis(4).n;
		int n3 = model->getHyper()->getAxis(3).n;
		int n2 = model->getHyper()->getAxis(2).n;
		int n1 = model->getHyper()->getAxis(1).n;
		tbb::parallel_for(tbb::blocked_range<int>(0,n3),
			[=](const tbb::blocked_range<int> &r3) {
				tbb::parallel_for(tbb::blocked_range<int>(0,n2),
					[=](const tbb::blocked_range<int> &r2) {
						tbb::parallel_for(tbb::blocked_range<int>(0,n1),
							[=](const tbb::blocked_range<int> &r1) {
								std::shared_ptr<complex3DReg> prev = std::make_shared<complex3DReg>(model->getHyper()->getAxis(1),model->getHyper()->getAxis(2),model->getHyper()->getAxis(3));
		for (int i4 = 0; i4 < n4; ++i4) {
			for (int i3 = r3.begin(); i3 < r3.end(); ++i3) {

				for (int i2 = r2.begin(); i2 < r2.end(); ++i2) {

					for (int i1 = r1.begin(); i1 < r1.end(); ++i1) {
						(*data->_mat)[i4][i3][i2][i1] += (*model->_mat)[i4][i3][i2][i1] + _alpha_ * (*prev->_mat)[i3][i2][i1];
						(*prev->_mat)[i3][i2][i1] = (*data->_mat)[i4][i3][i2][i1] ;
					}
				}
			}
		}
		});
		});
		});
	};
	void adjoint(std::shared_ptr<complex4DReg> model, std::shared_ptr<complex4DReg> data, bool add) {
		if (!add) model->scale(0);
		int n4 = model->getHyper()->getAxis(4).n;
		int n3 = model->getHyper()->getAxis(3).n;
		int n2 = model->getHyper()->getAxis(2).n;
		int n1 = model->getHyper()->getAxis(1).n;
		for (int i4 = n4-1; i4 > 0; --i4) {
		tbb::parallel_for(tbb::blocked_range<int>(0,n3),
			[=](const tbb::blocked_range<int> &r3) {
			for (int i3 = r3.begin(); i3 < r3.end(); ++i3) {
			tbb::parallel_for(tbb::blocked_range<int>(0,n2),
				[=](const tbb::blocked_range<int> &r2) {
				for (int i2 = r2.begin(); i2 < r2.end(); ++i2) {
				tbb::parallel_for(tbb::blocked_range<int>(0,n1),
					[=](const tbb::blocked_range<int> &r1) {
					for (int i1 = r1.begin(); i1 < r1.end(); ++i1) {
						(*model->_mat)[i4][i3][i2][i1] += (*data->_mat)[i4][i3][i2][i1] + _alpha_ * (*prev->_mat)[i3][i2][i1];
						(*prev->_mat)[i3][i2][i1] = (*model->_mat)[i4][i3][i2][i1] ;
					}});
				}});
			}});
		}
	};

private:
	float _alpha_;
	std::shared_ptr<complex3DReg> prev;

};

}
