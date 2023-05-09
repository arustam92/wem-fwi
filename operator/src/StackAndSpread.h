#pragma once

#include <complex1DReg.h>
#include <complex2DReg.h>
#include <complex3DReg.h>
#include <Operator.h>

namespace SEP {

class Stack : public Operator<complex3DReg,complex2DReg>
{
public:
	Stack()  {;}

	void forward(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) data->scale(0);
		std::vector<axis> ax = model->getHyper()->getAxes();

		for (int i3=0; i3 < ax[2].n; ++i3) {
			for (int i2=0; i2 < ax[1].n; ++i2) {
				for (int i1=0; i1 < ax[0].n; ++i1) {
					(*data->_mat)[i2][i1] += 1.f/ax[2].n * (*model->_mat)[i3][i2][i1].real() ;
				}
			}
		}

	};

	void adjoint(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) model->scale(0);
		std::vector<axis> ax = model->getHyper()->getAxes();

		for (int i3=0; i3 < ax[2].n; ++i3) {
			for (int i2=0; i2 < ax[1].n; ++i2) {
				for (int i1=0; i1 < ax[0].n; ++i1) {
					(*model->_mat)[i3][i2][i1] += 1.f/ax[2].n * (*data->_mat)[i2][i1].real();
				}
			}
		}

	};
};


}
