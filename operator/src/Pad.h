#pragma once

#include <complex1DReg.h>
#include <complex2DReg.h>
#include <complex3DReg.h>
#include <Operator.h>


namespace SEP {

class Pad {
public:
	Pad(const std::shared_ptr<hypercube>& inSpace, int beg, int end, bool extend = 0)  {
		_beg_ = beg;
		_end_ = end;
		_extend_ = extend;
		_ax_ = inSpace->getAxes();
		std::vector<axis> ax_pad(_ax_.size());
		for (int i=1; i < _ax_.size(); ++i)
			ax_pad[i] = _ax_[i];
 		int n1 = _ax_[0].n + beg + end;
		float o1 = _ax_[0].o - beg*_ax_[0].d;
		float d1 = _ax_[0].d;
		ax_pad[0] = axis(n1,o1,d1);
		_outHyper_ = std::make_shared<hypercube> (ax_pad);
	}

	inline std::shared_ptr<hypercube>& getPaddedHyper() {return _outHyper_;}

	void forward(const std::shared_ptr<float1DReg> model, std::shared_ptr<float1DReg> data, bool add) {
		if (!add) data->scale(0);
		axis ax1 = model->getHyper()->getAxis(1);
		float bval = 0;
		float eval = 0;

		for (int i1=0; i1 < ax1.n; ++i1) {
			(*data->_mat)[i1+_beg_] += (*model->_mat)[i1];
		}
		if (_extend_){
			bval = (*data->_mat)[_beg_];
			eval = (*data->_mat)[ax1.n-1+_beg_];
		}
		for (int i1=0; i1 < _beg_; ++i1)
			(*data->_mat)[i1] = bval;
		for (int i1=0; i1 < _end_; ++i1)
			(*data->_mat)[ax1.n + i1 + _beg_] = eval;

	};

	void forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
		if (!add) data->scale(0);
		axis ax1 = model->getHyper()->getAxis(1);
		std::complex<float> bval = 0;
		std::complex<float> eval = 0;

		for (int i1=0; i1 < ax1.n; ++i1) {
			(*data->_mat)[i1+_beg_] += (*model->_mat)[i1];
		}
		if (_extend_){
			bval = (*data->_mat)[_beg_];
			eval = (*data->_mat)[ax1.n-1+_beg_];
		}
		for (int i1=0; i1 < _beg_; ++i1)
			(*data->_mat)[i1] = bval;
		for (int i1=0; i1 < _end_; ++i1)
			(*data->_mat)[ax1.n + i1 + _beg_] = eval;
	};
	// not really adjoint but just un-padding
	void adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add) {
		if (!add) model->scale(0);
		axis ax1 = model->getHyper()->getAxis(1);
			for (int i1=0; i1 < ax1.n; ++i1) {
				 (*model->_mat)[i1] += (*data->_mat)[i1+_beg_];
			}
	};

	void forward(const std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) data->scale(0);
		auto ax = model->getHyper()->getAxes();
		std::complex<float> bval = 0;
		std::complex<float> eval = 0;

		for (int i2=0; i2 < ax[1].n; ++i2) {
			for (int i1=0; i1 < ax[0].n; ++i1) {
				(*data->_mat)[i2][i1+_beg_] += (*model->_mat)[i2][i1];
			}
			if (_extend_){
				bval = (*data->_mat)[i2][_beg_];
				eval = (*data->_mat)[i2][ax[0].n-1+_beg_];
			}
			for (int i1=0; i1 < _beg_; ++i1)
				(*data->_mat)[i2][i1] = bval;
			for (int i1=0; i1 < _end_; ++i1)
				(*data->_mat)[i2][ax[0].n + i1 + _beg_] = eval;
		}
	};

	void adjoint(std::shared_ptr<complex2DReg> model, const std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) model->scale(0);
		auto ax = model->getHyper()->getAxes();

		for (int i2=0; i2 < ax[1].n; ++i2) {
			for (int i1=0; i1 < ax[0].n; ++i1) {
				 (*model->_mat)[i2][i1] += (*data->_mat)[i2][i1+_beg_];
			}
		}
	};

	void forward(const std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, bool add) {
		if (!add) data->scale(0);
		auto ax = model->getHyper()->getAxes();
		std::complex<float> bval = 0;
		std::complex<float> eval = 0;

	for (int i3=0; i3 < ax[2].n; ++i3) {
		for (int i2=0; i2 < ax[1].n; ++i2) {
			for (int i1=0; i1 < ax[0].n; ++i1) {
				(*data->_mat)[i3][i2][i1+_beg_] += (*model->_mat)[i3][i2][i1];
			}
			if (_extend_){
				bval = (*data->_mat)[i3][i2][_beg_];
				eval = (*data->_mat)[i3][i2][ax[0].n-1+_beg_];
			}
			for (int i1=0; i1 < _beg_; ++i1)
				(*data->_mat)[i3][i2][i1] = bval;
			for (int i1=0; i1 < _end_; ++i1)
				(*data->_mat)[i3][i2][ax[0].n + i1 + _beg_] = eval;
		}
		}
	};

	void adjoint(std::shared_ptr<complex3DReg> model, const std::shared_ptr<complex3DReg> data, bool add) {
		if (!add) model->scale(0);
		auto ax = model->getHyper()->getAxes();

		for (int i3=0; i3 < ax[2].n; ++i3) {
		for (int i2=0; i2 < ax[1].n; ++i2) {
			for (int i1=0; i1 < ax[0].n; ++i1) {
				 (*model->_mat)[i3][i2][i1] += (*data->_mat)[i3][i2][i1+_beg_];
			}
		}
	}
	};

	void forward(const std::shared_ptr<float3DReg> model, std::shared_ptr<float3DReg> data, bool add) {
		if (!add) data->scale(0);
		auto ax = model->getHyper()->getAxes();
		float bval = 0;
		float eval = 0;

		for (int i3=0; i3 < ax[2].n; ++i3) {
			for (int i2=0; i2 < ax[1].n; ++i2) {
				for (int i1=0; i1 < ax[0].n; ++i1) {
					(*data->_mat)[i3][i2][i1+_beg_] += (*model->_mat)[i3][i2][i1];
				}
				if (_extend_){
					bval = (*data->_mat)[i3][i2][_beg_];
					eval = (*data->_mat)[i3][i2][ax[0].n-1+_beg_];
				}
				for (int i1=0; i1 < _beg_; ++i1)
					(*data->_mat)[i3][i2][i1] = bval;
				for (int i1=0; i1 < _end_; ++i1)
					(*data->_mat)[i3][i2][ax[0].n + i1 + _beg_] = eval;
			}
			}
	};

	void adjoint(std::shared_ptr<float3DReg> model, const std::shared_ptr<float3DReg> data, bool add) {
		if (!add) model->scale(0);
		auto ax = model->getHyper()->getAxes();

		for (int i3=0; i3 < ax[2].n; ++i3) {
		for (int i2=0; i2 < ax[1].n; ++i2) {
			for (int i1=0; i1 < ax[0].n; ++i1) {
				 (*model->_mat)[i3][i2][i1] += (*data->_mat)[i3][i2][i1+_beg_];
			}
		}
	}
	};

private:

	std::shared_ptr<hypercube> _outHyper_;
	int _beg_, _end_;
	bool _extend_;
	std::vector<axis> _ax_;
};


}
