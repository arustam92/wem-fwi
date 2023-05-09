#pragma once
#include "complex3DReg.h"
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <FFT1.h>
#include "Pad.h"

using namespace SEP;

namespace SEP {

// f(x) -> f(x) + iH(f(x))
class Hilbert : public Operator<complex3DReg,complex3DReg> {

public:
	Hilbert(int npad) : _npad_(npad) {;}

	void forward(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex3DReg> data, bool add) {
		if (!add) data->scale(0);
		auto ax3 = model->getHyper()->getAxis(3);
		auto ax2 = model->getHyper()->getAxis(2);
		auto ax1 = model->getHyper()->getAxis(1);

		// tbb::parallel_for(tbb::blocked_range<int>(0,ax3.n),
		// [=](const tbb::blocked_range<int> &r3) {
		// 	tbb::parallel_for(tbb::blocked_range<int>(0,ax2.n),
		// 	[=](const tbb::blocked_range<int> &r2) {
				auto sig = std::make_shared<complex1DReg>(model->getHyper()->getAxis(1));
				auto pad = std::make_shared<Pad>(sig->getHyper(),0,_npad_);
				auto sig_pad = std::make_shared<complex1DReg>(pad->getPaddedHyper());
				auto fft = std::make_shared<FFT1>(sig_pad,"in",1,1);
				for (int i3=0; i3 != ax3.n; ++i3) {
					for (int i2=0; i2 != ax2.n; ++i2) {
						int off = i3*ax1.n*ax2.n + i2*ax1.n;
						std::transform(model->_mat->data()+off,model->_mat->data()+off+ax1.n,sig->_mat->data(),
														[](const std::complex<float> &v) {return v.real(); } );
						pad->forward(sig, sig_pad, 0);
						fft->forward(sig_pad,sig_pad,1);
						for (int i1=1; i1<ax1.n/2; ++i1) {
							(*sig_pad->_mat)[i1] *= -2.f;
						}
						for (int i1=ax1.n/2; i1<ax1.n; ++i1) {
							(*sig_pad->_mat)[i1] = 0.f;
						}
						fft->adjoint(sig_pad,sig_pad,1);
						pad->adjoint(sig, sig_pad, 0);
						std::transform(sig->_mat->data(),sig->_mat->data()+ax1.n,
														data->_mat->data()+off, data->_mat->data()+off,
														[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );
					}
				}
			// });
		// });
	};

	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) model->scale(0);

	};

private:
	int _npad_;
};

}
