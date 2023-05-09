#include <tbb/tbb.h>
#include "Phshift.h"
#include <cmath>

using namespace SEP;

Phshift::Phshift(float dz, std::shared_ptr<hypercube> hyp) {
	_dz = dz;

	int nk = hyp->getAxis(1).n;
	dk = 2*pi/(hyp->getAxis(1).d*nk);

	// lookup table setup
	std_dev = 1.;
	min = -std_dev;
	max = 1.;
	ntable = 10001;
	dtable = (max-min)/(ntable-1);

	_sqrt.reset(new complex1D(boost::extents[ntable]));
	_amp.reset(new float1D(boost::extents[ntable]));

	float neg_percent = 1. - 1./(max-min);
	int i0 = -neg_percent*(ntable-1);
	_sqrt->reindex(i0);
	_amp->reindex(i0);

	k.resize(boost::extents[nk]);
	kz_prev.reset(new complex1D(boost::extents[nk]));
	kz.reset(new complex1D(boost::extents[nk]));
	std::fill(kz_prev->data(),kz_prev->data()+nk,0);
	for (int ik=1 - k.size()%2; ik<k.size()/2; ik++) {
		k[ik] = ik*dk;
		k[nk-ik - k.size()%2] = k[ik];
	}
	k[nk/2] = (nk/2)*dk;	// nyquist

	lookupTable();
}

void Phshift::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if (!add) data->scale(0.);

	std::complex<float> T, R;
	// float kz;
	float eps = 1e-5;
	std::complex<float> I = {0,1};
	std::complex<float> ph = {0,0},s={0,0},f={0,0};
	float re,im,amp,re2,im2;
	int ind;
	float sine;

		for (int ik=0; ik < k.size(); ik++) {

			// kz = std::conj(std::sqrt((-I*w+eps)*(-I*w+eps)*_s*_s + k[ik]*k[ik]));
			// kz = std::conj(std::sqrt(w*w*(_s + I*eps)*(_s + I*eps) - k[ik]*k[ik]));

			(*kz)[ik] = std::conj(std::sqrt(w*w*(_s)*(_s) - k[ik]*k[ik] + I*eps));
			R = ((*kz_prev)[ik] - (*kz)[ik])/((*kz_prev)[ik] + (*kz)[ik]);
			T = float(2)*(*kz)[ik] / ((*kz_prev)[ik] + (*kz)[ik]);

			// (*data->_mat)[ik] = std::exp(-kz*_dz)*(*model->_mat)[ik];

			// kz = std::max(min,1 - k[ik]*k[ik]/(w*w*_s*_s));
			//
			// ind = kz/dtable;
			// ph = _dz*w*_s*(*_sqrt)[ind];
			// amp = std::exp(ph.real());
			// re = amp*std::cos(ph.imag());
			// im = -amp*std::sin(ph.imag());
			// re2 = (*_amp)[ind]*(*model->_mat)[ik].real();
			// im2 = (*_amp)[ind]*(*model->_mat)[ik].imag();

			// f = {re*re2 - im*im2, re*im2 + im*re2};
			// float j = std::abs(k[ik])/k.size();
			// s = float(1.*j/(j+1))*s + float(1./(j+1)) * f;


			// (*data->_mat)[ik] = f;
			(*data->_mat)[ik] = (*model->_mat)[ik] * std::exp(-I*_dz*(*kz)[ik]);

	 	}
		// kz_prev = kz;
}

void Phshift::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if(!add) model->scale(0.);

		std::complex<float> kz;
		// float kz;
		int ind;

		std::complex<float> ph;
		float amp,re,im;
		float re2,im2;
		float eps = 1e-5;
		std::complex<float> I = {0,1};

		for (int ik=0; ik < k.size(); ik++) {

			// kz = (std::sqrt((-I*w+eps)*(-I*w+eps)*_s*_s + k[ik]*k[ik]));
			kz = (std::sqrt(w*w*(_s)*(_s) - k[ik]*k[ik] + I*eps));
			(*model->_mat)[ik] += std::exp(I*kz*_dz)*(*data->_mat)[ik];


			// kz = std::max(min,1 - k[ik]*k[ik]/(w*w*_s*_s));
			//
			// ind = kz/dtable;
			// ph = _dz*w*_s*(*_sqrt)[ind];
			// amp = std::exp(ph.real());
			// re = amp*std::cos(ph.imag());
			// im = amp*std::sin(ph.imag());
			// re2 = (*model->_mat)[ik].real();
			// im2 = (*model->_mat)[ik].imag();
			//
			// re2 += re*(*data->_mat)[ik].real() - im*(*data->_mat)[ik].imag();
			// im2 += re*(*data->_mat)[ik].imag() + im*(*data->_mat)[ik].real();
			//
			// (*model->_mat)[ik] = {re2,im2};

	 	}
 	// });

}


void Phshift::lookupTable() {
	auto i0 = _sqrt->index_bases()[0];
	float pi = 4.*std::atan(1.);
	float angle = 80;
	angle = angle*pi/180; // in radians
	int cutoff = std::cos(angle)*std::cos(angle) / dtable;
	int n = i0+ntable;

	std::complex<float> x;
	// fill the negative square root for evanescent energy
	for (auto i=i0; i<0; i++) {
		float val = -i*dtable;
		// x = {val,0};
		// (*_sqrt)[i] = -std::sqrt(val);
		(*_sqrt)[i] = {-std::sqrt(val),0};
		(*_amp)[i] = 1;
	}
	// fill the positive square root for propagating energy
	for (auto i=0; i<n; i++) {
		float val = i*dtable;
		(*_sqrt)[i] = {0,std::sqrt(val)};
		(*_amp)[i] = 1;
	}

	int tap = cutoff;
	for (auto i=0; i<cutoff; ++i) {
		float val = i*dtable;
		(*_amp)[i] = (-std::cos(pi/tap*i)+1)/2;
		// (*_amp)[i] = 1;
		// (*_sqrt)[i] = {eps*(-std::cos(pi/cutoff*i)-1)/2,std::sqrt(val)};
	}


}
