#include <tbb/tbb.h>
#include "Phshift.h"
#include <cmath>

using namespace SEP;

Phshift::Phshift(float dz, std::shared_ptr<hypercube> hyp, int nref, float eps) {
	_dz = dz;
	_eps = eps;

	nk = hyp->getAxis(1).n;
	dk = 2*pi/(hyp->getAxis(1).d*nk);
	xmin = hyp->getAxis(1).o;
	xmax = hyp->getAxis(1).o + (hyp->getAxis(1).n - 1)*hyp->getAxis(1).d;

	k.resize(boost::extents[nk]);
	sinc.resize(boost::extents[nk]);

	for (int i=0; i<nref; ++i) {
		kz_prev.push_back(std::shared_ptr<complex1D> (new complex1D(boost::extents[nk])));
		kz.push_back(std::shared_ptr<complex1D> (new complex1D(boost::extents[nk])));
		std::fill(kz_prev[i]->data(),kz_prev[i]->data()+nk,0);
	}
	for (int ik=1 - k.size()%2; ik<k.size()/2; ik++) {
		k[ik] = ik*dk;
		k[nk-ik - k.size()%2] = k[ik];
	}
	k[nk/2] = (nk/2)*dk;	// nyquist
	kN = k[nk/2];

	// lookup table setup
	min = std::abs(k[nk/2]) * 5000 / (2*pi*0.5);
	min = 1 - min * min;
	min = -1.f;
	max = 1.;
	ntable = 10001;
	dtable = (max-min)/(ntable-1);

	_sqrt.resize(boost::extents[ntable]);
	_amp.reset(new float1D(boost::extents[ntable]));

	float neg_percent = 1. - 1./(max-min);
	int i0 = -neg_percent*(ntable-1);
	_sqrt.reindex(i0);
	_amp->reindex(i0);


	lookupTable();
}

void Phshift::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if (!add) data->scale(0.);

	std::complex<float> T, R;
	float kz;
	float eps = 1e-5;
	std::complex<float> I = {0,1};
	std::complex<float> ph = {0,0},s={0,0},f={0,0};
	float re,im,amp,re2,im2;
	int ind;
	float sine;

		for (int ik=0; ik < k.size(); ik ++) {

			// kz = std::conj(std::sqrt((-I*w+eps)*(-I*w+eps)*_s*_s + k[ik]*k[ik]));
			// kz = std::conj(std::sqrt(w*w*(_s + I*eps)*(_s + I*eps) - k[ik]*k[ik]));

			// (*kz[_iref])[ik] = std::conj(std::sqrt(w*w*_s2 - k[ik]*k[ik] + I*eps));
			// R = ((*kz_prev[_iref])[ik] - (*kz[_iref])[ik])/((*kz_prev[_iref])[ik] + (*kz[_iref])[ik]);
			// T = float(2)*(*kz[_iref])[ik] / ((*kz_prev[_iref])[ik] + (*kz[_iref])[ik]);

			kz = std::max(min, 1 - k[ik]*k[ik]/(w*w*_s.real()));

			ind = kz/dtable;
			ph = _dz*w*std::sqrt(_s.real())*_sqrt[ind];
			// amp =  std::exp(ph.real());
			// re = amp * std::cos(ph.imag());
			// im = -amp * std::sin(ph.imag());
			// re2 = (*model->_mat)[ik].real();
			// im2 = (*model->_mat)[ik].imag();

			// f = {re*re2 - im*im2, re*im2 + im*re2};

			// (*data->_mat)[ik] += f;
			(*data->_mat)[ik] += (*model->_mat)[ik] * std::exp(-I*ph);
			// (*model->_mat)[ik] *= R;

	 	}
}

void Phshift::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if(!add) model->scale(0.);

		float kz;
		int ind;

		std::complex<float> ph;
		float amp,re,im;
		float re2,im2;
		float eps = 1e-5;
		std::complex<float> I = {0,1};

		for (int ik=0; ik < k.size(); ik ++) {

			// kz = (std::sqrt((-I*w+eps)*(-I*w+eps)*_s*_s + k[ik]*k[ik]));
			// kz = (std::sqrt(w*w*_s*_s - k[ik]*k[ik] + I*eps));
			// (*model->_mat)[ik] += std::exp(I*kz*_dz)*(*data->_mat)[ik];


			kz = std::max(0.f,1 - k[ik]*k[ik]/(w*w*_s.real()));
			//
			ind = kz/dtable;
			ph = _dz*w*std::sqrt(_s.real())*std::conj(_sqrt[ind]);
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
			(*model->_mat)[ik] += (*data->_mat)[ik] * std::exp(I*ph);

	 	}
 	// });

}


void Phshift::lookupTable() {
	auto i0 = _sqrt.index_bases()[0];
	float pi = 4.*std::atan(1.);
	float angle_max = 80;
	float angle_min = 90;
	angle_min = angle_min * pi/180; // in radians
	angle_max = angle_max * pi/180; // in radians
	int amin = std::cos(angle_min)*std::cos(angle_min) / dtable;
	int amax = std::cos(angle_max)*std::cos(angle_max) / dtable;
	int n = i0+ntable;

	for (auto i=amax; i<n; i++) {
		float val = i*dtable;
		_sqrt[i] = {std::sqrt(val),0};
	}

	int tap = amax - amin - 1;
	for (auto i=amin; i < amax; ++i) {
		float amp = _eps * (std::cos(pi/tap*(i-amin))+1)/2;
		std::complex<float> val = {i*dtable, -amp};
		_sqrt[i] = std::sqrt(val);
	}
	for (auto i=i0; i < amin; ++i) {
		std::complex<float> val = {i*dtable, -_eps};
		_sqrt[i] = std::sqrt(val);
	}

	float T = (xmax - xmin)/2;
	float S = (xmax + xmin)/2;
	std::complex<float> I = {0,1};
	float s;
	for (int i=0; i<nk; ++i) {
		if (k[i] == 0) {
			s = 1.f;
		}
		else {
			s = 2.f*std::sin(k[i]*T) / k[i];
		}
		sinc[i] = s * std::exp(-I*k[i]*S);
	}

}
