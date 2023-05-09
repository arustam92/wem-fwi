#include <tbb/tbb.h>
#include "Phshift.h"
#include "phase_shift.h"
#include <cmath>

using namespace SEP;

Phshift::Phshift(float dz, std::shared_ptr<hypercube> hyp, int nref, float eps) {
	_dz = dz;
	_eps = eps;

	nk = hyp->getAxis(1).n;
	dk = 2*pi/(hyp->getAxis(1).d*(nk-1));

	k.resize(boost::extents[nk]);
	for (int i=0; i<nref; ++i) {
		kz_prev.push_back(std::shared_ptr<complex1D> (new complex1D(boost::extents[nk])));
		kz.push_back(std::shared_ptr<complex1D> (new complex1D(boost::extents[nk])));
		std::fill(kz_prev[i]->data(),kz_prev[i]->data()+nk,0);
	}
	for (int ik=0; ik<k.size()/2; ik++) 
		k[ik] = ik*dk;
	for (int ik=1; ik<k.size()/2; ik++) 
		k[nk-ik] = -k[ik];

	// lookup table setup
	// min = std::abs(k[nk/2]) * 5000 / (2*pi*0.5);
	// min = 1 - min * min;
	min = 0.f;
	max = 1.f;
	ntable = 10001;
	dtable = (max-min)/(ntable-1);

	_sqrt.resize(boost::extents[ntable]);
	_amp.reset(new float1D(boost::extents[ntable]));
	epsilon.resize(boost::extents[ntable]);

	float neg_percent = 1. - 1./(max-min);
	int i0 = -neg_percent*(ntable-1);
	_sqrt.reindex(i0);
	_amp->reindex(i0);


	lookupTable();
}

void Phshift::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	float* mre = &reinterpret_cast<float(&)[2]>(model->_mat->data()[0])[0];
	float* dre = &reinterpret_cast<float(&)[2]>(data->_mat->data()[0])[0];

	if (!add) data->scale(0.);

	float a, b, re, im, c, cos2, simag;
	std::complex<float> I = {0,1.f};
	simag = _s.imag();
	

	ispc::ps_fwd(k.size(),k.data(),w,_dz,_s.real(),simag,_eps,mre, dre);

		// for (int ik=0; ik < k.size(); ik ++) {

		// 	cos2 = std::max(0.f, 1.f - k[ik]*k[ik] / (w*w*_s.real()*_s.real()));
		// 	int ind = cos2 / dtable;
		// 	float eps = _eps;
		// 	// std::cerr << eps << '\n';

		// 	a = w*w*(_s.real()*_s.real() - simag*simag) - k[ik]*k[ik];

		// 	b = 2*w*w*_s.real()*(simag-eps*_s.real()/(w*w));
		// 	c = std::sqrt(a*a + b*b);
		// 	re = std::sqrt((c+a)/2);
		// 	im = -std::sqrt((c-a)/2);
		// 	// (*data->_mat)[ik] += f;
		// 	(*data->_mat)[ik] += (*model->_mat)[ik] * std::exp((-I*re+im)*_dz);
		// 	// (*model->_mat)[ik] *= R;

	 	// }
}

void Phshift::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
	
	float* mre = &reinterpret_cast<float(&)[2]>(model->_mat->data()[0])[0];
	float* dre = &reinterpret_cast<float(&)[2]>(data->_mat->data()[0])[0];
	
	if(!add) model->scale(0.);

	float a, b, re, im, c, cos2, simag;
	std::complex<float> I = {0,1.f};
	simag = _s.imag();
	
	ispc::ps_adj(k.size(),k.data(),w,_dz,_s.real(),simag,_eps,mre, dre);

		// for (int ik=0; ik < k.size(); ik ++) {

		// 	cos2 = std::max(0.f, 1.f - k[ik]*k[ik] / (w*w*_s.real()*_s.real()));
		// 	int ind = cos2 / dtable;
		// 	float eps = _eps;
		// 	// std::cerr << eps << '\n';

		// 	a = w*w*(_s.real()*_s.real() - simag*simag) - k[ik]*k[ik];
		// 	b = 2*w*w*_s.real()*(simag-eps*_s.real()/(w*w));
		// 	c = std::sqrt(a*a + b*b);
		// 	re = std::sqrt((c+a)/2);
		// 	im = -std::sqrt((c-a)/2);
		// 	// (*data->_mat)[ik] += f;
		// 	(*model->_mat)[ik] += (*data->_mat)[ik] * std::exp((I*re+im)*_dz);
		// 	// (*model->_mat)[ik] *= R;

	 	// }
 	// });

}


void Phshift::lookupTable() {
	// auto i0 = _sqrt.index_bases()[0];
	float pi = 4.*std::atan(1.);
	std::vector<float> angles(4,0);
	std::vector<int> cos(4,0);
	angles[0] = 10;
	angles[1] = 30;
	angles[2] = 60;
	angles[3] = 70;

	for (int i=0; i<4; ++i) {
		float a = angles[i] * pi/180; // in radians
		cos[i] = std::cos(a)*std::cos(a) / dtable;
	}
	for (int i=0; i < cos[3]; ++i)
		epsilon[i] = _eps;

	int tap1 = cos[2]-cos[3]-1;
	for (int i=cos[3]; i < cos[2]; ++i) {
		epsilon[i] = _eps * (std::cos(pi/tap1 * (i-cos[3])) + 1)/2;
	}

	int tap2 = cos[0] - cos[1] - 1;
	for (int i=cos[1]; i < cos[0]; ++i)
		epsilon[i] = _eps * (-std::cos(pi/tap2 * (i-cos[1])) + 1)/2;

	for (int i=cos[0]; i < ntable; ++i)
		epsilon[i] = _eps;

	// for (auto i=amax; i<n; i++) {
	// 	float val = i*dtable;
	// 	_sqrt[i] = {std::sqrt(val),0};
	// }
	//
	// int tap = amax - amin - 1;
	// for (auto i=amin; i < amax; ++i) {
	// 	float amp = _eps * (std::cos(pi/tap*(i-amin))+1)/2;
	// 	std::complex<float> val = {i*dtable, -amp};
	// 	_sqrt[i] = std::sqrt(val);
	// }
	// for (auto i=i0; i < amin; ++i) {
	// 	std::complex<float> val = {i*dtable, -_eps};
	// 	_sqrt[i] = std::sqrt(val);
	// }

}
