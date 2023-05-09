#include <Scatter.h>
#include <cmath>

using namespace SEP;

Scatter::Scatter(std::shared_ptr<complex2DReg> slow, int ntaylor, std::shared_ptr<paramObj> par) {

	_ntaylor = ntaylor;
	_slow = slow;
	_dz = slow->getHyper()->getAxis(2).d;
	nx = slow->getHyper()->getAxis(1).n;
	_eps = par->getFloat("eps",0.04);
	// for reference velocities
	int tap = par->getInt("tap",0);
	_wfld_ref = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
	pad = std::make_shared<Pad>(_wfld_ref->getHyper(), tap,tap,false);
	model_pad = std::make_shared<complex1DReg>(pad->getPaddedHyper());
	_wfld_ref_pad = std::make_shared<complex1DReg>(model_pad->getHyper()->getAxis(1));

	_wfld_kx = std::make_shared<complex1DReg>(model_pad->getHyper()->getAxis(1));
	fft_out = std::make_shared<FFT1>(model_pad,"out",1,1);
	fft_in = std::make_shared<FFT1>(model_pad,"in",1,1);

	int nk = model_pad->getHyper()->getAxis(1).n;
	float dk = 2*pi/(model_pad->getHyper()->getAxis(1).d*(nk-1));
	k.resize(boost::extents[nk]);

	for (int ik=0; ik<k.size()/2; ik++) 
		k[ik] = ik*dk;
	for (int ik=1; ik<k.size()/2; ik++) 
		k[nk-ik] = -k[ik];

}

void Scatter::forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if(!add) data->scale(0.);

	pad->forward(model,model_pad,0);
	fft_out->forward(model_pad,_wfld_kx,false);

	for (int it=0; it<=_ntaylor; it++) {
		scatter(_wfld_kx,_wfld_ref_pad,it); //scaling by k2/w2
		fft_in->adjoint(_wfld_ref_pad,_wfld_ref_pad,1);
		pad->adjoint(_wfld_ref, _wfld_ref_pad, 0);
		slow_scale_fwd(_wfld_ref,_wfld_ref,_iz,it); //scaling by slownesses and coefs
		set_next_wfld(_wfld_ref,data,-1); //and multiply by -iwdz
	}
}


void Scatter::adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add) {
	if(!add) model->scale(0.);

	pad->forward(data,model_pad,0);
	fft_out->forward(model_pad,_wfld_kx,false);

	for (int it=0; it<=_ntaylor; it++) {
		scatter(_wfld_kx,_wfld_ref_pad,it); //scaling by k2/w2
		fft_in->adjoint(_wfld_ref_pad,_wfld_ref_pad,1);
		pad->adjoint(_wfld_ref, _wfld_ref_pad, 0);
		slow_scale_adj(_wfld_ref,_wfld_ref,_iz,it); //scaling by slownesses and coefs
		set_next_wfld(_wfld_ref,model,1); //and multiply by -iwdz
	}

}

void Scatter::scatter(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, int it) {

	float pow;
	for (int ik=0; ik < k.size(); ik++) {
		pow = std::pow(k[ik]/w,2*it);
		(*data->_mat)[ik] = (*model->_mat)[ik]*pow;
	}
}

void Scatter::slow_scale_fwd(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, int iz, int it) {

	std::complex<float> c;
	std::complex<float> I = {0.f, 1.f};
	for (int ix=0; ix<nx; ix++) {
		c = coef[it] * std::pow((*_slow->_mat)[iz][ix]-_eps*(*_slow->_mat)[iz][ix].real(),-.5f-it);
		(*data->_mat)[ix] = c * (*model->_mat)[ix];
	}

}

void Scatter::slow_scale_adj(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, int iz, int it) {

	std::complex<float> c;
	std::complex<float> I = {0.f, 1.f};
	for (int ix=0; ix<nx; ix++) {
		c = coef[it] * std::pow((*_slow->_mat)[iz][ix]-_eps*(*_slow->_mat)[iz][ix].real(),-.5f-it);
		(*model->_mat)[ix] = std::conj(c) * (*data->_mat)[ix];
	}

}

void Scatter::set_next_wfld(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, int sign) {

	float re, im;
	float re2, im2;
	for (int ix=0; ix<nx; ix++) {

		re = -sign*w*_dz*(*model->_mat)[ix].imag();
		im = sign*w*_dz*(*model->_mat)[ix].real();
		re2 = (*data->_mat)[ix].real() + re;
		im2 = (*data->_mat)[ix].imag() + im;
		(*data->_mat)[ix] = {re2,im2};
	}
}
