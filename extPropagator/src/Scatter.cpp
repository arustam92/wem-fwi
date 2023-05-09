#include <Scatter.h>
#include <cmath>

using namespace SEP;

Scatter::Scatter(std::shared_ptr<float2DReg> slow, int ntaylor) {

	_ntaylor = ntaylor;
	_slow = slow;
	_dz = slow->getHyper()->getAxis(2).d;
	nx = slow->getHyper()->getAxis(1).n;

	int nk = slow->getHyper()->getAxis(1).n;
	float dk = 2*pi/(slow->getHyper()->getAxis(1).d*(nk));
	k.resize(boost::extents[nk]);

	for (int ik=0; ik<k.size()/2; ik++) {
		k[ik] = ik*dk;
		k[nk-ik-1] = k[ik];
		// k[ik+nk/2] = (-nk/2 + ik)*dk;
	}

	// for reference velocities
	_wfld_ref = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
	_wfld_kx = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
	fft_out = std::make_shared<FFT1>(_wfld_ref,"out",1,1);
	fft_in = std::make_shared<FFT1>(_wfld_ref,"in",1,1);

}

void Scatter::forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if(!add) data->scale(0.);

	fft_out->forward(model,_wfld_kx,false);

	for (int it=0; it<=_ntaylor; it++) {
		scatter(_wfld_kx,_wfld_ref,it); //scaling by k2/w2
		fft_in->adjoint(_wfld_ref,_wfld_ref,1);
		slow_scale(_wfld_ref,_wfld_ref,_iz,it); //scaling by slownesses and coefs
		set_next_wfld(_wfld_ref,data,-1); //and multiply by -iwdz
	}
}


void Scatter::adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add) {
	if(!add) model->scale(0.);

	fft_out->forward(data,_wfld_kx,false);

	for (int it=0; it<=_ntaylor; it++) {
		scatter(_wfld_kx,_wfld_ref,it); //scaling by k2/w2
		fft_in->adjoint(_wfld_ref,_wfld_ref,1);
		slow_scale(_wfld_ref,_wfld_ref,_iz,it); //scaling by slownesses and coefs
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

void Scatter::slow_scale(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, int iz, int it) {

	for (int ix=0; ix<nx; ix++) {
		(*data->_mat)[ix] = coef[it]*(*model->_mat)[ix] *
							float((std::pow((*_slow->_mat)[iz][ix],-2*it)));

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
