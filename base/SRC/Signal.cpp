#include <cmath>
#include "Signal.h"

using namespace SEP;

std::shared_ptr<Signal> Signal::xcorr(const std::shared_ptr <float2DReg> x) {

	std::vector<axis> a1 = x->getHyper()->getAxes();
	std::vector<axis> a2 = this->getHyper()->getAxes();
	assert(a1[1].n == a2[1].n);
	assert(a1[0].d == a2[0].d);

	axis a_out1 = a1[1];	// spatial axis is the same
	int shift = a1[0].n - (int((a2[0].o - a1[0].o)/a1[0].d));		// for acausal signals

	axis a_out0(a1[0].n + a2[0].n - 1, -(shift)*a1[0].d, a1[0].d);		// temporal axis
	std::shared_ptr <hypercube> hypOut (new hypercube(a_out0,a_out1));
	std::shared_ptr <Signal> out (new Signal(hypOut));

	//	Operator xcorr (model,data,adj,add);
	#pragma omp parallel for schedule (dynamic)
		for (int j = 0; j < a2[1].n; j++) {
			for (int i = 0; i < a2[0].n; i++) {
				float this_i = (*this->_mat)[j][i];
				for (int k=0; k < a1[0].n; k++) {
					//(*out->_mat)[j][k+i] += this_i * (*x->_mat)[j][a2[0].n-1-k];
					(*out->_mat)[j][k+i] += this_i * (*x->_mat)[j][a1[0].n-1-k];
				}
			}
		}

	return out;
	
}

std::shared_ptr<Signal> Signal::xcorr(const std::shared_ptr <float1DReg> x) {

		std::vector<axis> a1 = x->getHyper()->getAxes();
	std::vector<axis> a2 = this->getHyper()->getAxes();
	assert(a1[0].d == a2[0].d);

	axis a_out1 = a1[1];	// spatial axis is the same

	axis a_out0(a1[0].n + a2[0].n - 1, -(a2[0].n-1)*a2[0].d, a2[0].d);		// temporal axis
	std::shared_ptr <hypercube> hypOut (new hypercube(a_out0,a_out1));
	std::shared_ptr <Signal> out (new Signal(hypOut));

	#pragma omp parallel for
		for (int j = 0; j < a1[1].n; j++) {
			for (int i = 0; i < a1[0].n; i++) {
				float this_i = (*this->_mat)[j][i];
				for (int k=0; k < a2[0].n; k++) {
					(*out->_mat)[j][k+i] += this_i * (*x->_mat)[a2[0].n-1-k];
				}
			}
		}

	return out;
}

std::shared_ptr<Signal> Signal::conv(const std::shared_ptr <float2DReg> x) {

	std::vector<axis> a1 = this->getHyper()->getAxes();
	std::vector<axis> a2 = x->getHyper()->getAxes();
	assert(a1[1].n == a2[1].n);
	assert(a1[0].d == a2[0].d);

	axis a_out1 = a1[1];	// spatial axis is the same
	int shift = ((a1[0].o - a2[0].o)/a2[0].d);		// for acausal signals

	axis a_out0(a1[0].n + a2[0].n - 1, a2[0].o, a1[0].d);		// temporal axis
	std::shared_ptr <hypercube> hypOut (new hypercube(a_out0,a_out1));
	std::shared_ptr <Signal> out (new Signal(hypOut));

	#pragma omp parallel for schedule (dynamic)
		for (int j = 0; j < a2[1].n; j++) {
			for (int i = 0; i < a2[0].n; i++) {
				float x_i = (*x->_mat)[j][i];
				for (int k=0; k < a1[0].n; k++) {
					(*out->_mat)[j][k+i] += x_i * (*this->_mat)[j][k];
				}
			}
		}

	return out;
}


void Signal::levinson(const std::shared_ptr <float2DReg> acorr, const std::shared_ptr <float2DReg> xcorr, float eps) {

	

	int n1=(this->getHyper()->getAxis(1)).n;
	int n2=(this->getHyper()->getAxis(2)).n;
	int nlag = -(this->getHyper()->getAxis(1)).o/(this->getHyper()->getAxis(1)).d;

	int zeroA = ((acorr->getHyper()->getAxis(1)).n+1)/2-1;
	// s1+s2-1 = n;
	// s1 = n+1-s2;
	int zeroX = std::abs((xcorr->getHyper()->getAxis(1)).o)/(xcorr->getHyper()->getAxis(1)).d;
	// int zeroX = zeroA+1 + lenx-(zeroA+1);
	// zeroX-zeroA

	for (int j=0; j<n2; j++) {
		(*acorr->_mat)[j][zeroA] *= (1+eps);
	}

	#pragma omp parallel for schedule (dynamic)
	for (int j=0; j<n2; j++) {

		float e, v, c, ew, vw, cw;

			std::vector<float> acf(n1);
	std::vector<float> xcf(n1);
	std::vector<float> w(n1);
	std::vector<float> old_w(n1);
	std::vector<float> filt(n1);
	std::vector<float> old_filt(n1);



		for (int i=0; i<n1; i++) {
			acf[i] = (*acorr->_mat)[j][zeroA+i]/(*acorr->_mat)[j][zeroA];
			xcf[i] = (*xcorr->_mat)[j][zeroX-nlag+i]/(*acorr->_mat)[j][zeroA];
		}
		w[0] = 1;	
		v = 1;
		vw = 1;
		filt[0] = xcf[0];
		for (int iter=1; iter<n1; iter++) {
			old_w[iter] = 0;
			old_filt[iter] = 0;
			e=0;
			ew=0;
			for (int i=0; i<iter;i++){
				old_w[i] = w[i];
				old_filt[i] = filt[i];
			}
			for (int i=0; i<iter; i++) {
				ew += acf[i+1]*w[iter-i-1];
				e += acf[i+1]*filt[iter-i-1];
			}
			cw=ew/vw;
			vw = vw-cw*ew;
			c=(e-xcf[iter])/vw;			
			for (int i=1; i<=iter; i++) {
				w[i] = old_w[i] - cw*old_w[iter-i];	
			}
			for (int i=0; i<=iter; i++) {
				filt[i] = old_filt[i] - c*w[iter-i];			
			}
		}

		for(int i=0; i<n1; i++) {
			(*this->_mat)[j][i] = filt[i];
		}
	}
	return;
}

std::shared_ptr<float2DReg> Signal::kolmogorov(float eps, bool wantwave) {

	float val;
	
	std::vector<axis> ax = this->getHyper()->getAxes();
	std::vector<float> max(ax[1].n);
	double scale = 1.0/ax[0].n;

	std::shared_ptr<float2DReg> out = this->clone();
	std::shared_ptr <complex2DReg> s(new complex2DReg(this->getHyper()));
	std::shared_ptr <complex2DReg> spec = s->clone();
	std::shared_ptr<Operator> fft1 (new FFT1(out, s));
	fft1->forward();

	#pragma omp parallel for schedule (dynamic)
	for (int ix=0; ix<ax[1].n; ix++) {
		val = std::abs((*s->_mat)[ix][0]);
		for (int it=1; it<ax[0].n; it++) {
   			max[ix] = std::max(val, ((*s->_mat)[ix][it]).real());
		}
		for(int iw=0;iw<ax[0].n;iw++) {
			if (max[ix]<1e-12) break;
			(*spec->_mat)[ix][iw] = ((*s->_mat)[ix][iw]);
			(*s->_mat)[ix][iw] = std::log(std::abs((*s->_mat)[ix][iw])+eps*max[ix]);
		}

	}

	fft1->adjoint();
	out->scale(scale);

	#pragma omp parallel for schedule (dynamic)
	for (int ix=0; ix<ax[1].n; ix++) {
		for(int iw=1;iw<ax[0].n/2;iw++) {
			(*out->_mat)[ix][iw] *= 2;
		}
		for(int iw=ax[0].n/2;iw<ax[0].n;iw++) {
			(*out->_mat)[ix][iw] = 0;
		}

	}

	fft1->forward();

	if (!wantwave) {
		#pragma omp parallel for schedule (dynamic)
		for (int ix=0; ix<ax[1].n; ix++) {
			for(int iw=0;iw<ax[0].n;iw++) {
				(*s->_mat)[ix][iw] = (*spec->_mat)[ix][iw]*
									std::exp(s->_I*((*s->_mat)[ix][iw]).imag());
			}

		}
	}
	else {
		#pragma omp parallel for schedule (dynamic)
		for (int ix=0; ix<ax[1].n; ix++) {
			for(int iw=0;iw<ax[0].n;iw++) {
				(*s->_mat)[ix][iw] = std::abs((*spec->_mat)[ix][iw])*
									std::exp(s->_I*((*s->_mat)[ix][iw]).imag());
			}

		}
	}


	fft1->adjoint();
	out->scale(scale);

	return out;

}
	























