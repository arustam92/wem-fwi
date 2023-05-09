#pragma once

#include <Operator.h>
#include <float2DReg.h>
#include <complex2DReg.h>
#include <FFT1.h>
#include "paramObj.h"
#include <OneStep.h>

namespace SEP {

class OneWay : public Operator<complex2DReg,complex2DReg>
{
public:
	OneWay(const std::shared_ptr<float2DReg> slow,const std::shared_ptr<paramObj> par, std::shared_ptr<RefSampler> ref);
	inline void setFreq(float freq) {prop->setFreq(freq);};

	inline std::shared_ptr<OneStep> getProp() {return prop;};

	inline std::vector<int> getBounds() {return bounds;}

protected:
	std::shared_ptr<complex1DReg>  _wfld_next;
	std::shared_ptr<complex1DReg> _wfld_prev, _wfld_temp;

	std::shared_ptr<FFT1> fft_k_prev;
	std::shared_ptr<OneStep> prop;
	// std::shared_ptr<PSPI> pspi;
	std::vector<int> bounds = std::vector<int>(2,0);

	int oz,orec;
	int nz,nx;

};

class Down : public OneWay
{
public:
	Down(std::shared_ptr<float2DReg> slow,std::shared_ptr<paramObj> par,std::shared_ptr<RefSampler> ref)
	: OneWay(slow,par,ref) {
		bounds[0] = oz;
		bounds[1] = nz-1;
	}

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);

};

class Up : public OneWay
{
public:
	Up(std::shared_ptr<float2DReg> slow,std::shared_ptr<paramObj> par,std::shared_ptr<RefSampler> ref)
	: OneWay(slow,par,ref) {
		bounds[0] = nz-1;
		bounds[1] = orec;
	};

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);

};

}
