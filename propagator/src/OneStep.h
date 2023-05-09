#pragma once
#include <Operator.h>
#include <paramObj.h>
#include <float2DReg.h>
#include <complex1DReg.h>
#include <Phshift.h>
#include <Taper.h>
#include <Selector.h>
#include <FFT1.h>
#include <RefSampler.h>
#include "Pad.h"

namespace SEP{

class OneStep : public Operator<complex1DReg,complex1DReg>

{
public:
	OneStep (std::shared_ptr<complex2DReg> slow, std::shared_ptr<paramObj> par, std::shared_ptr<RefSampler> _ref) {

		_dz = slow->getHyper()->getAxis(2).d;

		nref = par->getInt("nref",1);
		ref = _ref;

		int tap = par->getInt("tap",0);
		int npad = par->getInt("pad",0);
		_wfld_ref = std::make_shared<complex1DReg>(slow->getHyper()->getAxis(1));
		pad = std::make_shared<Pad>(_wfld_ref->getHyper(), tap,tap,false);
		model_pad = std::make_shared<complex1DReg>(pad->getPaddedHyper());
		_wfld_ref_pad = std::make_shared<complex1DReg>(model_pad->getHyper()->getAxis(1));

		fft_in = std::make_shared<FFT1>(model_pad,"in",1,1);
		fft_out = std::make_shared<FFT1>(model_pad,"out",1,1);
		taper = std::make_shared<Taper>(_wfld_ref->getHyper()->getAxis(1).n,npad);
		ps = std::make_shared<Phshift>(_dz,model_pad->getHyper(),nref,par->getFloat("eps",0.04));
		select = std::make_shared<Selector>();
	}

	virtual void setFreq(float freq) {};
	virtual void setDepth(int iz) {};

	virtual void forward(const std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {};
	virtual void adjoint(std::shared_ptr<complex1DReg> model, const std::shared_ptr<complex1DReg> data, bool add) {};

	// void reflect(std::shared_ptr<complex1DReg> model,
	// 					std::shared_ptr<complex1DReg> data, bool add) {
	//
	// 	  if(!add) data->scale(0.);
	//
	// 	  fft_out->forward(model,model_kx,0);
	//
	// 		for (int iref=0; iref < nref; ++iref) {
	//
	// 			ps->setRef(iref);
	// 			ps->reflect(model_kx,_wfld_ref,0);
	//
	// 			fft_in->adjoint(_wfld_ref,_wfld_ref,1);
	// 			taper->forward(_wfld_ref,_wfld_ref,1);
	//
	// 			select->setLocation(ref->getRefLoc(_iz,iref));
	// 			select->forward(_wfld_ref,data,1);
	// 		}
	// 	}

		inline void reset() {ps->swapKz();}

		virtual ~OneStep() {};

protected:

	int nref, _iz;
	float _dz;

	std::shared_ptr<complex1DReg> _wfld_ref , model_pad, _wfld_ref_pad;

	std::shared_ptr<FFT1> fft_in ,fft_out;
	std::shared_ptr<Phshift> ps;
	std::shared_ptr<Pad> pad;
	std::shared_ptr<Taper> taper;
	std::shared_ptr<RefSampler> ref;
	std::shared_ptr<Selector> select;

};
}
