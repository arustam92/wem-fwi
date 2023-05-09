#pragma once
#include <Acquisition.h>
#include <Born_refl.h>
#include <Born_up.h>
#include <Born_down.h>
#include <illumination.h>

namespace SEP {

class Born_full : public Operator<float2DReg,float3DReg>, Acquisition
{
public:
	// Born_full();

	Born_full(std::shared_ptr<float1DReg> wave, std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
					std::shared_ptr<float2DReg> scoord = nullptr, std::shared_ptr<float2DReg> rcoord = nullptr);
	// Born_full(std::shared_ptr<float3DReg> wave, std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
	// 				std::shared_ptr<float2DReg> scoord = nullptr, std::shared_ptr<float2DReg> rcoord = nullptr);
	void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add);
	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add) {
		(*this.*ptr_adjoint)(model,data,add);
	};

	Born_full(std::shared_ptr<complex1DReg> wave, std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
					std::shared_ptr<float2DReg> scoord = nullptr, std::shared_ptr<float2DReg> rcoord = nullptr);
	void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {
		(*this.*c_ptr_adjoint)(model,data,add);
	};

	void setBgSlow(std::shared_ptr<float2DReg> slow) {
		_slow = slow;
		ref.reset(new RefSampler(slow,nref));
		reflect.reset(new Reflect(slow));
		// ref = std::make_shared<RefSampler>(slow,nref);
		// reflect = std::make_shared<Reflect>(slow);
	}

	void setWfld(std::shared_ptr<complex3DReg> wavefield) {_full_wfld_w = wavefield;}
	std::shared_ptr<float3DReg> getWfld() {
		_full_wfld_t = std::make_shared<float3DReg> (_full_wfld_w->getHyper());
		FFT1 fft (_full_wfld_t,_full_wfld_w,1,1);
		fft.adjoint(_full_wfld_t,_full_wfld_w,0);
		return _full_wfld_t;
	}
	void snap_wfld(std::shared_ptr<complex2DReg> wfld,int i) {
		std::vector<axis> axes = _full_wfld_w->getHyper()->getAxes();
		for (size_t iz = 0; iz < axes[2].n; iz++) {
			for (size_t ix = 0; ix < axes[1].n; ix++) {
				(*_full_wfld_w->_mat)[iz][ix][i] += (*wfld->_mat)[iz][ix];
			}
		}
	}

private:
	void construct() {

		float df = 1./(wave_f->getHyper()->getAxis(1).d*wave_f->getHyper()->getAxis(1).n);
		float fMIN = param->getFloat("fmin",1);
		float fMAX = param->getFloat("fmax",0);
		if (fMAX == 0) {
			fmin = 1;
			fmax = param->getInt("nfreq",0);
			if (fmax==0) fmax = wave_f->getHyper()->getAxis(1).n/2+1;
		}
		else {
			// calculate indexes
			fmin = fMIN/df;
			fmax = fMAX/df;
		}

		freq.resize(boost::extents[fmax+1]);
		for (int i=0; i<freq.size(); i++) {
			freq[i] = (i*df);
		}

		if (param->getBool("illum",false)) {
			// illum = std::make_shared<Illumination>(_slow->getHyper());
			ptr_adjoint = &Born_full::adj_illum;
		}
		else {
			ptr_adjoint = &Born_full::adj_reg;
		}

		nshots = _src->getZ().size();
		std::cerr << "fmin = " << fMIN << "; fmax = " << fMAX << '\n';
		std::cerr << "nfreq = " << fmax-fmin << "; from imin = " << fmin << " to imax = " << fmax << '\n';
		std::cerr << "nshots = " << _src->getZ().size() << "; nrec = " << _rec->getZ().size() << '\n';
	}

protected:

	std::shared_ptr<complex3DReg> _full_wfld_w;
	std::shared_ptr<float3DReg> _full_wfld_t;

	std::shared_ptr<Reflect> reflect;
	std::shared_ptr<float2DReg> _slow;
	std::shared_ptr<RefSampler> ref;

	std::shared_ptr<paramObj> param;

	std::shared_ptr<complex1DReg> wave_f;

	int fmin, fmax;
	boost::multi_array<float,1> freq;

	int nref;
	int nshots;

	// std::shared_ptr<Illumination> illum;

	void (Born_full::*ptr_adjoint)(std::shared_ptr<float2DReg>, std::shared_ptr<float3DReg>, bool);
	void adj_reg(std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add);
	void adj_illum(std::shared_ptr<float2DReg> model, std::shared_ptr<float3DReg> data, bool add);

	void (Born_full::*c_ptr_adjoint)(std::shared_ptr<float2DReg>, std::shared_ptr<complex3DReg>, bool);
	void c_adj_reg(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);
	void c_adj_illum(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add);

};

}
