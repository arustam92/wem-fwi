#pragma once
#include <float1DReg.h>
#include <float2DReg.h>
#include <complex1DReg.h>
#include <complex2DReg.h>
#include <Acquisition.h>
#include "ConformalMap.h"
#include "ShiftAndScale.h"
#include "ComplexLog.h"
#include "MapChain.h"
#include "Identity.h"
#include "Rotation.h"
#include "FFT1.h"
#include "Pad.h"
#include "StackAndSpread.h"
#include <fstream>


namespace SEP {

class PropParam : public Acquisition
{
public:
	PropParam(std::shared_ptr<float1DReg> wave, std::shared_ptr<hypercube> modspace, std::shared_ptr<paramObj> par) :
	Acquisition(modspace,par) {
		calcPropParam(wave,par);
		prepareConformalMaps(modspace,par);
		
		int pad = par->getInt("pad",0);
		pad_slow = std::make_shared<Pad>(modspace,pad,pad,true);
		pad_model = std::make_shared<Pad>(modspace,pad,pad,false);

		NxNz = pad_model->getPaddedHyper()->getAxis(1).n * pad_model->getPaddedHyper()->getAxis(2).n;
		stack = std::make_shared<Stack>();
	}

	PropParam(std::shared_ptr<complex1DReg> wave, std::shared_ptr<hypercube> modspace, std::shared_ptr<paramObj> par) :
	Acquisition(modspace,par) {
		throw SEPException("Constructor not implemented yet!");
	}

	int get_nfreq() {return freq.size();};

private:
	void calcPropParam(std::shared_ptr<float1DReg> wave, std::shared_ptr<paramObj> par) {

		param = par;

		Pad pad_t(wave->getHyper(),0,par->getInt("tpad",0),false);
		auto wave_pad = std::make_shared<float1DReg>(pad_t.getPaddedHyper());
		pad_t.forward(wave,wave_pad,0);
		_wave_f = std::make_shared<complex1DReg>(wave_pad->getHyper());
		FFT1 fft_wave (wave_pad,_wave_f, 1, 1);
		fft_wave.forward(wave_pad,_wave_f,false);

		nref = par->getInt("nref",1);

		float tmax = wave_pad->getHyper()->getAxis(1).o + (wave_pad->getHyper()->getAxis(1).d*(wave_pad->getHyper()->getAxis(1).n-1));
		float fMIN = par->getFloat("fmin",1);
		float fMAX = par->getFloat("fmax",0);
		int fmin, fmax;
		if (fMAX == 0) {
			fmin = 1;
			fmax = par->getInt("nfreq",0);
			if (fmax==0) fmax = wave_pad->getHyper()->getAxis(1).n/2+1;
		}
		else {
			// calculate indexes
			fmin = fMIN * tmax;
			fmax = fMAX * tmax + .5f;
		}

		freq.resize(boost::extents[fmax-fmin+1]);
		index.resize(boost::extents[fmax-fmin+1]);
		for (int i=0; i<freq.size(); i++) {
			freq[i] = fMIN + float(i)/tmax;
			index[i] = (i+fmin);
		}
		nshots = _src->getZ().size();

		std::cerr << "fmin = " << freq[0] << "; fmax = " << freq[freq.size()-1] << '\n';
		std::cerr << "nfreq = " << freq.size() << "; from imin = " << fmin << " to imax = " << fmax << '\n';
		std::cerr << "nshots = " << _src->getZ().size() << "; nrec = " << _rec->getZ().size() << '\n';

}


void prepareConformalMaps(std::shared_ptr<hypercube> modspace,std::shared_ptr<paramObj> par) {
	if (par->getInt("cmap",0) == 1) {
		std::complex<float> a = 0, b;
		axis ax1 = modspace->getAxis(1);
		axis ax2 = modspace->getAxis(2);
		float max1 = ax1.o + (ax1.n - 1)*ax1.d;
		float max2 = ax2.o + (ax2.n - 1)*ax2.d;

		for (int is=0; is < nshots; ++is) {
			float z = ax2.o + getSrc()->getZ(is)*ax2.d;
			float x = ax1.o + getSrc()->getX(is)*ax1.d;
			b = {z,x};
			for (int i2=0; i2 < ax2.n; i2 += ax2.n-1) {
				for (int i1=0; i1 < ax1.n; i1 += ax1.n-1) {
					std::complex<float> corner = {ax2.o + i2*ax2.d, ax1.o + i1*ax1.d};
					a = std::max(a.real(),std::abs(corner - b));
				}
			}
			a = 1.f/a;
			std::shared_ptr<ConformalMap> sas = std::make_shared<ShiftAndScale>(modspace,a,b);
			std::shared_ptr<ConformalMap> log = std::make_shared<ComplexLog> (sas->getOutHyper(),par->getFloat("eps_log",0.01));
	 		std::shared_ptr<ConformalMap> chain = std::make_shared<MapChain>(sas,log);
			cmap.push_back(chain);
		}
	}
	else if (par->getInt("cmap",0) == 2) {
		axis ax1 = modspace->getAxis(1);
		axis ax2 = modspace->getAxis(2);
		std::complex<float> a = {1.f, 0.f};
		std::complex<float> b = {ax2.o, ax1.o};
		std::shared_ptr<ConformalMap> sas = std::make_shared<ShiftAndScale>(modspace,a,b);
		// read angles
		Json::Value jsonArg;
		Json::Reader reader;
		std::ifstream file(par->getString("cmap_angles"));
		if (!reader.parse( file, jsonArg )) throw SEPException(reader.getFormattedErrorMessages());
		auto angles = jsonArg["a"];
		for (int is=0; is < nshots; ++is) {
			std::shared_ptr<ConformalMap> rot = std::make_shared<Rotation>(modspace,angles[is].asFloat());
			std::shared_ptr<ConformalMap> chain = std::make_shared<MapChain>(sas,rot);
			cmap.push_back(rot);
		}
		file.close();
	}
	else {
		auto id = std::make_shared<Identity>();
		for (int is=0; is < nshots; ++is) {
			cmap.push_back(id);
		}
	}
	}

protected:

	typedef boost::multi_array<int,1> int1D;
	std::vector<std::shared_ptr<ConformalMap>> cmap;

	std::shared_ptr<complex1DReg> _wave_f; //for nonlinear op
	std::shared_ptr<paramObj> param;
	float1D freq;
	int1D index;
	int nref;
	int nshots;
	int NxNz;
	std::shared_ptr<Pad> pad_model, pad_slow;
	std::shared_ptr<Stack> stack;

	inline std::shared_ptr<Device> getSrc() {return _src;};
	inline std::shared_ptr<Device> getRec() {return _rec;};
	inline std::shared_ptr<ConformalMap> getCmap(int &is) {return cmap[is];};
	inline std::shared_ptr<complex1DReg> getWavelet() {return _wave_f;};
	inline std::shared_ptr<paramObj> getParObj() {return param;};
	inline int& getNshots() {return nshots;};
	inline int& getNref() {return nref;};
	inline int& getNxNz() {return NxNz;};

	void sliceModel(std::shared_ptr<complex3DReg> input,std::shared_ptr<complex2DReg> output, int ifreq) {
		// auto view4d = (c4Dview)(*input->_mat)[indices[ifreq][range()][range()][range()]];
		// (*output->_mat) = view4d;
		int NxNz = getNxNz();
		std::copy_n(input->_mat->data()+ifreq*NxNz,NxNz,output->_mat->data());
	}
	void sliceModel(std::shared_ptr<complex2DReg> input,std::shared_ptr<complex3DReg> output, int ifreq) {
		// auto view4d = (c4Dview)(*output->_mat)[indices[ifreq][range()][range()][range()]];
		// view4d = (*input->_mat);
		int NxNz = getNxNz();
		// std::copy_n(input->_mat->data(),NxNz,output->_mat->data()+ifreq*NxNz);
		std::transform(input->_mat->data(),input->_mat->data()+NxNz,
										output->_mat->data()+ifreq*NxNz, output->_mat->data()+ifreq*NxNz,
										[](const std::complex<float> &i, const std::complex<float> &j) {return i+j; } );
	}

};

}
