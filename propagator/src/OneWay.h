#pragma once

#include <Operator.h>
#include <float2DReg.h>
#include <complex2DReg.h>
#include <FFT1.h>
#include "paramObj.h"
#include <OneStep.h>
#include <Scatter.h>
#include <boost/range/irange.hpp>

namespace SEP {

class OneWay : public Operator<complex2DReg,complex2DReg>
{
public:
	OneWay(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par, const std::shared_ptr<RefSampler>& ref);
	virtual inline void setFreq(float freq) {prop->setFreq(freq);};

	virtual inline std::shared_ptr<OneStep>& getProp() {return prop;};
	virtual inline std::vector<int>& getBounds() {return bounds;}

	virtual void forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add) {};
	virtual void adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add) {};

	virtual ~OneWay() {};

	// inline void reverse() {zrange = boost::irange(zrange[0],zrange[nz-1]);}

protected:
	// boost::integer_range<int> zrange;
	std::shared_ptr<complex1DReg>  _wfld_next;
	std::shared_ptr<complex1DReg> _wfld_prev, _wfld_temp;

	typedef complex2D::array_view<1>::type complex2D_view;
	typedef boost::multi_array_types::index_range range;

	std::shared_ptr<FFT1> fft_k_prev;
	std::shared_ptr<OneStep> prop;
	std::shared_ptr<Scatter> src_filter;
	std::vector<int> bounds = std::vector<int>(2,0);

	int oz,orec;
	int nz,nx;

};

class Down : public OneWay
{
public:
	Down(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,const std::shared_ptr<RefSampler>& ref)
	: OneWay(slow,par,ref) {
		bounds[0] = 0;
		bounds[1] = nz-1;
		// zrange = boost::irange(0,nz-1);
	}

	void forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add);
	void adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add);

};

class DownSource : public Down
{
public:
	DownSource(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,const std::shared_ptr<RefSampler>& ref)
	: Down(slow,par,ref) {
		int ntaylor = par->getInt("ntaylor",1);
		src_filter = std::make_shared<Scatter> (slow,ntaylor,par);
	}

	inline void setFreq(float freq) {prop->setFreq(freq); src_filter->setFreq(freq);};
	void forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add);
	void adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add);

};

class Up : public OneWay
{
public:
	Up(const std::shared_ptr<complex2DReg>& slow, const std::shared_ptr<paramObj>& par, const std::shared_ptr<RefSampler>& ref)
	: OneWay(slow,par,ref) {
		bounds[0] = nz-1;
		bounds[1] = 0;
	};

	void forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add);
	void adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add);

};

class UpSource : public Up
{
public:
	UpSource(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,const std::shared_ptr<RefSampler>& ref)
	: Up(slow,par,ref) {
		int ntaylor = par->getInt("ntaylor",1);
		src_filter = std::make_shared<Scatter> (slow,ntaylor,par);
	}

	inline void setFreq(float freq) {prop->setFreq(freq); src_filter->setFreq(freq);};
	void forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add);
	void adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add);

};

}
