#pragma once
#include <OneWay.h>
#include <Operator.h>
#include <OneStep.h>
#include <Scatter.h>
#include <paramObj.h>
#include <Reflect.h>
#include <ic.h>


namespace SEP {

class LinOneWay : public Operator<complex2DReg,complex2DReg>
{
public:

	LinOneWay(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,
				const std::shared_ptr<complex2DReg>& bg_wfld, const std::shared_ptr<OneWay>& oneway);

	inline void setFreq(float freq){sc->setFreq(freq);prop->setFreq(freq);}
	

protected:
	const std::shared_ptr<complex2DReg>& _bg_wfld;
	std::shared_ptr<Scatter> sc;
	std::shared_ptr<OneStep> prop;
	std::shared_ptr<Transmission> trans;
	std::shared_ptr<IC> ic;

	const std::shared_ptr<complex2DReg>& _slow;

	int ntaylor;
	float w, _dz;

	std::shared_ptr<complex2DReg> _wfld_sc;
	std::shared_ptr<complex1DReg>  _wfld_next;
	std::shared_ptr<complex1DReg> _wfld_prev, _wfld_ref;

	int nz,nx;
	std::vector<int> bounds = std::vector<int> (2,0);
};

class LinDown : public LinOneWay
{
public:
	LinDown(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,
				const std::shared_ptr<complex2DReg>& bg_wfld, const std::shared_ptr<OneWay>& oneway)
	: LinOneWay(slow,par,bg_wfld,oneway)  {};

	void forward(const std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, const std::shared_ptr<complex2DReg> data, bool add);
};

class LinUp : public LinOneWay
{
public:
	LinUp(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,
				const std::shared_ptr<complex2DReg>& bg_wfld, const std::shared_ptr<OneWay>& oneway)
	: LinOneWay(slow,par,bg_wfld,oneway) {};

	void forward(const std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, const std::shared_ptr<complex2DReg> data, bool add);
};

}
