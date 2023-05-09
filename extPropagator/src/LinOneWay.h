#pragma once
#include <OneWay.h>
#include <Operator.h>
#include <OneStep.h>
#include <Scatter.h>
#include <paramObj.h>
#include <ic.h>


namespace SEP {

class LinOneWay : public Operator<float2DReg,complex2DReg>
{
public:

	LinOneWay(std::shared_ptr<float2DReg> slow,std::shared_ptr<paramObj> par,
				std::shared_ptr<complex2DReg> bg_wfld, std::shared_ptr<OneWay> oneway);

	inline void setFreq(float freq){sc->setFreq(freq);prop->setFreq(freq);}

protected:
	std::shared_ptr<complex2DReg> _bg_wfld;
	std::shared_ptr<Scatter> sc;
	std::shared_ptr<OneStep> prop;
	std::shared_ptr<IC> ic;

	std::shared_ptr<float2DReg> _slow;
	std::shared_ptr<float2DReg> temp;

	int ntaylor;
	float w, _dz;

	std::shared_ptr<complex2DReg> _wfld_sc;
	std::shared_ptr<complex1DReg>  _wfld_next;
	std::shared_ptr<float1DReg>  _image_z;
	std::shared_ptr<complex1DReg> _wfld_prev, _wfld_ref;

	int nz,nx;
	std::vector<int> bounds = std::vector<int> (2,0);
};

class LinDown : public LinOneWay
{
public:
	LinDown(std::shared_ptr<float2DReg> slow,std::shared_ptr<paramObj> par,
			std::shared_ptr<complex2DReg> bg_wfld, std::shared_ptr<OneWay> oneway)
	: LinOneWay(slow,par,bg_wfld,oneway)  {};

	void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
};

class LinUp : public LinOneWay
{
public:
	LinUp(std::shared_ptr<float2DReg> slow,std::shared_ptr<paramObj> par,
		std::shared_ptr<complex2DReg> bg_wfld, std::shared_ptr<OneWay> oneway)
	: LinOneWay(slow,par,bg_wfld,oneway) {};

	void forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
};

}
