#pragma once
#include <OneWay.h>
#include <Operator.h>
#include <OneStep.h>
#include <Scatter.h>
#include "LinOneWay.h"
#include "dReflect.h"
#include <paramObj.h>
#include <ic.h>


namespace SEP {

class LinTransDown : public LinOneWay
{
public:
	LinTransDown(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,
				const std::shared_ptr<complex2DReg>& bg_wfld, const std::shared_ptr<OneWay>& oneway)
	: LinOneWay(slow,par,bg_wfld,oneway)  {
		dtrans = std::make_shared<dTransmission>(slow);
		tmp_wfld = std::make_shared<complex2DReg>(slow->getHyper());
	};

	void forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add);
	void adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add);

private:
	std::shared_ptr<dTransmission> dtrans;
	std::shared_ptr<complex2DReg> tmp_wfld;
};

class LinTransUp : public LinOneWay
{
public:
	LinTransUp(const std::shared_ptr<complex2DReg>& slow,const std::shared_ptr<paramObj>& par,
				const std::shared_ptr<complex2DReg>& bg_wfld, const std::shared_ptr<OneWay>& oneway)
	: LinOneWay(slow,par,bg_wfld,oneway) {
		dtrans = std::make_shared<dTransmission>(slow);
	};

	void forward(const std::shared_ptr<complex2DReg>& model, std::shared_ptr<complex2DReg>& data, bool add);
	void adjoint(std::shared_ptr<complex2DReg>& model, const std::shared_ptr<complex2DReg>& data, bool add);

private:
	std::shared_ptr<dTransmission> dtrans;
};

}
