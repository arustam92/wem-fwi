#pragma once

#include <Operator.h>
#include <complex1DReg.h>

namespace SEP {

class Selector : public Operator<complex1DReg,complex1DReg>
{

public:

	Selector() {};

	void setLocation(std::vector<int>& loc) {_loc = loc;}

	void forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add);

private:
	std::vector<int> _loc;

};

}
