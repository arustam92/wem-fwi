#include <vector>
#include <floatHyper.h>
#include <complexHyper.h>

namespace SEP {

class Operator 
{
public:
	Operator() { };

	virtual void adjoint(std::shared_ptr<floatHyper> model, std::shared_ptr<floatHyper> data, bool add) {};
	virtual void forward(std::shared_ptr<floatHyper> model, std::shared_ptr<floatHyper> data, bool add) {};

	virtual void adjoint(std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data, bool add) {};
	virtual void forward(std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data, bool add) {};
	
	virtual void forward(std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add) {};
	virtual void adjoint(std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add) {};


	void getModel();

	std::vector<double> dotTest();

};

}
