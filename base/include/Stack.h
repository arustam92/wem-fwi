

namespace SEP {


class Stack : public Operator
{
public:
	Stack(std::shared_ptr<floatHyper> model, std::shared_ptr<floatHyper> data);
	Stack(std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data);

	~Stack();

	void forward();
	void adjoint();

	
	
};


}