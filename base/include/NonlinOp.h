
namespace SEP {


class NonlinOp 
{
public:
	NonlinOp(model, data);
	~NonlinOp();

	virtual void forward();
	std::shared_ptr<Operator> linOp; // linearization
	
};
}

