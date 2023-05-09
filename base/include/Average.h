

namespace SEP {

class Average : public Operator
{
public:
	Average(std::shared_ptr<floatHyper> model, std::shared_ptr<floatHyper> data) : mod(model),dat(data) {};
	~Average();
	
	void forward() {
		int ndim = dat->getHyper()->getNdimG1();
		std::vector<int> n;
		for (int i=ndim; i>0; i--) n.push_back(dat->getHyper->getAxis(i).n);

		for (int i=0; i<n.size(), i++) {


		}

	};
private:
	std::shared_ptr<floatHyper> mod, dat;

};

}