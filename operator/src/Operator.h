#pragma once

#include <vector>
#include <floatHyper.h>
#include <complexHyper.h>

namespace SEP {

template <class M, class D>
class Operator
{
public:
	Operator() { };

	virtual void adjoint(std::shared_ptr<M> model, std::shared_ptr<D> data, bool add) {};
	virtual void forward(std::shared_ptr<M> model, std::shared_ptr<D> data, bool add) {};

	void setDomainRange(std::shared_ptr<M> model, std::shared_ptr<D> data) {
		_model = model->clone();
		_data = data->clone();
	}

	void dotTest() {
		std::shared_ptr<M> m1 = _model->clone();
		std::shared_ptr<D> d1 = _data->clone();
		_model->random();
		_data->random();

		forward(_model,d1,false);
		adjoint(m1,_data,false);

		// std::cerr << typeid(*this).name() << '\n';
		std::cerr << "********** ADD = FALSE **********" << '\n';
		std::cerr << "<m,A'd>: " << _model->dot(m1) << std::endl;
		std::cerr << "<Am,d>: " << std::conj(_data->dot(d1)) << std::endl;
		std::cerr << std::conj(_data->dot(d1))/_model->dot(m1) -1. << std::endl;
		std::cerr << "*********************************" << '\n';


		forward(_model,d1,true);
		adjoint(m1,_data,true);

		std::cerr << "********** ADD = TRUE **********" << '\n';
		std::cerr << "<m,A'd>: " << _model->dot(m1) << std::endl;
		std::cerr << "<Am,d>: " << std::conj(_data->dot(d1)) << std::endl;
		std::cerr << std::conj(_data->dot(d1))/_model->dot(m1) -1. << std::endl;
		std::cerr << "*********************************" << '\n';
};

protected:
	std::shared_ptr<M> _model;
	std::shared_ptr<D> _data;

};

}
