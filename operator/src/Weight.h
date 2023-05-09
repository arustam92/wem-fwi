#pragma once
#include <float2DReg.h>
#include <complex2DReg.h>
#include <Operator.h>

namespace SEP {

template <class T>
class Weight : public Operator<T,T>
{

public:
	Weight() {}
	Weight(std::shared_ptr<T> weight) {_w = weight;}

	void forward(std::shared_ptr<T> model, std::shared_ptr<T> data, bool add) {
			if (!add) data->scale(0.);
			int len = model->getHyper()->getN123();
			for (int i=0; i < len; i++)	data->getVals()[i] = model->getVals()[i] * _w->getVals()[i];
	};

	void adjoint(std::shared_ptr<T> model, std::shared_ptr<T> data, bool add) {
		if (!add) model->scale(0.);
		int len = model->getHyper()->getN123();
		for (int i=0; i < len; i++)	model->getVals()[i] = data->getVals()[i] * std::conj(_w->getVals()[i]);
	};

private:
	std::shared_ptr<T> _w;

};

}
