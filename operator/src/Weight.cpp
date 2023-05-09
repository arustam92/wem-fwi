#include <Weight.h>

using namespace SEP;

void Weight<T>::forward(std::shared_ptr<T> model, std::shared_ptr<T> data, bool add) {
		if (!add) data->scale(0.);
		int len = model->getHyper()->getN123();
		for (int i=0; i < len; i++)	data->getVals()[i] += model->getVals()[i] * _w->getVals()[i];
}

void Weight<T>::adjoint(std::shared_ptr<T> model, std::shared_ptr<T> data, bool add) {
	if (!add) model->scale(0.);
	int len = model->getHyper()->getN123();
	for (int i=0; i < len; i++)	model->getVals()[i] += data->getVals()[i] * _w->getVals()[i];
}
