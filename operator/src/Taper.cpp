#include <Taper.h>
#include <algorithm>
#include <execution>

using namespace SEP;

void Taper::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if(!add) data->scale(0.);

	int n = model->getHyper()->getAxis(1).n;

	float pi = 4*std::atan(1.);

	std::transform(
		std::execution::par_unseq,
		model->getVals(),model->getVals()+n,
								filter.data(),data->getVals(),
								[](const std::complex<float> &m, const float &f) {
									return m*f; 
									});
}

void Taper::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if(!add) model->scale(0.);

	int n = data->getHyper()->getAxis(1).n;

	float pi = 4*std::atan(1.);

	forward(data,model,add);
}
