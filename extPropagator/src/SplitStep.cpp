#include <SplitStep.h>

using namespace SEP;

void SplitStep::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
	if (!add) data->scale(0.);
	
	std::complex<float> ss_corr;
	for (int i=0; i < _loc.size(); ++i) {
		ss_corr = {0., w*_dz*((*_slow->_mat)[_iz][_loc[i]] - _sref)};
		(*data->_mat)[_loc[i]] = (*model->_mat)[_loc[i]]*std::exp(-ss_corr);
	}
}

void SplitStep::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
	if (!add) model->scale(0.);

	std::complex<float> ss_corr;
	for (int i=0; i < _loc.size(); ++i) {
		ss_corr = {0., w*_dz*((*_slow->_mat)[_iz][_loc[i]] - _sref)};
		(*data->_mat)[_loc[i]] = (*model->_mat)[_loc[i]]*std::exp(ss_corr);
	}
}
