#include <SplitStep.h>

using namespace SEP;

void SplitStep::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
	if (!add) data->scale(0.);

	std::complex<float> I = {0,1.f};
	std::complex<float> ss_true;
	std::complex<float> ss_corr;
	for (int i=0; i < _loc.size(); ++i) {
		ss_true = (*_slow->_mat)[_iz][_loc[i]];
		ss_corr = -I * w*_dz*(ss_true - _sref) / std::sqrt(_sref);
		(*data->_mat)[_loc[i]] = (*model->_mat)[_loc[i]]*std::exp(ss_corr);
	}
}

void SplitStep::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
	if (!add) model->scale(0.);

	std::complex<float> I = {0,1.f};
	std::complex<float> ss_true;
	std::complex<float> ss_corr;
	for (int i=0; i < _loc.size(); ++i) {
		ss_true = (*_slow->_mat)[_iz][_loc[i]];
		ss_corr = I * w*_dz*std::conj((ss_true - _sref) / std::sqrt(_sref));
		(*model->_mat)[_loc[i]] = (*data->_mat)[_loc[i]]*std::exp(ss_corr);
	}
}
