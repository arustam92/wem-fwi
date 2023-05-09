#include <ic.h>

using namespace SEP;

void IC::forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	if (!add) data->scale(0.);

	int nx = data->getHyper()->getAxis(1).n;
	int nz = data->getHyper()->getAxis(2).n;
	
	for (int iz=0; iz<nz; iz++) {
		for (int ix=0; ix<nx; ix++) {
			(*data->_mat)[iz][ix] += (*_bg_wfld->_mat)[iz][ix] * (*model->_mat)[iz][ix];
		}
	}
}
void IC::adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if (!add) model->scale(0.);
	int nx = model->getHyper()->getAxis(1).n;
	int nz = model->getHyper()->getAxis(2).n;

	for (int iz=0; iz<nz; iz++) {
		for (int ix=0; ix<nx; ix++) {
			(*model->_mat)[iz][ix] += std::conj((*_bg_wfld->_mat)[iz][ix]) * (*data->_mat)[iz][ix];
		}
}
}

void IC::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
	if (!add) data->scale(0.);

		for (int ix=0; ix<_nx; ix++) {
			data->_mat->data()[ix] += (*_bg_wfld->_mat)[_iz][ix] * model->_mat->data()[ix];
		}
}
void IC::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {

	if (!add) model->scale(0.);

	for (int ix=0; ix<_nx; ix++) {
		model->_mat->data()[ix] += std::conj((*_bg_wfld->_mat)[_iz][ix]) * data->_mat->data()[ix];
	}

}
