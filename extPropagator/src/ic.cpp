#include <ic.h>

using namespace SEP;

void IC::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	if (!add) data->scale(0.);

	int nx = data->getHyper()->getAxis(1).n;
	int nz = data->getHyper()->getAxis(2).n;

	for (int iz=0; iz<nz; iz++) {
		for (int ix=0; ix<nx; ix++) {
			(*data->_mat)[iz][ix] += (*_bg_wfld->_mat)[iz][ix] * (*model->_mat)[iz][ix];
		}
	}
}
void IC::adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {

	if (!add) model->scale(0.);

	int nx = model->getHyper()->getAxis(1).n;
	int nz = model->getHyper()->getAxis(2).n;
	// float pi = std::atan(1)*4;

	for (int iz=0; iz<nz; iz++) {
		for (int ix=0; ix<nx; ix++) {
			(*model->_mat)[iz][ix] += ((*_bg_wfld->_mat)[iz][ix].real()*(*data->_mat)[iz][ix].real()+
									(*_bg_wfld->_mat)[iz][ix].imag()*(*data->_mat)[iz][ix].imag());
			// (*_illum->_mat)[iz][ix] += std::abs((*_bg_wfld)->_mat[iz][ix]);
		}
}
}
