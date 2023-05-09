#include <Average.h>

using namespace SEP;


void Average::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<float1DReg> data, bool add) {
	int nz = model->getHyper()->getAxis(2).n;
	int nx = model->getHyper()->getAxis(1).n;

	if (!add) data->scale(0);
	
	for (int j=0; j<nz; j++) {
		for (int i=0; i<nx; i++) {
			(*data->_mat)[j] += (*model->_mat)[j][i]/nx;
		}
	}
}
 

void Average::adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<float1DReg> data, bool add) {
	int n1 = model->getHyper()->getAxis(1).n;
	int n2 = model->getHyper()->getAxis(2).n;

	if (!add) model->scale(0);

	for (int j=0; j<n2; j++) {
		for (int i=0; i<n1; i++) {
			(*model->_mat)[j][i] += (*data->_mat)[i]/n2;
		}
	}
}