#include <WeightC.h>

using namespace SEP;

// void WeightC::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<float2DReg> data, bool add) {
// 	std::shared_ptr<float2DReg> temp = model->clone();
// 		if (!add) data->scale(0.);
// 		for (long long j=0; j<model->getHyper()->getAxis(2).n; j++) {
// 		  for (long long i=0; i<model->getHyper()->getAxis(1).n; i++) {
// 		  	(*data->_mat)[j][i] += (*temp->_mat)[j][i] * (*_w->_mat)[j][i];
// 		  }
// 		}
// }

// void WeightC::adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<float2DReg> data, bool add) {
// 	std::shared_ptr<float2DReg> temp = data->clone();
// 		if (!add) model->scale(0.);
// 		for (long long j=0; j<data->getHyper()->getAxis(2).n; j++) {
// 		  for (long long i=0; i<data->getHyper()->getAxis(1).n; i++) {		  
// 		  	(*model->_mat)[j][i] += (*temp->_mat)[j][i] * (*_w->_mat)[j][i];
// 		  }
// 		}
// }

void WeightC::forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	std::shared_ptr<complex2DReg> temp = model->clone();
		if (!add) data->scale(0.);
		for (long long j=0; j<model->getHyper()->getAxis(2).n; j++) {
		  for (long long i=0; i<model->getHyper()->getAxis(1).n; i++) {
		  	(*data->_mat)[j][i] += (*temp->_mat)[j][i] * (*_w->_mat)[i][j];
		  }
		}
}

void WeightC::adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	std::shared_ptr<complex2DReg> temp = data->clone();
		if (!add) model->scale(0.);
		for (long long j=0; j<data->getHyper()->getAxis(2).n; j++) {
		  for (long long i=0; i<data->getHyper()->getAxis(1).n; i++) {
		  	(*model->_mat)[j][i] += (*temp->_mat)[j][i] * (*_w->_mat)[i][j];
		  }
		}
}