#include <Selector.h>

using namespace SEP;

void Selector::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
		if(!add) data->scale(0.);
		for (int i=0; i<_loc.size(); ++i) {
			(*data->_mat)[_loc[i]] = (*model->_mat)[_loc[i]];
		}
	};


void Selector::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex1DReg> data, bool add) {
		if(!add) model->scale(0.);
		for (int i=0; i<_loc.size(); ++i) {

			(*model->_mat)[_loc[i]] = (*data->_mat)[_loc[i]];
		}
};
