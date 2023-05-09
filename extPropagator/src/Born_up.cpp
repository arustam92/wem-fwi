#include <Born_up.h>

using namespace SEP;

Born_up::Born_up(std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
	std::shared_ptr<Injection> inj_rec, std::shared_ptr<Up> up,
	std::shared_ptr<complex2DReg> bg_wfld, std::shared_ptr<complex2DReg> sc_wfld) {

	_up = up;
	_inj_rec = inj_rec;

	_sc_wfld = sc_wfld;
	lin_up = std::make_shared<LinUp>(slow,par,bg_wfld,up);

}

void Born_up::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	lin_up->forward(model,_sc_wfld,0);

	_up->forward(_sc_wfld,_sc_wfld,1);

	_inj_rec->adjoint(data,_sc_wfld,1);

}

void Born_up::adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) model->scale(0.);

	_inj_rec->forward(data,_sc_wfld,0);

	_up->adjoint(_sc_wfld,_sc_wfld,1);

	lin_up->adjoint(model,_sc_wfld,1);

}
