#include <Born_down.h>

using namespace SEP;

Born_down::Born_down(std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par,
	std::shared_ptr<Injection> inj_rec,
	std::shared_ptr<Down> down, std::shared_ptr<Up> up, std::shared_ptr<Reflect> reflect,
	std::shared_ptr<complex2DReg> bg_wfld, std::shared_ptr<complex2DReg> sc_wfld) {

	_down = down;
	_up = up;
	_reflect = reflect;
	_inj_rec = inj_rec;

	_sc_wfld = sc_wfld;
	lin_down = std::make_shared<LinDown>(slow,par,bg_wfld,down);

}

void Born_down::forward(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	lin_down->forward(model,_sc_wfld,0);

	_down->forward(_sc_wfld,_sc_wfld,1);

	_reflect->forward(_sc_wfld,_sc_wfld,1);

	_up->forward(_sc_wfld,_sc_wfld,1);

	_inj_rec->adjoint(data,_sc_wfld,1);

}

void Born_down::adjoint(std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) model->scale(0.);

	_inj_rec->forward(data,_sc_wfld,0);

	_up->adjoint(_sc_wfld,_sc_wfld,1);

	_reflect->adjoint(_sc_wfld,_sc_wfld,1);

	_down->adjoint(_sc_wfld,_sc_wfld,1);

	lin_down->adjoint(model,_sc_wfld,1);

}
