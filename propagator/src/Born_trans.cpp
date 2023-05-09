#include <Born_trans.h>

using namespace SEP;

Born_trans::Born_trans (std::shared_ptr<complex2DReg>& slow, std::shared_ptr<Injection>& inj_rec,
	std::shared_ptr<Down>& down, std::shared_ptr<complex2DReg>& bg_wfld,std::shared_ptr<complex2DReg>& sc_wfld) : _bg_wfld(bg_wfld), _down(down), _inj_rec(inj_rec),_sc_wfld(sc_wfld) {

	dtrans = std::make_shared<dTransmission> (slow);
	dt = std::make_shared<complex2DReg> (bg_wfld->getHyper());
	ic = std::make_shared<IC>(bg_wfld);
}

void Born_trans::forward (std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	ic->forward(model,dt,0);
	dtrans->forward(dt,_sc_wfld,0);
	_down->forward(_sc_wfld,_sc_wfld,1);
	_inj_rec->adjoint(data,_sc_wfld,1);
}

void Born_trans::adjoint (std::shared_ptr<complex2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) model->scale(0.);

	_inj_rec->forward(data,_sc_wfld,0);
	_down->adjoint(_sc_wfld,_sc_wfld,1);
	dtrans->adjoint(dt, _sc_wfld, false);
	ic->adjoint(model, dt, 1);
}
