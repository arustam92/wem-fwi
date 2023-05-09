#include <Born_refl.h>

using namespace SEP;

Born_refl::Born_refl (std::shared_ptr<float2DReg> slow, std::shared_ptr<Injection> inj_rec,
	std::shared_ptr<Up> up, std::shared_ptr<complex2DReg> bg_wfld,std::shared_ptr<complex2DReg> sc_wfld) {

	drefl = std::make_shared<dReflect> (slow);
	_bg_wfld = bg_wfld;
	_up = up;
	_inj_rec = inj_rec;
	_sc_wfld = sc_wfld;

	dr = std::make_shared<float2DReg> (slow->getHyper());
	ic = std::make_shared<IC>(bg_wfld);


}

void Born_refl::forward (std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	drefl->forward(model,dr,0);

	ic->forward(dr,_sc_wfld,0);

	_up->forward(_sc_wfld,_sc_wfld,1);

	_inj_rec->adjoint(data,_sc_wfld,1);
}

void Born_refl::adjoint (std::shared_ptr<float2DReg> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) model->scale(0.);

	// float sc = 0.01;

	_inj_rec->forward(data,_sc_wfld,0);

	_up->adjoint(_sc_wfld,_sc_wfld,1);

	// cross correlation with bg wfld --> dr
	ic->adjoint(dr, _sc_wfld,false);

	// gradient mixing scheme
	// dr->scale(sc);

	// sum because of the imaging condition
	drefl->adjoint(model,dr,add);
}
