#include <Born_refl.h>
#include <algorithm>

using namespace SEP;

Born_refl::Born_refl (std::shared_ptr<Operator<complex2DReg,complex2DReg>> refl, std::shared_ptr<Injection>& inj_rec,
	std::shared_ptr<OneWay>& oneway, std::shared_ptr<complex2DReg>& bg_wfld,std::shared_ptr<complex2DReg>& sc_wfld) : 
	_bg_wfld(bg_wfld), _oneway(oneway), _inj_rec(inj_rec),_sc_wfld(sc_wfld) {

	std::shared_ptr<Reflect> R = std::static_pointer_cast<Reflect>(refl);
	drefl = std::make_shared<dReflect> (R->getBg());

	dr = std::make_shared<complex2DReg> (bg_wfld->getHyper());
	ic = std::make_shared<IC>(bg_wfld);


}

void Born_refl::forward (std::vector<std::shared_ptr<complex2DReg>> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) data->scale(0.);

	drefl->forward(model,dr,0);
	ic->forward(dr,_sc_wfld,0);
	_oneway->forward(_sc_wfld,_sc_wfld,1);
	_inj_rec->adjoint(data,_sc_wfld,1);
}

void Born_refl::adjoint (std::vector<std::shared_ptr<complex2DReg>> model, std::shared_ptr<complex3DReg> data, bool add) {

	if(!add) std::for_each(model.begin(), model.end(), [](std::shared_ptr<complex2DReg>& v) { v->scale(0); });
	
	_inj_rec->forward(data,_sc_wfld,0);
	_oneway->adjoint(_sc_wfld,_sc_wfld,1);
	// cross correlation with bg wfld --> dr
	ic->adjoint(dr, _sc_wfld, 0);
	drefl->adjoint(model,dr,1);
	// sum because of the imaging condition
	
}
