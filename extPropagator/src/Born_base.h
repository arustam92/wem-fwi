#pragma once
#include <WEM.h>

namespace SEP {

template <class M, class D>
class Born_base : public Operator<M,D>
{

public:
	Born_base() {};
	Born_base(std::shared_ptr<float1DReg> wave, std::shared_ptr<float2DReg> slow, std::shared_ptr<paramObj> par) {
		wem = std::make_shared<WEM> (slow,par);
		_sc_wfld = std::make_shared<complex2DReg> (wem->_wfld[0]->getHyper());
 	};

	void setThread(int ithread) {_ithread = ithread;}

	void setWfld(std::shared_ptr<complex3DReg> wavefield) {wem->setWfld(wavefield);}

	std::shared_ptr<WEM> wem;

protected:

	std::shared_ptr<complex2DReg> _sc_wfld;	
	
	int _ithread;

};



}