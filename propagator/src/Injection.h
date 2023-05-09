#pragma once

#include <complex1DReg.h>
#include <complex2DReg.h>
#include <complex3DReg.h>
#include <Operator.h>

/* TODO
 Include the case if gaussian is outside of the boundary
*/

namespace SEP {

class Injection :
public Operator<complex1DReg,complex2DReg>,
public Operator<complex2DReg,complex2DReg>,
public Operator<complex3DReg,complex2DReg> {
public:
	Injection(int iz, int ix, int ng = 0, int tap = 0) {
		_iz.push_back(iz);
		_ix.push_back(ix);
		_ng = ng;
		_tap = tap;
		gauss.reset(new boost::multi_array<float, 1>(boost::extents[2*ng+1]));
		compGauss();
	};
	Injection(const std::vector<int>& iz, const std::vector<int>& ix, int ng = 0, int tap = 0)  {
		_ix=ix;
		_iz=iz;
		_ng = ng;
		_tap = tap;
		gauss.reset(new boost::multi_array<float, 1>(boost::extents[2*ng+1]));
		compGauss();
	}

	inline void setStep(int step) {_istep=step;}
	inline void setShot(int shot) {_ishot=shot;}


	void forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex2DReg> data, bool add);

	void forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add);

	void forward(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex2DReg> data, bool add);
	void adjoint(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex2DReg> data, bool add);

	void adjoint(std::shared_ptr<complex4DReg> model, std::shared_ptr<complex2DReg> data, bool add);

private:

	void compGauss() {
		// its actually cosine squared
		float pi = 4*std::atan(1.);
		if (_ng==0) (*gauss)[0] = 1.;
		else {
			for (int i=0; i<gauss->size()/2; i++) {
				(*gauss)[i] = (-std::cos(pi/_ng*i)+1)/2;
			}
			for (int i=gauss->size()/2; i<gauss->size(); i++) {
				(*gauss)[i] = (std::cos(pi/_ng*(i-gauss->size()+_ng+1))+1)/2;
			}
		}
	}

	std::vector<int> _ix, _iz;
	int _istep, _ng, _ishot, _tap;
	std::shared_ptr<boost::multi_array<float, 1>> gauss;

};

}
