#pragma once
#include <Device.h>
#include <hypercube.h>

namespace SEP {

class Acquisition
{
public:

	Acquisition() {} ;

	Acquisition(std::shared_ptr<hypercube> hypVel,std::shared_ptr<paramObj> par) {

			// std::cerr << "Creating irregular geometry by binary files..." << '\n';
			// std::cerr << "Receivers...:" << '\n';
			_rec = std::make_shared<Device>(par,hypVel,"r");
			// std::cerr << "Sources..." << '\n';
			_src = std::make_shared<Device>(par,hypVel,"s");

	}
	// Acquisition(std::shared_ptr<float2DReg> scoord, std::shared_ptr<float2DReg> rcoord, std::shared_ptr<hypercube> hypVel) {
	// 	_rec = std::make_shared<Device>(rcoord,hypVel);
	// 	_src = std::make_shared<Device>(scoord,hypVel);
	// }

	std::shared_ptr<Device> getRec() {return _rec;}
	std::shared_ptr<Device> getSrc() {return _src;}


protected:
	std::shared_ptr<Device> _rec, _src;

};

}
