#pragma once
#include "paramObj.h"
#include <float2DReg.h>
#include <hypercube.h>

namespace SEP {

class Device
{
public:
	Device(std::shared_ptr<paramObj> par, std::shared_ptr<hypercube> hypVel, std::string name) {

		readModelParams(hypVel);
		createRegGeom(par,name);

	}

	Device(std::shared_ptr<float2DReg> coord, std::shared_ptr<hypercube> hypVel) {

		readModelParams(hypVel);
		createIrregGeom(coord);

	}

	int getX(int is) {return _ix[is];}
	int getZ(int is) {return _iz[is];}

	std::vector<int> getX() {return _ix;}
	std::vector<int> getZ() {return _iz;}

private:
	std::vector<int> _ix, _iz;
	int n;
	float beg_x, beg_z,dx,dz;
	float dx_vel,dz_vel, ox, oz;

	void readModelParams(std::shared_ptr<hypercube> hyper) {
		dx_vel = hyper->getAxis(1).d;
		dz_vel = hyper->getAxis(2).d;
		ox = hyper->getAxis(1).o;
		oz = hyper->getAxis(2).o;
	}

	void createRegGeom(std::shared_ptr<paramObj> par, std::string name) {

		if (name=="src") {
				n = par->getInt("ns",1);
				beg_x = par->getFloat("osx",0.);
				beg_z = par->getFloat("osz",0.);
				dx = par->getFloat("dsx",0.);
				dz = par->getFloat("dsz",0.);
			}
		if (name=="rec") {
				n = par->getInt("nr",1);
				beg_x = par->getFloat("orx",0.);
				beg_z = par->getFloat("orz",0.);
				dx = par->getFloat("drx",0.);
				dz = par->getFloat("drz",0.);
		}

		int indx, indz;
		for (int i=0; i<n; i++) {
			indx = (beg_x + i*dx - ox)/dx_vel;
			indz = (beg_z + i*dz - oz)/dz_vel;
			if (indx < 0 || indz < 0) std::cerr << "Device " << i << " is out of bounds!" << '\n';;
			_ix.push_back(indx);
			_iz.push_back(indz);
		}
	}

	void createIrregGeom(std::shared_ptr<float2DReg> coord) {
		int indx, indz;
		for (int i=0; i < coord->getHyper()->getAxis(1).n; ++i) {
			indx = ((*coord->_mat)[0][i] - ox)/dx_vel;
			indz = ((*coord->_mat)[1][i] - oz)/dz_vel;
			_ix.push_back(indx);
			_iz.push_back(indz);
			if (indx < 0 || indz < 0) std::cerr << "Device " << i << " is out of bounds!" << '\n';;
		}
	}

};

}
