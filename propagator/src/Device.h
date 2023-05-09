#pragma once
#include "paramObj.h"
#include <float2DReg.h>
#include <hypercube.h>
#include "ioModes.h"
#include <fstream>

namespace SEP {

class Device
{
public:
	Device(std::shared_ptr<paramObj> par, std::shared_ptr<hypercube> hypVel, const std::string& name) {
		if (name != "s" && name != "r") throw(SEPException("Device must be specified to be r or s!"));
		readModelParams(hypVel);
		if (!par->getString("coord","").empty()) createIrregGeom(par,name);
		else createRegGeom(par,name);
	}

	int& getX(int is) {return _ix[is];}
	int& getZ(int is) {return _iz[is];}

	std::vector<int>& getX() {return _ix;}
	std::vector<int>& getZ() {return _iz;}

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

	void createRegGeom(std::shared_ptr<paramObj> par, const std::string& name) {
		n = par->getInt("n"+name,0);
		if (!n) throw(SEPException("Number of devices " + name + " is zero!"));
		beg_x = par->getFloat("o"+name+"x",0.);
		if (beg_x < ox) throw(SEPException("Devices " + name + " start outside of the model!"));
		beg_z = par->getFloat("o"+name+"z",0.);
		if (beg_z < oz) throw(SEPException("Devices " + name + " start outside of the model!"));
		dx = par->getFloat("d"+name+"x",0.);
		dz = par->getFloat("d"+name+"z",0.);
		if (!dx && !dz) throw(SEPException("Devices " + name + " have spacing equal to zero!"));

		int indx, indz;
		for (int i=0; i<n; i++) {
			float x = beg_x + i*dx;
			float z = beg_z + i*dz;
			indx = (x - ox)/dx_vel + .5f;
			indz = (z - oz)/dz_vel + .5f;
			if (indx < 0 || indz < 0)
			throw(SEPException(std::string("Device ") +
	                       std::to_string(z) + std::string(",") + std::to_string(x) +
											 	std::string("} is out of bounds!")));
			_ix.push_back(indx);
			_iz.push_back(indz);
		}
	}

	void createIrregGeom(std::shared_ptr<paramObj> par, const std::string& name) {
    Json::Value jsonArg;   // will contains the root value after parsing.
    Json::Reader reader;
    std::ifstream file(par->getString("coord"));
    if (!reader.parse( file, jsonArg )) throw SEPException(reader.getFormattedErrorMessages());
		auto arrX = jsonArg[name+"x"];
		auto arrZ = jsonArg[name+"z"];
		if (arrX.size() != arrZ.size()) throw SEPException("z- and x-coordinates not equal size!");
		for (auto i=0; i < arrX.size(); ++i) {
			int indx = (arrX[i].asFloat() - ox)/dx_vel + .5f;
			int indz = (arrZ[i].asFloat() - oz)/dz_vel + .5f;
			if (indx < 0 || indz < 0)
				throw(SEPException(std::string("Device ") +
	                       std::to_string(arrZ[i].asFloat()) +
												 std::string(",") + std::to_string(arrX[i].asFloat()) +
											 	std::string("} is out of bounds!")));
			_ix.push_back(indx);
			_iz.push_back(indz);
		}
		file.close();
	}

};

}
