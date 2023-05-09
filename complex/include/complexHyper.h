#ifndef COMPLEX_HYPER_H
#define COMPLEX_HYPER_H

#pragma once
#include <hypercube.h>
#include <cassert>
#include <cstdint>
#include <sstream>
#include "Vector.h"
#include <complex>

namespace SEP {

class complexHyper : public Vector 
{

	private:
		std::shared_ptr<SEP::hypercube> _hyper;
 		std::complex<float> *_vals;
  		uint64_t _checkSum;
  		std::string norm;

	public:
  	
		/*********************************** Accessors **********************************/
  		bool getSpaceOnly() const { return _spaceOnly; }
		std::complex<float> *getVals() { return _vals; }
		const std::complex<float> *getCVals() const { return _vals; }
		std::string getNorm() const 
		{
			if (this->norm == std::string("")) return "L2";
			return this->norm;
		}		  

		std::shared_ptr<SEP::hypercube> getHyper() const { return _hyper; }  	
  	
		/*********************************** Constructors *******************************/
  		complexHyper() { ; }

		/*********************************** Modifiers **********************************/
  		void setHyper(std::shared_ptr<SEP::hypercube> h) { _hyper = h; }
  		void setData(std::complex<float> *ptr)
  		{
  			_vals = ptr;
    		setNotSpace();
  		}
  		
 		void setNorm(std::string nrm) 
		{
			if (nrm == std::string("L2"))
				{
		  			this->norm = nrm;
		  			return;
				} 
			else if (nrm == std::string("L1")) 
				{
		  			this->norm = nrm;
		  			return;
				}
			assert(1 == 2);
		} 		

		/******************************** Useful functions ******************************/
  		virtual void add(std::shared_ptr<complexHyper> vec);
  		virtual void scale(const std::complex<float> val);
  		virtual void scaleAdd(const std::complex<float> sc1, std::shared_ptr<complexHyper> vec2, const std::complex<float> sc2);
  		virtual void random();
  		//virtual void signum();
  		virtual double dot(std::shared_ptr<complexHyper> vec2) const;
  // 		float min() const;
		// float max() const;
 		virtual float absMax() const;		

  		void createMask(const float zero, const float err);

		double calcNorm() const
		{
			if (this->norm == std::string("L2"))
				return this->L2Obj();
			else 
			{
		  		return this->L1Obj();
			}
		}
  		void calcCheckSum();
  		void setCheckSum(const uint64_t x) { _checkSum = x; }
		bool isDifferent(std::shared_ptr<complexHyper> vec2) 
		{
			calcCheckSum();
			if (vec2->getCheckSum() != getCheckSum()) return true;
			return false;
		}

		double L2Obj() const;
		double L1Obj() const;

		// virtual void softClip(const float val);

		virtual void infoStream(const int lev, std::stringstream &x);

		virtual bool checkSame(const std::shared_ptr<SEP::complexHyper> vec2) const;
		uint64_t getCheckSum() { return _checkSum; }

		const std::complex<float> _I = std::complex<float>(0.0,1.0);
		

};

}  // namespace SEP

#endif
