#include <iostream>
#include "float1DReg.h"
#include "float2DReg.h"
#include "ioModes.h"
#include <typeinfo>
#include <string>
#include <vector>

using namespace SEP;

int main(int argc, char **argv)

{

	/************************************* random stuff *********************************/	
	std::vector<float> test = {2.2,3.5,5.6,40.2};
	std::cout << "vector: " << test[1] << std::endl;
	float dz = 3.5;
	test.push_back(2.55/dz);
	std::cout << "vector: " << test[4] << std::endl;

	/************************************* MAIN IO **************************************/

	ioModes modes(argc,argv); // modes is an object
	std::shared_ptr <genericIO> ptrIo = modes.getDefaultIO(); // ptrIo is a (shared) pointer to a genericIO object
	
	/********************************* parameter file ***********************************/

	std::shared_ptr <paramObj> ptrPar = ptrIo->getParamObj(); // ptrPar is a (shared) pointer to a paramObj object1
	float dts = ptrPar->getFloat("dts",0.4);
	int nts = ptrPar->getInt("nts",413);
	float ots = ptrPar->getInt("ots",0.0);

	/************************************* wavelet **************************************/

	std::shared_ptr <genericRegFile> ptrWavelet = ptrIo->getRegFile(std::string("wavelet"),usageIn); 
	std::shared_ptr <hypercube> waveletHyper = ptrWavelet->getHyper();

	hypercube hyper(5);

	/* What I would like to do */ 
	//std::shared_ptr <float1DReg> wavelet; // Declare wavelet as a (shared) pointer to a float1DReg object
	//	wavelet = new float1DReg[1];

	/* Dynamic allocation for the wavelet */
	std::shared_ptr<float1DReg> wavelet(new float1DReg(waveletHyper)); 
	ptrWavelet->readFloatStream(wavelet); // Call the readFloatStream member function of the genericRegFile to read in the file

	// Wavelet is a shared pointer that points to a float1DReg object 
	// Wavelet has a member variable called "_mat" which is a pointer to a float1D array
	std::cout << "address of the array that _mat is pointing to: " << wavelet->_mat << std::endl; // Here, _mat is also a pointer 	
	std::cout << "data type of wavelet: " << typeid(wavelet).name() << std:: endl;
	std::cout << "data type of wavelet->_mat: " << typeid(wavelet->_mat).name() << std:: endl;	
	std::cout << "data type of *wavelet->_mat: " << typeid(*wavelet->_mat).name() << std:: endl; 		

	std::shared_ptr<float1DReg> waveletOut = wavelet->clone(); // copy wavelet

	for (int its = 0; its < nts; its++)
	{
		(*waveletOut->_mat)[its] = its;
	}

	/************************************* velocity model *******************************/

	std::shared_ptr <genericRegFile> ptrVel = ptrIo->getRegFile(std::string("vel"),usageIn); 
	std::shared_ptr <hypercube> velHyper = ptrVel->getHyper();	
	
	/* Test function of hypercube class */
	// getNdim
	std::cout << "nb of dim velHyper: " << velHyper->getNdim() << std::endl;
	std::cout << "nb of dim velHyper: " << velHyper->getNdimG1() << std::endl;	

	// Get Axes and display info
	std::vector<SEP::axis> velAxisVector = velHyper->getAxes();
	std::cout << "size of velAxisVector: " << velAxisVector.size() << std::endl;		

	for (int ivect = 0; ivect < velAxisVector.size(); ivect++)
	{
		std::cout << "ivect: " << ivect << std::endl;
		std::cout << "n: " << velAxisVector[ivect].n << std::endl;
		std::cout << "n: " << velHyper->getAxis(1).n << std::endl;		
		std::cout << "o: " << velHyper->getAxis(1).o << std::endl;				
		std::cout << "d: " << velHyper->getAxis(1).d << std::endl;					
		std::cout << "getN123: " << velHyper->getN123() << std::endl;						
		std::cout << "same size?: " << velHyper->sameSize(waveletHyper) << std::endl;							
	}
	
	// Allocate and read velocity model - there are 2 ways to do it
	std::shared_ptr<float2DReg> vel(new float2DReg(velHyper)); 
	std::shared_ptr<float2DReg> vel2 = std::make_shared<float2DReg>(velHyper); 	
	ptrVel->readFloatStream(vel); 
	ptrVel->readFloatStream(vel2); 	
	
	// Copy velocity model
	std::shared_ptr<float2DReg> velOut = vel->clone();

	// Testing functions of float2DReg class
	velOut->add(vel);
	velOut->scale(0.0);
	velOut->scaleAdd(3,vel,4); // velOut = 3*velOut + 4*vel
	velOut->random();

	std::cout << "dot: " << typeid((vel->_mat)[1][0]).name() << std::endl;
	std::cout << "absMax: " << velOut->absMax() << std::endl;
	std::cout << "min: " << velOut->min() << std::endl;
	std::cout << "max: " << velOut->max() << std::endl;	
	std::cout << "iz=0, ix=0: " << (*vel->_mat)[0][0] << std::endl;
	std::cout << "n1: " << velHyper->getAxis(1).n << std::endl;	
	std::cout << "n2: " << velHyper->getAxis(2).n << std::endl;	
 	std::cout << "norm type: " << vel->getNorm() << std::endl;
			
	// Apply something to the velocity model
	for (int ix = 0; ix < velHyper->getAxis(2).n; ix++)
	{
		for (int iz = 0; iz < velHyper->getAxis(1).n; iz++)
		{
			(*velOut->_mat)[ix][iz] = iz * ix;
		}
	}
	
	/************************************* write output *********************************/	

	std::shared_ptr <genericRegFile> ptrVelOut = ptrIo->getRegFile(std::string("velOut"),usageOut);	// Declare (the shared pointer to) the generic Reg file (out)
	ptrVelOut->setHyper(velHyper); 	// Copy the pointer to the hypercube object
	ptrVelOut->writeDescription(); // write description on SEP header file 
	ptrVelOut->writeFloatStream(velOut);

	return 0;
}
