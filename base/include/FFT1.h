#ifndef FFT1_H
#define FFT1_H

#include <complexHyper.h>
#include <floatHyper.h>
#include <Operator.h>
#include <fftw3.h>
#include <axis.h>

namespace SEP{

class FFT1 : public Operator {

public:
	FFT1(std::shared_ptr<complexHyper> model, int rank, int a);
	// void forward(std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data, bool add);
	// void adjoint(std::shared_ptr<floatHyper> model, std::shared_ptr<complexHyper> data, bool add);

	void forward (std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add);
	void adjoint (std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add);

private:
	int _axis, _rank;
	int idist, odist, istride, ostride, howmany, n1;
	bool _fwd_exec, _adj_exec;
	axis new_axis;
};
}

#endif