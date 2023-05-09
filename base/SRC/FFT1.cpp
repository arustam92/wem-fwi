#include <FFT1.h>

using namespace SEP;

FFT1::FFT1 (std::shared_ptr<complexHyper> model, int rank, int a) {

	_rank = rank;
	_axis = a;

	int ndim = model->getHyper()->getNdimG1();
	axis ax = model->getHyper()->getAxis(_axis);
	n1 = ax.n;

	howmany = 1;
	for (int i=ndim;i>_rank;i--) howmany *= model->getHyper()->getAxis(i).n;  
	
	/* 
	j-th element of k-th transform is
	out = in + j*stride + k*idist
	*/
 	idist = n1;
 	odist = n1;
	istride = 1; 
	ostride = 1;

	float of = -1.0/(2*ax.d);
	float df = 1.0/(ax.d*n1);

	new_axis = axis(n1,of,df); 	

 	_fwd_exec = false;
 	_adj_exec = false;
}


void FFT1::forward (std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add) {

	fftwf_complex *_in = reinterpret_cast<fftwf_complex*>(model->getVals());
	fftwf_complex *_out = reinterpret_cast<fftwf_complex*>(data->getVals());

 	fftwf_plan _fwd_plan = fftwf_plan_many_dft(_rank, &n1, howmany,
                             _in, NULL, istride, idist,
                             _out, NULL, ostride, odist,
                             FFTW_FORWARD,FFTW_ESTIMATE);

	fftwf_execute(_fwd_plan);

	if (_adj_exec) {
		data->scale(1.0/n1);
		_fwd_exec = false;
	}
	else {
		_fwd_exec = true;
	}

	data->getHyper()->setAxis(_axis,new_axis);

	fftwf_destroy_plan(_fwd_plan);
}

void FFT1::adjoint (std::shared_ptr<complexHyper> model, std::shared_ptr<complexHyper> data, bool add) {

	fftwf_complex *_in = reinterpret_cast<fftwf_complex*>(model->getVals());
	fftwf_complex *_out = reinterpret_cast<fftwf_complex*>(data->getVals());

	fftwf_plan _adj_plan = fftwf_plan_many_dft(_rank, &n1, howmany,
                             _out, NULL, istride, idist,
                             _in, NULL, ostride, odist,
                             FFTW_BACKWARD,FFTW_ESTIMATE);

	fftwf_execute(_adj_plan);
	
	if (_fwd_exec) {
		model->scale(1.0/n1);
		_adj_exec = false;
	}
	else {
		_adj_exec = true;
	}

	fftwf_destroy_plan(_adj_plan);
}

