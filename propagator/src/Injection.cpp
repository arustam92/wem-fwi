#include <Injection.h>
#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

using namespace SEP;

void Injection::forward(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) data->scale(0.);
		int nx = data->getHyper()->getAxis(1).n;

		for (int i=0; i<_ix.size(); i++) {
			for (int ig=-_ng; ig<=_ng; ig++) {
				int off = _ix[i]+ig+_tap;
				if (off >= 0 && off < nx)
					(*data->_mat)[_iz[i]][off] += (*gauss)[ig+_ng]*(*model->_mat)[_istep];
			}
		}
}

void Injection::adjoint(std::shared_ptr<complex1DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
		if (!add) model->scale(0.);
		int nx = data->getHyper()->getAxis(1).n;

		for (int i=0; i<_ix.size(); i++) {
			for (int ig=-_ng; ig<=_ng; ig++) {
				int off = _ix[i]+ig+_tap;
				if (off >= 0 && off < nx)
					(*model->_mat)[_istep] += (*gauss)[ig+_ng]*(*data->_mat)[_iz[i]][off];
			}
		}
}

void Injection::forward(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	if(!add) data->scale(0.);

	tbb::parallel_for(tbb::blocked_range<int>(0,_ix.size()),
		[&] (const tbb::blocked_range<int> &r) {
		for (int i=r.begin(); i<r.end(); i++) {
			 (*data->_mat)[_iz[i]][_ix[i]] += (*model->_mat)[i][_istep];
		}
	});
}

void Injection::adjoint(std::shared_ptr<complex2DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	if(!add) model->scale(0.);

	tbb::parallel_for(tbb::blocked_range<int>(0,_ix.size()),
		[&] (const tbb::blocked_range<int> &r) {
		for (int i=r.begin(); i<r.end(); i++) {
			(*model->_mat)[i][_istep] += (*data->_mat)[_iz[i]][_ix[i]];
		}
	});
}

void Injection::forward(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	if (!add) data->scale(0.);
	int nx = data->getHyper()->getAxis(1).n;

	for (int ig=-_ng; ig<=_ng; ig++) {
		for (int i=0; i<_ix.size(); i++) {
				int off = _ix[i]+ig+_tap;
				if (off >= 0 && off < nx)
					(*data->_mat)[_iz[i]][off] += (*gauss)[ig+_ng] * (*model->_mat)[_ishot][i][_istep];
			}
		}
}

void Injection::adjoint(std::shared_ptr<complex3DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	if (!add) model->scale(0.);
	int nx = data->getHyper()->getAxis(1).n;
	
	for (int ig=-_ng; ig<=_ng; ig++) {
		for (int i=0; i<_ix.size(); i++) {		
				int off = _ix[i]+ig+_tap;
				if (off >= 0 && off < nx)
					(*model->_mat)[_ishot][i][_istep] += (*gauss)[ig+_ng] * (*data->_mat)[_iz[i]][off];
			}
		}
}

void Injection::adjoint(std::shared_ptr<complex4DReg> model, std::shared_ptr<complex2DReg> data, bool add) {
	if (!add) model->scale(0.);

	axis ax2 = model->getHyper()->getAxis(3);
	axis ax1 = model->getHyper()->getAxis(2);
		for (int i2=0; i2<ax2.n; ++i2) {
			for (int i1=0; i1<ax1.n; ++i1) {
				(*model->_mat)[_ishot][i2][i1][_istep] += (*data->_mat)[i2][i1+_tap];
		}
	}
}
