#include <complex2DReg.h>
using namespace SEP;

std::shared_ptr<complex2DReg> complex2DReg::clone() const {
  if (getSpaceOnly()) {

    std::shared_ptr<complex2DReg> x(new complex2DReg(getHyper()));
    x->setNorm(getNorm());

    return x;
  } else {
    std::shared_ptr<complex2DReg> x(new complex2DReg(getHyper(), *_mat));
    x->setNorm(getNorm());

    return x;
  }
}
std::shared_ptr<complex2DReg> complex2DReg::cloneSpace() const {
  std::shared_ptr<complex2DReg> x(new complex2DReg(getHyper()));
  x->_mat = 0;
  x->setSpace();
  x->setNorm(getNorm());

  return x;
}

void complex2DReg::initNoData(std::shared_ptr<SEP::hypercube> hyp) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  assert(axes.size() == 2);

  _mat.reset(new complex2D(boost::extents[axes[1].n][axes[0].n]));
  setData(_mat->data());
}
void complex2DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                          const complex2D &vals) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  assert(axes.size() == 2);
  assert(axes[0].n == vals.shape()[1] && axes[1].n == vals.shape()[0]);
  _mat.reset(new complex2D(boost::extents[axes[1].n][axes[0].n]));
  setData(_mat->data());
  for (long long j = 0; j < axes[1].n; j++) {
    for (long long i = 0; i < axes[0].n; i++) {
      (*_mat)[j][i] = vals[j][i];
    }
  }
}

std::shared_ptr<float2DReg> complex2DReg::real() {
  std::shared_ptr<float2DReg> x (new float2DReg(this->getHyper()));
  for (long long i=0; i<this->getHyper()->getAxis(1).n; i++) {
    for (long long j=0; j<this->getHyper()->getAxis(2).n; j++) {
      (*x->_mat)[j][i] = (*this->_mat)[j][i].real();
    }
  }
  return x;
}

std::shared_ptr<float2DReg> complex2DReg::imag() {
  std::shared_ptr<float2DReg> x (new float2DReg(this->getHyper()));
  for (long long i=0; i<this->getHyper()->getAxis(1).n; i++) {
    for (long long j=0; j<this->getHyper()->getAxis(2).n; j++) {
      (*x->_mat)[j][i] = (*this->_mat)[j][i].imag();
    }
  }
  return x;
}
