#include <complex3DReg.h>
using namespace SEP;

std::shared_ptr<complex3DReg> complex3DReg::clone() const {
  if (getSpaceOnly()) {
    std::shared_ptr<complex3DReg> x(new complex3DReg(getHyper()));
    x->setNorm(getNorm());
    return x;
  } else {
    std::shared_ptr<complex3DReg> x(new complex3DReg(getHyper(), *_mat));
    x->setNorm(getNorm());

    return x;
  }
}
std::shared_ptr<complex3DReg> complex3DReg::cloneSpace() const {
  std::shared_ptr<complex3DReg> x(new complex3DReg(getHyper()));
  x->_mat = 0;
  x->setSpace();
  x->setNorm(getNorm());

  return x;
}

void complex3DReg::initNoData(std::shared_ptr<SEP::hypercube> hyp) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  assert(axes.size() == 3);

  _mat.reset(new complex3D(boost::extents[axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
}
void complex3DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                          const complex3D &vals) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  assert(axes.size() == 3);
  assert(axes[0].n == vals.shape()[2] && axes[1].n == vals.shape()[1] &&
         axes[2].n == vals.shape()[0]);
  _mat.reset(new complex3D(boost::extents[axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
  for (long long k = 0; k < axes[2].n; k++) {
    for (long long j = 0; j < axes[1].n; j++) {
      for (long long i = 0; i < axes[0].n; i++) {
        (*_mat)[k][j][i] = vals[k][j][i];
      }
    }
  }
}
void complex3DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                          std::shared_ptr<complex3D> vals) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  assert(axes.size() == 3);
  assert(axes[0].n == vals->shape()[2] && axes[1].n == vals->shape()[1] &&
         axes[2].n == vals->shape()[0]);
  _mat.reset(new complex3D(boost::extents[axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
  for (long long k = 0; k < axes[2].n; k++) {
    for (long long j = 0; j < axes[1].n; j++) {
      for (long long i = 0; i < axes[0].n; i++) {
        (*_mat)[k][j][i] = (*vals)[k][j][i];
      }
    }
  }
}
