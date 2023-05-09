#include <float4DReg.h>

using namespace SEP;

std::shared_ptr<float4DReg> float4DReg::clone() const {
  if (getSpaceOnly()) {
    std::shared_ptr<float4DReg> x(new float4DReg(getHyper()));
    x->setNorm(getNorm());

    return x;
  } else {
    std::shared_ptr<float4DReg> x(new float4DReg(getHyper(), *_mat));
    x->setNorm(getNorm());

    return x;
  }
}
std::shared_ptr<float4DReg> float4DReg::cloneSpace() const {
  std::shared_ptr<float4DReg> x(new float4DReg(getHyper()));
  x->_mat = 0;
  x->setNorm(getNorm());

  x->setSpace();
  return x;
}

void float4DReg::initNoData(std::shared_ptr<SEP::hypercube> hyp) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  assert(axes.size() == 4);

  _mat.reset(
      new float4D(boost::extents[axes[3].n][axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
}
void float4DReg::initData(std::shared_ptr<SEP::hypercube> hyp,
                          const float4D &vals) {
  const std::vector<SEP::axis> axes = hyp->getAxes();
  setHyper(hyp);

  assert(axes.size() == 4);
  assert(axes[0].n == vals.shape()[3] && axes[1].n == vals.shape()[2] &&
         axes[2].n == vals.shape()[1] && axes[3].n == vals.shape()[0]);
  _mat.reset(
      new float4D(boost::extents[axes[3].n][axes[2].n][axes[1].n][axes[0].n]));
  setData(_mat->data());
  for (long long l = 0; l < axes[3].n; l++) {
    for (long long k = 0; k < axes[2].n; k++) {
      for (long long j = 0; j < axes[1].n; j++) {
        for (long long i = 0; i < axes[0].n; i++) {
          (*_mat)[l][k][j][i] = vals[l][k][j][i];
        }
      }
    }
  }
}
