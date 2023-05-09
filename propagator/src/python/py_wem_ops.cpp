#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "SplitStep.h"
#include "Scatter.h"
#include "SSF.h"
#include "Reflect.h"
#include "RefSampler.h"
#include "Phshift.h"
#include "OneWay.h"
#include "LinOneWay.h"
#include "dReflect.h"
#include "Injection.h"

namespace py = pybind11;

// namespace gp257 {

using namespace SEP;

PYBIND11_MODULE(py_wem_ops, clsOps) {
  py::class_<SplitStep, std::shared_ptr<SplitStep>>(clsOps, "SplitStep")
      .def(py::init<std::shared_ptr<complex2DReg> &>())

      .def("forward",
            (void (SplitStep::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            SplitStep::forward)
      .def("adjoint",
            (void (SplitStep::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            SplitStep::adjoint)

      .def("setFreq",
            (void (SplitStep::*)(float&)) &
            SplitStep::setFreq)
      .def("setLocation",
            (void (SplitStep::*)(std::vector<int>& )) &
            SplitStep::setLocation)
      .def("setSlow",
            (void (SplitStep::*)(std::complex<float>&)) &
            SplitStep::setSlow)
      .def("setDepth",
            (void (SplitStep::*)(int)) &
                  SplitStep::setDepth);


  py::class_<Scatter, std::shared_ptr<Scatter>>(clsOps, "Scatter")
      .def(py::init<std::shared_ptr<complex2DReg>, int, std::shared_ptr<paramObj> &>())
      .def("forward",
            (void (Scatter::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            Scatter::forward)
      .def("adjoint",
            (void (Scatter::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            Scatter::adjoint)
      .def("setFreq",
            (void (Scatter::*)(float)) &
            Scatter::setFreq)
      .def("setDepth",
            (void (Scatter::*)(int)) &
                  Scatter::setDepth);

  py::class_<SSF, std::shared_ptr<SSF>>(clsOps, "SSF")
      .def(py::init<std::shared_ptr<complex2DReg>, std::shared_ptr<paramObj>, std::shared_ptr<RefSampler> &>())
      .def("forward",
            (void (SSF::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            SSF::forward)
      .def("adjoint",
            (void (SSF::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            SSF::adjoint)
      .def("setFreq",
            (void (SSF::*)(float)) &
            SSF::setFreq)
      .def("setDepth",
            (void (SSF::*)(int)) &
                  SSF::setDepth);

  py::class_<Reflect, std::shared_ptr<Reflect>>(clsOps, "Reflect")
      .def(py::init<std::vector<std::shared_ptr<complex2DReg>> &>())
      .def("forward",
            (void (Reflect::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            Reflect::forward)
      .def("adjoint",
            (void (Reflect::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            Reflect::adjoint);

  py::class_<RefSampler, std::shared_ptr<RefSampler>>(clsOps, "RefSampler")
      .def(py::init<std::shared_ptr<complex2DReg>, int &>())
      .def("getRefSlow",
            (std::complex<float>& (RefSampler::*)(int, int)) &
            RefSampler::getRefSlow)
      ;

  py::class_<Phshift, std::shared_ptr<Phshift>>(clsOps, "Phshift")
      .def(py::init<float,std::shared_ptr<hypercube>,int,float &>())
      .def(py::init<float,std::shared_ptr<hypercube>,int &>())
      .def("forward",
            (void (Phshift::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            Phshift::forward)
      .def("adjoint",
            (void (Phshift::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            Phshift::adjoint)

      .def("setFreq",
            (void (Phshift::*)(float)) &
            Phshift::setFreq)
      .def("setSlow",
            (void (Phshift::*)(std::complex<float>)) &
            Phshift::setSlow)
      .def("setRef",
            (void (Phshift::*)(int)) &
                  Phshift::setRef);

  py::class_<OneWay, std::shared_ptr<OneWay>>(clsOps, "OneWay")
      .def(py::init<std::shared_ptr<complex2DReg>,std::shared_ptr<paramObj>, std::shared_ptr<RefSampler> &>());

  py::class_<Down, std::shared_ptr<Down>, OneWay>(clsOps, "Down")
      .def(py::init<std::shared_ptr<complex2DReg>,std::shared_ptr<paramObj>, std::shared_ptr<RefSampler> &>())
      .def("forward",
            (void (Down::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            Down::forward)
      .def("adjoint",
            (void (Down::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            Down::adjoint)
      .def("setFreq",
            (void (Down::*)(float)) &
            Down::setFreq);

  py::class_<Up, std::shared_ptr<Up>, OneWay>(clsOps, "Up")
      .def(py::init<std::shared_ptr<complex2DReg>,std::shared_ptr<paramObj>, std::shared_ptr<RefSampler> &>())
      .def("forward",
            (void (Up::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            Up::forward)
      .def("adjoint",
            (void (Up::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            Up::adjoint)
      .def("setFreq",
            (void (Up::*)(float)) &
            Up::setFreq);

  py::class_<IC, std::shared_ptr<IC>>(clsOps, "IC")
      .def(py::init<const std::shared_ptr<complex2DReg>&>())
      .def("forward",
            (void (IC::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            IC::forward)
      .def("adjoint",
            (void (IC::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            IC::adjoint)
      .def("forward",
            (void (IC::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            IC::forward)
      .def("adjoint",
            (void (IC::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            IC::adjoint)
      .def("setDepth",
            (void (IC::*)(int)) &
            IC::setDepth);

  py::class_<LinDown, std::shared_ptr<LinDown>>(clsOps, "LinDown")
      .def(py::init<std::shared_ptr<complex2DReg>&,std::shared_ptr<paramObj>&, std::shared_ptr<complex2DReg>&, std::shared_ptr<OneWay> &>())
      .def("forward",
            (void (LinDown::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            LinDown::forward)
      .def("adjoint",
            (void (LinDown::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            LinDown::adjoint)
      .def("setFreq",
            (void (LinDown::*)(float)) &
            LinDown::setFreq);

  py::class_<dReflect, std::shared_ptr<dReflect>>(clsOps, "dReflect")
      .def(py::init<std::vector<std::shared_ptr<complex2DReg>>&>())
      .def("forward",
            (void (dReflect::*)(std::vector<std::shared_ptr<complex2DReg>>,std::shared_ptr<complex2DReg>, bool)) &
            dReflect::forward)
      .def("adjoint",
            (void (dReflect::*)(std::vector<std::shared_ptr<complex2DReg>>,std::shared_ptr<complex2DReg>, bool)) &
            dReflect::adjoint);

      py::class_<Injection, std::shared_ptr<Injection>>(clsOps, "Injection")
                              .def(py::init<const std::vector<int>&, const std::vector<int>&, int, int>(), "Initialize Injection")
                              .def("forward",
                                                (void (Injection::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex2DReg>, bool)) &
                                                Injection::forward,
                                                "Forward operator of Injection")
                              .def("adjoint",
                                                (void (Injection::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex2DReg>, bool)) &
                                                Injection::adjoint,
                                                "Adjoint operator of Injection")
                              .def("forward",
                                                (void (Injection::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<complex2DReg>, bool)) &
                                                Injection::forward,
                                                "Forward operator of Injection")
                              .def("adjoint",
                                                (void (Injection::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<complex2DReg>, bool)) &
                                                Injection::adjoint,
                                                "Adjoint operator of Injection")
                              .def("setStep",
                                                (void (Injection::*)(int)) &
                                                Injection::setStep,
                                                "Set step")
                              .def("setShot",
                                                (void (Injection::*)(int)) &
                                                Injection::setShot,
                                                "Set shot");

   }
// }  // namespace gp257
