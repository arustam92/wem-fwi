#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "FFT1.h"
#include "Taper.h"
#include "LanczosInterpolation2D.h"
#include "LanczosInterpolation3D.h"
#include "Spline2D.h"
#include "Spline3D.h"
#include "Pad.h"
#include "Hilbert.h"

namespace py = pybind11;

// namespace gp257 {

using namespace SEP;

PYBIND11_MODULE(pyOp, clsOps) {
  py::class_<FFT1, std::shared_ptr<FFT1>>(clsOps, "FFT")
      .def(py::init<std::shared_ptr<floatHyper>, std::shared_ptr<complexHyper>, int, int>(), "Initialize FFT")
      .def("forward",
            (void (FFT1::*)(std::shared_ptr<floatHyper>,std::shared_ptr<complexHyper>, bool)) &
            FFT1::forward,
            "Forward operator of FFT")
      .def("adjoint",
            (void (FFT1::*)(std::shared_ptr<floatHyper>,std::shared_ptr<complexHyper>, bool)) &
            FFT1::adjoint,
            "Adjoint operator of FFT")
            
      .def(py::init<std::shared_ptr<complexHyper>, std::string, int, int>(), "Initialize FFT")
      .def("forward",
            (void (FFT1::*)(std::shared_ptr<complexHyper>,std::shared_ptr<complexHyper>, bool)) &
            FFT1::forward,
            "Forward operator of FFT")
      .def("adjoint",
            (void (FFT1::*)(std::shared_ptr<complexHyper>,std::shared_ptr<complexHyper>, bool)) &
            FFT1::adjoint,
            "Adjoint operator of FFT");     

  py::class_<LanczosInterpolation3D, std::shared_ptr<LanczosInterpolation3D>>(clsOps, "LanczosInterpolation3D")
      .def(py::init<std::shared_ptr<complex3DReg>, std::shared_ptr<complex3DReg>, std::vector<float>, std::vector<float>>(), "Initialize LanczosInterpolation3D")
      .def("forward",
            (void (LanczosInterpolation3D::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<complex3DReg>, bool)) &
            LanczosInterpolation3D::forward,
            "Forward operator of LanczosInterpolation3D")
      .def("adjoint",
            (void (LanczosInterpolation3D::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<complex3DReg>, bool)) &
            LanczosInterpolation3D::adjoint,
            "Adjoint operator of LanczosInterpolation3D");

  py::class_<LanczosInterpolation2D, std::shared_ptr<LanczosInterpolation2D>>(clsOps, "LanczosInterpolation2D")
      .def(py::init<std::shared_ptr<complex2DReg>, std::shared_ptr<complex2DReg>, std::vector<float>>(), "Initialize LanczosInterpolation2D")
      .def("forward",
            (void (LanczosInterpolation2D::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            LanczosInterpolation2D::forward,
            "Forward operator of LanczosInterpolation2D")
      .def("adjoint",
            (void (LanczosInterpolation2D::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            LanczosInterpolation2D::adjoint,
              "Adjoint operator of LanczosInterpolation2D");

  py::class_<Spline2D, std::shared_ptr<Spline2D>>(clsOps, "Spline2D")
      .def(py::init<std::shared_ptr<complex2DReg>, std::shared_ptr<complex2DReg>, float, float>(), "Initialize Spline2D")
      .def("forward",
            (void (Spline2D::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            Spline2D::forward,
            "Forward operator of Spline2D")
      .def("adjoint",
            (void (Spline2D::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex2DReg>, bool)) &
            Spline2D::adjoint,
            "Adjoint operator of LanczosInterpolation2D");

  py::class_<Spline3D, std::shared_ptr<Spline3D>>(clsOps, "Spline3D")
      .def(py::init<std::shared_ptr<complex3DReg>, std::shared_ptr<complex3DReg>, float, float, std::vector<float>>(), "Initialize Spline3D")
      .def("forward",
            (void (Spline3D::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<complex3DReg>, bool)) &
            Spline3D::forward,
            "Forward operator of Spline3D")
      .def("adjoint",
            (void (Spline3D::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<complex3DReg>, bool)) &
            Spline3D::adjoint,
            "Adjoint operator of Spline3D");

  py::class_<Pad, std::shared_ptr<Pad>>(clsOps, "Pad")
      .def(py::init<std::shared_ptr<hypercube>, int, bool>(), "Initialize Pad")
      .def("forward",
            (void (Pad::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            Pad::forward,
            "Forward operator of Pad")
      .def("adjoint",
            (void (Pad::*)(std::shared_ptr<complex1DReg>,std::shared_ptr<complex1DReg>, bool)) &
            Pad::adjoint,
            "Adjoint operator of Pad");

  py::class_<Hilbert, std::shared_ptr<Hilbert>>(clsOps, "Hilbert")
      .def(py::init<int>(), "Initialize Hilbert")
      .def("forward",
            (void (Hilbert::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<complex3DReg>, bool)) &
            Hilbert::forward,
            "Forward operator of Hilbert")
      .def("adjoint",
            (void (Hilbert::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<complex3DReg>, bool)) &
            Hilbert::adjoint,
            "Adjoint operator of Hilbert");

   }
// }  // namespace gp257
