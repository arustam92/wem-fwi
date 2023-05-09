#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "WEM.h"
#include "Born_full.h"

namespace py = pybind11;

// namespace gp257 {

using namespace SEP;

PYBIND11_MODULE(pyWEM, clsOps) {
  py::class_<WEM, std::shared_ptr<WEM>>(clsOps, "WEM")
      .def(py::init<std::shared_ptr<float1DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>>(),
          "Initialize NonLinear WEM")
      .def(py::init<std::shared_ptr<float1DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>,
            std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg>>(), "Initialize NonLinear WEM")
      .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>>(),
            "Initialize NonLinear WEM")
      .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>,
          std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg>>(), "Initialize NonLinear WEM")
      .def(py::init<std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>>(), "Initialize Linear WEM")
      .def(py::init<std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>,
            std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg>>(), "Initialize Linear WEM")

      .def("forward",
            (void (WEM::*)(std::shared_ptr<float1DReg>,std::shared_ptr<float3DReg>, bool)) &
            WEM::forward,
            "Forward operator of WEM")
      .def("adjoint",
            (void (WEM::*)(std::shared_ptr<float1DReg>,std::shared_ptr<float3DReg>, bool)) &
            WEM::adjoint,
            "Adjoint operator of WEM")
      .def("forward",
            (void (WEM::*)(std::shared_ptr<float2DReg>,std::shared_ptr<float3DReg>, bool)) &
            WEM::forward,
            "Forward operator of WEM_nl")
      .def("forward",
            (void (WEM::*)(std::shared_ptr<float2DReg>,std::shared_ptr<complex3DReg>, bool)) &
            WEM::forward,
            "Forward operator of complex WEM_nl")

      .def("setFminFmax", (void (WEM::*)(int fmin, int fmax)) & WEM::setFminFmax, "Set limits for extended modeling ");

  py::class_<Born_full, std::shared_ptr<Born_full>>(clsOps, "Born")

      .def(py::init<std::shared_ptr<float1DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>&>(),"Initialize Born")
      .def(py::init<std::shared_ptr<float1DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>,
          std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg> &>(),"Initialize Born")
      // .def(py::init<std::shared_ptr<float3DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>&>(),"Initialize Born")
      // .def(py::init<std::shared_ptr<float3DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>,
      //     std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg> &>(),"Initialize Born")
      .def("forward",
            (void (Born_full::*)(std::shared_ptr<float2DReg>,std::shared_ptr<float3DReg>, bool)) &
            Born_full::forward,
            "Forward operator of Born")
      .def("adjoint",
            (void (Born_full::*)(std::shared_ptr<float2DReg>,std::shared_ptr<float3DReg>, bool)) &
            Born_full::adjoint,
            "Adjoint operator of Born")
      .def("setBgSlow",
            (void (Born_full::*)(std::shared_ptr<float2DReg>)) & Born_full::setBgSlow,
            "Set background slowness")

      .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj> &>(), "Initialize ")
      .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>,
            std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg> &>(), "Initialize ")
      .def("forward",
          (void (Born_full::*)(std::shared_ptr<float2DReg>,std::shared_ptr<complex3DReg>, bool)) &
          Born_full::forward,
          "Forward operator of complex Born")
      .def("adjoint",
          (void (Born_full::*)(std::shared_ptr<float2DReg>,std::shared_ptr<complex3DReg>, bool)) &
          Born_full::adjoint,
         "Adjoint operator of complex Born");

   }
// }  // namespace gp257
