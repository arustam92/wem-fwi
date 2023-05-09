#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ConformalMap.h"
#include "ComplexLog.h"
#include "ShiftAndScale.h"
#include "MapChain.h"
#include "Rotation.h"

namespace py = pybind11;

using namespace SEP;

PYBIND11_MODULE(pyCMap, clsOps) {
  py::class_<ConformalMap, std::shared_ptr<ConformalMap>>(clsOps, "ConformalMap")
      .def(py::init<std::shared_ptr<hypercube>>(), "Initialize ConformalMap")
      .def("getInHyper",
            (std::shared_ptr<hypercube> (ConformalMap::*)())
            & ConformalMap::getInHyper,
            "Get input hypercube")
      .def("getOutHyper",
            (std::shared_ptr<hypercube> (ConformalMap::*)())
            & ConformalMap::getOutHyper,
            "Get output hypercube");

  py::class_<ComplexLog, std::shared_ptr<ComplexLog>, ConformalMap>(clsOps, "ComplexLog")
      .def(py::init<std::shared_ptr<hypercube>, float>(), "Initialize ComplexLog")
      .def("forward",
            (std::shared_ptr<complex2DReg> (ComplexLog::*)(std::shared_ptr<complex2DReg>))
            & ComplexLog::forward,
            "Forward ComplexLog")
      .def("inverse",
            (std::shared_ptr<complex2DReg> (ComplexLog::*)(std::shared_ptr<complex2DReg>))
            & ComplexLog::inverse,
            "Inverse ComplexLog");

  py::class_<ShiftAndScale, std::shared_ptr<ShiftAndScale>, ConformalMap>(clsOps, "ShiftAndScale")
      .def(py::init<std::shared_ptr<hypercube>, std::complex<float>&, std::complex<float>& >(), "Initialize ShiftAndScale")
      .def("forward",
            (std::shared_ptr<complex2DReg> (ShiftAndScale::*)(std::shared_ptr<complex2DReg>))
            & ShiftAndScale::forward,
            "Forward ShiftAndScale")
      .def("inverse",
            (std::shared_ptr<complex2DReg> (ShiftAndScale::*)(std::shared_ptr<complex2DReg>))
            & ShiftAndScale::inverse,
            "Inverse ShiftAndScale");

  py::class_<MapChain, std::shared_ptr<MapChain>>(clsOps, "MapChain")
      .def(py::init<std::shared_ptr<ConformalMap>, std::shared_ptr<ConformalMap>>(), "Initialize MapChain")
      .def("forward",
            (std::shared_ptr<complex2DReg> (MapChain::*)(std::shared_ptr<complex2DReg>))
            & MapChain::forward,
            "Forward MapChain")
      .def("inverse",
            (std::shared_ptr<complex2DReg> (MapChain::*)(std::shared_ptr<complex2DReg>))
            & MapChain::inverse,
            "Inverse MapChain");

  py::class_<Rotation, std::shared_ptr<Rotation>, ConformalMap>(clsOps, "Rotation")
      .def(py::init<std::shared_ptr<hypercube>, float>(), "Initialize Rotation")
      .def("forward",
            (std::shared_ptr<complex2DReg> (Rotation::*)(std::shared_ptr<complex2DReg>))
            & Rotation::forward,
            "Forward Rotation")
      .def("inverse",
            (std::shared_ptr<complex2DReg> (Rotation::*)(std::shared_ptr<complex2DReg>))
            & Rotation::inverse,
            "Inverse Rotation");

   };
