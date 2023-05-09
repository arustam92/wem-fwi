#include <pybind11/chrono.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "WEM.h"
#include "BornFull.h"
#include "BornRefl.h"
#include "BornTomo.h"
#include "extWEM.h"
#include "extBornFull.h"

namespace py = pybind11;

// namespace gp257 {

using namespace SEP;

PYBIND11_MODULE(pyWEM, clsOps) {
  py::class_<WEM, std::shared_ptr<WEM>>(clsOps, "WEM")
      .def(py::init<std::shared_ptr<float1DReg>, std::shared_ptr<hypercube>, std::shared_ptr<paramObj> &>(),
          "Initialize NonLinear WEM")
      .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<hypercube>, std::shared_ptr<paramObj> &>(),
            "Initialize NonLinear WEM")

      .def("forward",
            (void (WEM::*)(const std::vector<std::shared_ptr<complex2DReg>>&,std::shared_ptr<float3DReg>, bool)) &
            WEM::forward,
            "Forward operator of WEM_nl")
      .def("forward",
            (void (WEM::*)(const std::vector<std::shared_ptr<complex2DReg>>&,std::shared_ptr<complex3DReg>, bool)) &
            WEM::forward,
            "Forward operator of complex WEM_nl")

      .def("get_nfreq",
            (int (WEM::*)()) &
            WEM::get_nfreq,
            "Get number of propagated frequencies");
      // .def("wavefield",
      //       (void (WEM::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<float4DReg>, bool)) &
      //       WEM::wavefield,
      //       "Snap wavefield");

  py::class_<BornFull, std::shared_ptr<BornFull>>(clsOps, "Born")

      .def(py::init<std::shared_ptr<float1DReg>, std::vector<std::shared_ptr<complex2DReg>>, std::shared_ptr<paramObj> &>(),"Initialize Born")
      // .def(py::init<std::shared_ptr<float3DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>&>(),"Initialize Born")
      // .def(py::init<std::shared_ptr<float3DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>,
      //     std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg> &>(),"Initialize Born")
      .def("forward",
            (void (BornFull::*)(const std::vector<std::shared_ptr<complex2DReg>>&,std::shared_ptr<float3DReg>, bool)) &
            BornFull::forward,
            "Forward operator of Born")
      .def("adjoint",
            (void (BornFull::*)(std::vector<std::shared_ptr<complex2DReg>>&,std::shared_ptr<float3DReg>, bool)) &
            BornFull::adjoint,
            "Adjoint operator of Born")
      .def("forward",
            (void (BornFull::*)(const std::vector<std::shared_ptr<complex2DReg>>&,std::shared_ptr<complex3DReg>, bool)) &
            BornFull::forward,
            "Forward operator of Born")
      .def("adjoint",
            (void (BornFull::*)(std::vector<std::shared_ptr<complex2DReg>>&,std::shared_ptr<complex3DReg>, bool)) &
            BornFull::adjoint,
            "Adjoint operator of Born")
      .def("setBgSlow",
            (void (BornFull::*)(std::vector<std::shared_ptr<complex2DReg>>&)) & BornFull::setBgSlow,
            "Set background slowness");

      // .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<complex2DReg>, std::shared_ptr<paramObj> &>(), "Initialize ")
      // .def("forward",
      //     (void (BornFull::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex3DReg>, bool)) &
      //     BornFull::forward,
      //     "Forward operator of complex Born")
      // .def("adjoint",
      //     (void (BornFull::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex3DReg>, bool)) &
      //     BornFull::adjoint,
      //    "Adjoint operator of complex Born");



 py::class_<BornRefl, std::shared_ptr<BornRefl>>(clsOps, "BornRefl")

     .def(py::init<std::shared_ptr<float1DReg>, std::shared_ptr<complex2DReg>, std::shared_ptr<paramObj> &>(), "Initialize ")
     .def("forward",
         (void (BornRefl::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<float3DReg>, bool)) &
         BornRefl::forward,
         "Forward operator of complex Born")
     .def("adjoint",
         (void (BornRefl::*)(std::shared_ptr<complex3DReg>,std::shared_ptr<float3DReg>, bool)) &
         BornRefl::adjoint,
        "Adjoint operator of complex Born")

      .def("setBgSlow",
            (void (BornRefl::*)(std::shared_ptr<complex2DReg>)) & BornRefl::setBgSlow,
            "Set background slowness");

//      .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<complex2DReg>, std::shared_ptr<paramObj> &>(), "Initialize ")
//      .def("forward",
//          (void (BornRefl::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex3DReg>, bool)) &
//          BornRefl::forward,
//          "Forward operator of complex Born")
//      .def("adjoint",
//          (void (BornRefl::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex3DReg>, bool)) &
//          BornRefl::adjoint,
//         "Adjoint operator of complex Born");



py::class_<BornTomo, std::shared_ptr<BornTomo>>(clsOps, "BornTomo")

    .def(py::init<std::shared_ptr<float1DReg>, std::shared_ptr<complex2DReg>, std::shared_ptr<paramObj> &>(), "Initialize ")
    .def("setBgSlow",
          (void (BornTomo::*)(std::shared_ptr<complex2DReg>)) & BornTomo::setBgSlow,
          "Set background slowness")
     .def("setBgRefl",
          (void (BornTomo::*)(std::shared_ptr<complex3DReg>)) & BornTomo::setBgRefl,
          "Set background reflectivity")
    .def("forward",
        (void (BornTomo::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<float3DReg>, bool)) &
        BornTomo::forward,
        "Forward operator of complex Born")
    .def("adjoint",
        (void (BornTomo::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<float3DReg>, bool)) &
        BornTomo::adjoint,
       "Adjoint operator of complex Born");

    // .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<complex2DReg>, std::shared_ptr<paramObj>,bool &>(), "Initialize ")
    // .def(py::init<std::shared_ptr<complex1DReg>, std::shared_ptr<complex2DReg>, std::shared_ptr<paramObj>,bool,
    //       std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg> &>(), "Initialize ")
    // .def("forward",
    //     (void (BornTomo::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex3DReg>, bool)) &
    //     BornTomo::forward,
    //     "Forward operator of complex Born")
    // .def("adjoint",
    //     (void (BornTomo::*)(std::shared_ptr<complex2DReg>,std::shared_ptr<complex3DReg>, bool)) &
    //     BornTomo::adjoint,
    //    "Adjoint operator of complex Born");


 py::class_<extWEM, std::shared_ptr<extWEM>>(clsOps, "extWEM")
     .def(py::init<std::shared_ptr<float1DReg>, std::shared_ptr<hypercube>, std::shared_ptr<paramObj>>(),
         "Initialize NonLinear WEM")
     .def("forward",
           (void (extWEM::*)(const std::vector<std::shared_ptr<complex3DReg>>&,std::shared_ptr<float3DReg>, bool)) &
           extWEM::forward,
           "Forward operator of WEM_nl")
     .def("forward",
           (void (extWEM::*)(const std::vector<std::shared_ptr<complex3DReg>>&,std::shared_ptr<complex3DReg>, bool)) &
           extWEM::forward,
           "Forward operator of complex WEM_nl")
     .def("wavefield",
           (void (extWEM::*)(const std::vector<std::shared_ptr<complex3DReg>>&,std::shared_ptr<float4DReg>, bool)) &
           extWEM::wavefield,
           "Snap wavefield of WEM_nl");


 py::class_<extBornFull, std::shared_ptr<extBornFull>>(clsOps, "extBorn")

     .def(py::init<std::shared_ptr<float1DReg>, std::vector<std::shared_ptr<complex3DReg>>, std::shared_ptr<paramObj> &>(),"Initialize Born")

     // .def(py::init<std::shared_ptr<float3DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>&>(),"Initialize Born")
     // .def(py::init<std::shared_ptr<float3DReg>, std::shared_ptr<float2DReg>, std::shared_ptr<paramObj>,
     //     std::shared_ptr<float2DReg>, std::shared_ptr<float2DReg> &>(),"Initialize Born")
     .def("forward",
           (void (extBornFull::*)(const std::vector<std::shared_ptr<complex3DReg>>&,std::shared_ptr<float3DReg>, bool)) &
           extBornFull::forward,
           "Forward operator of Born")
     .def("adjoint",
           (void (extBornFull::*)(std::vector<std::shared_ptr<complex3DReg>>&,std::shared_ptr<float3DReg>, bool)) &
           extBornFull::adjoint,
           "Adjoint operator of Born")

      .def("forward",
           (void (extBornFull::*)(const std::vector<std::shared_ptr<complex3DReg>>&,std::shared_ptr<complex3DReg>, bool)) &
           extBornFull::forward,
           "Forward operator of Born")
     .def("adjoint",
           (void (extBornFull::*)(std::vector<std::shared_ptr<complex3DReg>>&,std::shared_ptr<complex3DReg>, bool)) &
           extBornFull::adjoint,
           "Adjoint operator of Born")

     .def("setBgSlow",
           (void (extBornFull::*)(const std::vector<std::shared_ptr<complex3DReg>>&)) & extBornFull::setBgSlow,
           "Set background slowness");

   }
// }  // namespace gp257
