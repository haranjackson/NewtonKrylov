#include <vector>

#include "../include/pybind11/eigen.h"
#include "../include/pybind11/functional.h"
#include "../include/pybind11/pybind11.h"
#include "../include/pybind11/stl.h"

#include "../newton_krylov.h"

namespace py = pybind11;


PYBIND11_PLUGIN(NewtonKrylov)
{
    py::module m("NewtonKrylov",
                 "Python bindings to the Newton-Krylov C++ implementation");

    m.def("solve",
          &nonlin_solve,
          py::arg("F"),
          py::arg("x"),
          py::arg("f_tol")=pow(mEPS,1./3),
          py::arg("f_rtol")=INF,
          py::arg("x_tol")=INF,
          py::arg("x_rtol")=INF);

    return m.ptr();
}
