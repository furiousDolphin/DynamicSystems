
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

#include "System.hpp"

namespace py = pybind11;

PYBIND11_MODULE(module_uno, m) 
{
    m.doc() = "";
    
    py::class_<System>(m, "System")
        .def(py::init<Eigen::VectorXd, Eigen::VectorXd>(),
            py::arg("b_coeffs"),
            py::arg("a_coeffs"))
        .def("step_response", 
             &System::step_response, 
             "oblicza przebieg odpowiedzi skokowej dla podanego wektora czasu",
             py::arg("t_dense"))
        .def("impulse_response", 
             &System::impulse_response, 
             "oblicza przebieg odpowiedzi impulsowej dla podanego wektora czasu",
             py::arg("t_dense"))
        .def("set_forcing_func", 
             &System::set_forcing_func, 
             "",
             py::arg("func") )
        .def("do_RK4_step", 
             &System::do_RK4_step, 
             "",
             py::arg("dt") );
}