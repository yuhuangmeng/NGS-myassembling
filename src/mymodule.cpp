#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/stl.h>

#include "myassembling.hpp"


PYBIND11_MODULE(myassembling, m)
{
  cout << "Loading myassembling library" << endl;

  py::class_<ngcomp::LocalSupportMatrix, shared_ptr<ngcomp::LocalSupportMatrix>>
    (m, "LocalSupportMatrix")
    .def_readonly("mat", &ngcomp::LocalSupportMatrix::mat)
    .def_readonly("core_elements", &ngcomp::LocalSupportMatrix::core_elements)
    .def_readonly("support_elements", &ngcomp::LocalSupportMatrix::support_elements)
    .def_readonly("core_dofs", &ngcomp::LocalSupportMatrix::core_dofs)
    .def_readonly("support_dofs", &ngcomp::LocalSupportMatrix::support_dofs)
    .def_readonly("core_in_support", &ngcomp::LocalSupportMatrix::core_in_support)
    ;
  
  m.def ("MyAssembleMatrix",
         &ngcomp::MyAssembleMatrix,
         py::arg("fes"),py::arg("integrator"));

  m.def ("MyAssembleGivenLocalSupportMatrix",
         &ngcomp::MyAssembleGivenLocalSupportMatrix,
         py::arg("fes"), py::arg("integrator"),
         py::arg("core_elements"), py::arg("support_elements"),
         py::arg("core_dofs"), py::arg("support_dofs"),
         py::arg("core_in_support"),
         R"raw_string(Assemble a Python-supplied local support matrix.

Python constructs the support patch and index maps. C++ only assembles over
support_elements using support_dofs as the local numbering.)raw_string");
}    
