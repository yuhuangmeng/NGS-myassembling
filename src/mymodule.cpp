#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/stl.h>

#include "myassembling.hpp"
#include "myassembling_nonlinear.hpp"


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

  py::class_<ngcomp::LocalSupportVector, shared_ptr<ngcomp::LocalSupportVector>>
    (m, "LocalSupportVector")
    .def_readonly("vec", &ngcomp::LocalSupportVector::vec)
    .def_readonly("core_elements", &ngcomp::LocalSupportVector::core_elements)
    .def_readonly("support_elements", &ngcomp::LocalSupportVector::support_elements)
    .def_readonly("core_dofs", &ngcomp::LocalSupportVector::core_dofs)
    .def_readonly("support_dofs", &ngcomp::LocalSupportVector::support_dofs)
    .def_readonly("core_in_support", &ngcomp::LocalSupportVector::core_in_support)
    ;
  
  m.def ("MyAssembleMatrix",
         &ngcomp::MyAssembleMatrix,
         py::arg("fes"),py::arg("integrator"));

  m.def ("MyAssembleVector",
         &ngcomp::MyAssembleVector,
         py::arg("fes"), py::arg("integrator"));

  m.def ("MyAssembleGivenLocalSupportMatrix",
         &ngcomp::MyAssembleGivenLocalSupportMatrix,
         py::arg("fes"), py::arg("integrator"),
         py::arg("core_elements"), py::arg("support_elements"),
         py::arg("core_dofs"), py::arg("support_dofs"),
         py::arg("core_in_support"),
         R"raw_string(Assemble a Python-supplied local support matrix.

Python constructs the support patch and index maps. C++ only assembles over
support_elements using support_dofs as the local numbering.)raw_string");

  m.def ("MyAssembleGivenLocalSupportVector",
         &ngcomp::MyAssembleGivenLocalSupportVector,
         py::arg("fes"), py::arg("integrator"),
         py::arg("core_elements"), py::arg("support_elements"),
         py::arg("core_dofs"), py::arg("support_dofs"),
         py::arg("core_in_support"),
         R"raw_string(Assemble a Python-supplied local support vector.

Python constructs the support patch and index maps. C++ only assembles over
support_elements using support_dofs as the local numbering.)raw_string");

  m.def ("NonlinearAssemblyStatus",
         &ngcomp::NonlinearAssemblyStatus,
         R"raw_string(Return the current status of nonlinear local assembly support.)raw_string");
}    
