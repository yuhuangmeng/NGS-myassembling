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

  m.def ("MyBuildLocalSupportPatch",
         &ngcomp::MyBuildLocalSupportPatch,
         py::arg("fes"), py::arg("core_elements"),
         R"raw_string(Build the exact support closure for a core element set.

support_elements are all elements whose dof list intersects core_dofs.
The invariant support_dofs[core_in_support[i]] == core_dofs[i] holds.)raw_string");

  m.def ("MyAssembleLocalSupportMatrix",
         &ngcomp::MyAssembleLocalSupportMatrix,
         py::arg("fes"), py::arg("integrator"), py::arg("core_elements"),
         R"raw_string(Assemble the support-closure local matrix.

The returned support matrix is assembled only over support_elements. Its
core-core block equals the corresponding globally assembled core-core block up
to roundoff.)raw_string");

  m.def ("MyAssembleLocalSupportMatrices",
         &ngcomp::MyAssembleLocalSupportMatrices,
         py::arg("fes"), py::arg("integrator"), py::arg("partition"));
}    
