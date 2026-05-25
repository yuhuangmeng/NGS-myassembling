#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/stl.h>

#include "myassembling_linear.hpp"
#include "myassembling_nonlinear.hpp"


PYBIND11_MODULE(myassembling, m)
{
  cout << "Loading myassembling library" << endl;

  py::class_<ngcomp::LocalMatrix, shared_ptr<ngcomp::LocalMatrix>>
    (m, "LocalMatrix")
    .def_readonly("mat", &ngcomp::LocalMatrix::mat)
    .def_readonly("core_elements", &ngcomp::LocalMatrix::core_elements)
    .def_readonly("support_elements", &ngcomp::LocalMatrix::support_elements)
    .def_readonly("core_dofs", &ngcomp::LocalMatrix::core_dofs)
    .def_readonly("support_dofs", &ngcomp::LocalMatrix::support_dofs)
    .def_readonly("core_in_support", &ngcomp::LocalMatrix::core_in_support)
    ;

  py::class_<ngcomp::LocalVector, shared_ptr<ngcomp::LocalVector>>
    (m, "LocalVector")
    .def_readonly("vec", &ngcomp::LocalVector::vec)
    .def_readonly("core_elements", &ngcomp::LocalVector::core_elements)
    .def_readonly("support_elements", &ngcomp::LocalVector::support_elements)
    .def_readonly("core_dofs", &ngcomp::LocalVector::core_dofs)
    .def_readonly("support_dofs", &ngcomp::LocalVector::support_dofs)
    .def_readonly("core_in_support", &ngcomp::LocalVector::core_in_support)
    ;

  m.def ("MyAssembleMatrix",
         &ngcomp::MyAssembleMatrix,
         py::arg("fes"),py::arg("integrator"));

  m.def ("MyAssembleVector",
         &ngcomp::MyAssembleVector,
         py::arg("fes"), py::arg("integrator"));

  m.def ("MyAssembleLocalMatrix",
         &ngcomp::MyAssembleLocalMatrix,
         py::arg("fes"), py::arg("integrator"),
         py::arg("core_elements"), py::arg("support_elements"),
         py::arg("core_dofs"), py::arg("support_dofs"),
         py::arg("core_in_support"),
         R"raw_string(Assemble a Python-supplied local core matrix.

Python constructs the support patch and index maps. C++ assembles over
support_elements and returns the core-core block using core_dofs as the local
numbering.)raw_string");

  m.def ("MyAssembleLocalVector",
         &ngcomp::MyAssembleLocalVector,
         py::arg("fes"), py::arg("integrator"),
         py::arg("core_elements"), py::arg("support_elements"),
         py::arg("core_dofs"), py::arg("support_dofs"),
         py::arg("core_in_support"),
         R"raw_string(Assemble a Python-supplied local core vector.

Python constructs the support patch and index maps. C++ assembles over
support_elements and returns the core entries using core_dofs as the local
numbering.)raw_string");

  m.def ("MyAssembleLocalNonlinearResidual",
         &ngcomp::MyAssembleLocalNonlinearResidual,
         py::arg("fes"), py::arg("a"), py::arg("u"),
         py::arg("core_elements"), py::arg("support_elements"),
         py::arg("core_dofs"), py::arg("support_dofs"),
         py::arg("core_in_support"),
         R"raw_string(Assemble a local nonlinear residual over supplied VOL patch elements.

Python supplies the same support patch metadata as the linear local assembler.
The current solution vector u is a global NGSolve vector. C++ assembles over
support_elements and returns the core entries using core_dofs as the local
numbering.)raw_string");

  m.def ("MyAssembleLocalNonlinearJacobian",
         &ngcomp::MyAssembleLocalNonlinearJacobian,
         py::arg("fes"), py::arg("a"), py::arg("u"),
         py::arg("core_elements"), py::arg("support_elements"),
         py::arg("core_dofs"), py::arg("support_dofs"),
         py::arg("core_in_support"),
         R"raw_string(Assemble a local nonlinear Jacobian over supplied VOL patch elements.

Python supplies the same support patch metadata as the linear local assembler.
The current solution vector u is a global NGSolve vector. C++ assembles over
support_elements and returns the core-core Jacobian using core_dofs as the
local numbering.)raw_string");
}
