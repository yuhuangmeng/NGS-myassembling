#ifndef FILE_MYASSEMBLING_NONLINEAR_HPP
#define FILE_MYASSEMBLING_NONLINEAR_HPP

#include "myassembling_common.hpp"

/*
  Local nonlinear assembler for verification and patch/subdomain experiments.

  This code reproduces the element contribution path of NGSolve's nonlinear
  BilinearForm::Apply and BilinearForm::AssembleLinearization on a selected
  patch. It reuses the global FESpace, global finite elements, global element
  transformations, global DoF numbering, and FESpace DoF transformations.

  It does not define a true local FESpace or a true local GridFunction. Element
  contributions are computed with global NGSolve objects and scattered into the
  local core-DoF numbering supplied by the same patch metadata as the linear
  local assembler.

  First version limitations:
    - VOL element integrators only.
    - No facet/skeleton integrators.
    - No boundary/co-dimension terms.
    - No special elements.
    - No static condensation.
    - No mixed trial/test spaces.
    - No parallel DoF handling.

  Unsupported features throw clear Exceptions instead of silently producing a
  wrong local residual or Jacobian.
*/

namespace ngcomp
{
  shared_ptr<LocalMatrix>
  MyAssembleLocalNonlinearJacobian(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> core_elements,
                                   std::vector<int> support_elements,
                                   std::vector<int> core_dofs,
                                   std::vector<int> support_dofs,
                                   std::vector<int> core_in_support);

  shared_ptr<LocalVector>
  MyAssembleLocalNonlinearResidual(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> core_elements,
                                   std::vector<int> support_elements,
                                   std::vector<int> core_dofs,
                                   std::vector<int> support_dofs,
                                   std::vector<int> core_in_support);
}

#endif
