#ifndef FILE_MYASSEMBLING_NONLINEAR_HPP
#define FILE_MYASSEMBLING_NONLINEAR_HPP

#include "myassembling_common.hpp"

/*
  Local nonlinear assembler for DoF-based local systems.

  Python supplies only local_dofs, which are the overlapping DoFs and therefore
  the rows/columns of the local nonlinear system. C++ computes the finite
  element support closure:

    support_elements = all volume elements touching local_dofs
    support_dofs     = all DoFs appearing on support_elements

  Residual and Jacobian are integrated over support_elements and scattered only
  to local_dofs. The implementation reuses global NGSolve finite elements,
  transformations, DoF numbering, and FESpace DoF transformations.

  First version limitations:
    - VOL element integrators only.
    - No facet/skeleton integrators.
    - No boundary/co-dimension terms.
    - No special elements.
    - No static condensation.
    - No mixed trial/test spaces.
    - No parallel DoF handling.
*/

namespace ngcomp
{
  class LocalNonlinearOperator
  {
  public:
    LocalNonlinearOperator(shared_ptr<FESpace> fes,
                           shared_ptr<BilinearForm> a,
                           std::vector<int> local_dofs,
                           std::vector<int> boundary_dofs,
                           std::vector<double> boundary_values);

    shared_ptr<LocalVector>
    Residual(const BaseVector & u) const;

    shared_ptr<LocalMatrix>
    Jacobian(const BaseVector & u) const;

  private:
    shared_ptr<FESpace> fes_;
    shared_ptr<BilinearForm> a_;
    shared_ptr<MeshAccess> ma_;

    std::vector<int> core_dofs_;
    std::vector<int> support_elements_;
    std::vector<int> support_dofs_;
    std::vector<int> core_in_support_;

    std::vector<int> global_to_core_;
    std::vector<int> global_to_support_;
    std::vector<int> global_to_boundary_;

    std::vector<int> boundary_dofs_;
    std::vector<double> boundary_values_;

    int dim_ = 1;
    int local_size_ = 0;

    int BoundaryIndex(int gdof) const;
    bool IsBoundaryDof(int gdof) const;
  };

  shared_ptr<LocalVector>
  MyAssembleLocalNonlinearResidual(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> local_dofs,
                                   std::vector<int> boundary_dofs,
                                   std::vector<double> boundary_values);

  shared_ptr<LocalMatrix>
  MyAssembleLocalNonlinearJacobian(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> local_dofs,
                                   std::vector<int> boundary_dofs,
                                   std::vector<double> boundary_values);
}

#endif
