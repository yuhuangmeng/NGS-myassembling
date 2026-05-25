#ifndef FILE_MYASSEMBLING_COMMON_HPP
#define FILE_MYASSEMBLING_COMMON_HPP

#include <comp.hpp>
#include <string>
#include <vector>

namespace ngcomp
{
  struct LocalMatrix
  {
    shared_ptr<BaseSparseMatrix> mat;
    std::vector<int> core_elements;
    std::vector<int> support_elements;
    std::vector<int> core_dofs;
    std::vector<int> support_dofs;
    std::vector<int> core_in_support;
  };

  struct LocalVector
  {
    shared_ptr<BaseVector> vec;
    std::vector<int> core_elements;
    std::vector<int> support_elements;
    std::vector<int> core_dofs;
    std::vector<int> support_dofs;
    std::vector<int> core_in_support;
  };

  namespace myassembling_detail
  {
    void CheckElementNumbers(shared_ptr<MeshAccess> ma,
                             const std::vector<int> & elements);

    void CheckSupportMetadata(const std::vector<int> & core_dofs,
                              const std::vector<int> & support_dofs,
                              const std::vector<int> & core_in_support);

    std::vector<int> BuildGlobalToLocal(shared_ptr<FESpace> fes,
                                        const std::vector<int> & dofs);

    int LocalDof(int gdof,
                 const std::vector<int> & global_to_local);

    void FillPatchMetadata(LocalMatrix & result,
                           const std::vector<int> & core_elements,
                           const std::vector<int> & support_elements,
                           const std::vector<int> & core_dofs,
                           const std::vector<int> & support_dofs,
                           const std::vector<int> & core_in_support);

    void FillPatchMetadata(LocalVector & result,
                           const std::vector<int> & core_elements,
                           const std::vector<int> & support_elements,
                           const std::vector<int> & core_dofs,
                           const std::vector<int> & support_dofs,
                           const std::vector<int> & core_in_support);
  }
}

#endif
