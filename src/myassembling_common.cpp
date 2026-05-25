#include "myassembling_common.hpp"

namespace ngcomp::myassembling_detail
{
  void CheckElementNumbers(shared_ptr<MeshAccess> ma,
                           const std::vector<int> & elements,
                           const std::string & name)
  {
    for (int elnr : elements)
      if (elnr < 0 || elnr >= ma->GetNE(VOL))
        throw Exception(name + " contains invalid volume element number: " + ToString(elnr));
  }

  void CheckSupportMetadata(const std::vector<int> & core_dofs,
                            const std::vector<int> & support_dofs,
                            const std::vector<int> & core_in_support)
  {
    if (core_in_support.size() != core_dofs.size())
      throw Exception("core_in_support and core_dofs must have the same length");

    for (size_t i = 0; i < core_dofs.size(); i++)
      {
        int local = core_in_support[i];
        if (local < 0 || local >= int(support_dofs.size()))
          throw Exception("core_in_support index out of range at index " + ToString(i));
        if (support_dofs[local] != core_dofs[i])
          throw Exception("support_dofs[core_in_support[i]] must equal core_dofs[i] at index " + ToString(i));
      }
  }

  std::vector<int> BuildGlobalToLocal(shared_ptr<FESpace> fes,
                                      const std::vector<int> & dofs)
  {
    std::vector<int> global_to_local(fes->GetNDof(), -1);
    for (size_t i = 0; i < dofs.size(); i++)
      {
        int gdof = dofs[i];
        if (gdof < 0 || gdof >= fes->GetNDof())
          throw Exception("local dof list contains invalid global dof: " + ToString(gdof));
        if (global_to_local[gdof] != -1)
          throw Exception("local dof list contains duplicate global dof: " + ToString(gdof));
        global_to_local[gdof] = int(i);
      }
    return global_to_local;
  }

  int LocalDof(int gdof,
               const std::vector<int> & global_to_local)
  {
    if (gdof < 0 || size_t(gdof) >= global_to_local.size())
      return -1;
    return global_to_local[gdof];
  }

  void FillPatchMetadata(LocalMatrix & result,
                         const std::vector<int> & core_elements,
                         const std::vector<int> & support_elements,
                         const std::vector<int> & core_dofs,
                         const std::vector<int> & support_dofs,
                         const std::vector<int> & core_in_support)
  {
    result.core_elements = core_elements;
    result.support_elements = support_elements;
    result.core_dofs = core_dofs;
    result.support_dofs = support_dofs;
    result.core_in_support = core_in_support;
  }

  void FillPatchMetadata(LocalVector & result,
                         const std::vector<int> & core_elements,
                         const std::vector<int> & support_elements,
                         const std::vector<int> & core_dofs,
                         const std::vector<int> & support_dofs,
                         const std::vector<int> & core_in_support)
  {
    result.core_elements = core_elements;
    result.support_elements = support_elements;
    result.core_dofs = core_dofs;
    result.support_dofs = support_dofs;
    result.core_in_support = core_in_support;
  }
}
