#include "myassembling_common.hpp"

#include <set>
#include <unordered_map>
#include <unordered_set>

namespace ngcomp
{
  shared_ptr<LocalSupportInfo>
  BuildLocalSupportInfo(shared_ptr<FESpace> fes,
                        std::vector<int> core_dofs)
  {
    if (!fes)
      throw Exception("BuildLocalSupportInfo: fes is null");

    std::set<int> unique_core_dofs;
    for (int dof : core_dofs)
      {
        if (dof < 0)
          continue;
        if (dof >= fes->GetNDof())
          throw Exception("BuildLocalSupportInfo: core dof is outside the FESpace dof range: " + ToString(dof));
        unique_core_dofs.insert(dof);
      }

    auto info = make_shared<LocalSupportInfo>();
    info->core_dofs = std::vector<int>(unique_core_dofs.begin(), unique_core_dofs.end());
    info->support_elements = myassembling_detail::FindElementsTouchingDofs(fes, info->core_dofs);
    info->support_dofs = myassembling_detail::DofsOfElements(fes, info->support_elements);
    info->core_in_support = myassembling_detail::CoreInSupport(info->core_dofs, info->support_dofs);
    return info;
  }
}

namespace ngcomp::myassembling_detail
{
  void CheckElementNumbers(shared_ptr<MeshAccess> ma,
                           const std::vector<int> & elements)
  {
    for (int elnr : elements)
      if (elnr < 0 || elnr >= ma->GetNE(VOL))
        throw Exception("invalid volume element number: " + ToString(elnr));
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

  std::vector<int> FindElementsTouchingDofs(shared_ptr<FESpace> fes,
                                            const std::vector<int> & dofs)
  {
    std::unordered_set<int> dof_set;
    dof_set.reserve(dofs.size());
    for (int dof : dofs)
      {
        if (dof < 0 || dof >= fes->GetNDof())
          throw Exception("local dof list contains invalid global dof: " + ToString(dof));
        dof_set.insert(dof);
      }

    std::vector<int> elements;
    Array<int> dnums;
    auto ma = fes->GetMeshAccess();

    for (int elnr = 0; elnr < ma->GetNE(VOL); elnr++)
      {
        fes->GetDofNrs(ElementId(VOL, elnr), dnums);
        bool touches = false;
        for (int i = 0; i < dnums.Size(); i++)
          if (dnums[i] >= 0 && dof_set.count(dnums[i]))
            {
              touches = true;
              break;
            }
        if (touches)
          elements.push_back(elnr);
      }

    return elements;
  }

  std::vector<int> DofsOfElements(shared_ptr<FESpace> fes,
                                  const std::vector<int> & elements)
  {
    CheckElementNumbers(fes->GetMeshAccess(), elements);

    std::set<int> dofs;
    Array<int> dnums;
    for (int elnr : elements)
      {
        fes->GetDofNrs(ElementId(VOL, elnr), dnums);
        for (int i = 0; i < dnums.Size(); i++)
          if (dnums[i] >= 0)
            dofs.insert(dnums[i]);
      }

    return std::vector<int>(dofs.begin(), dofs.end());
  }

  std::vector<int> CoreInSupport(const std::vector<int> & core_dofs,
                                 const std::vector<int> & support_dofs)
  {
    std::unordered_map<int, int> support_index;
    support_index.reserve(support_dofs.size());
    for (size_t i = 0; i < support_dofs.size(); i++)
      {
        int dof = support_dofs[i];
        if (support_index.count(dof))
          throw Exception("support_dofs contains duplicate global dof: " + ToString(dof));
        support_index[dof] = int(i);
      }

    std::vector<int> core_in_support;
    core_in_support.reserve(core_dofs.size());
    for (size_t i = 0; i < core_dofs.size(); i++)
      {
        auto it = support_index.find(core_dofs[i]);
        if (it == support_index.end())
          throw Exception("core dof is missing from support_dofs at index " + ToString(i));
        core_in_support.push_back(it->second);
      }

    CheckSupportMetadata(core_dofs, support_dofs, core_in_support);
    return core_in_support;
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
