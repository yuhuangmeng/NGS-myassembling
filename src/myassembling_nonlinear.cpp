#include "myassembling_nonlinear.hpp"

namespace ngcomp
{
  namespace
  {
    void CheckSupportedProblem(shared_ptr<FESpace> fes,
                               shared_ptr<BilinearForm> a)
    {
      if (!fes)
        throw Exception("AssembleNonlinearLocal: fes is null");
      if (!a)
        throw Exception("AssembleNonlinearLocal: bilinear form is null");
      if (a->MixedSpaces())
        throw Exception("AssembleNonlinearLocal: mixed trial/test spaces are not implemented");
      if (a->GetFESpace().get() != fes.get())
        throw Exception("AssembleNonlinearLocal: bilinear form does not use the supplied FESpace");
      if (fes->IsParallel())
        throw Exception("AssembleNonlinearLocal: parallel DoFs are not implemented");
      if (a->GetSpecialElements().Size())
        throw Exception("AssembleNonlinearLocal: special elements are not implemented");
      if (a->UsesEliminateInternal() || a->UsesEliminateHidden() ||
          a->UsesKeepInternal() || a->UsesStoreInner())
        throw Exception("AssembleNonlinearLocal: static condensation/internal elimination is not implemented");

      for (auto bfi : a->Integrators())
        {
          if (bfi->VB() != VOL)
            throw Exception("AssembleNonlinearLocal: only VOL element integrators are implemented");
          if (dynamic_pointer_cast<FacetBilinearFormIntegrator>(bfi))
            throw Exception("AssembleNonlinearLocal: facet/skeleton integrators are not implemented");
        }
    }

    void CheckElementNumber(shared_ptr<MeshAccess> ma, int elnr)
    {
      if (elnr < 0 || elnr >= ma->GetNE(VOL))
        throw Exception("AssembleNonlinearLocal: volume element number out of range: " + ToString(elnr));
    }

    std::vector<int> BuildGlobalToLocal(shared_ptr<FESpace> fes,
                                        const std::vector<int> & dofs)
    {
      std::vector<int> global_to_local(fes->GetNDof(), -1);
      for (size_t i = 0; i < dofs.size(); i++)
        {
          int gdof = dofs[i];
          if (gdof < 0 || gdof >= fes->GetNDof())
            throw Exception("AssembleNonlinearLocal: local dof list contains an invalid global dof");
          global_to_local[gdof] = int(i);
        }
      return global_to_local;
    }

    void CheckSupportMetadata(const std::vector<int> & core_dofs,
                              const std::vector<int> & support_dofs,
                              const std::vector<int> & core_in_support)
    {
      if (core_in_support.size() != core_dofs.size())
        throw Exception("AssembleNonlinearLocal: core_in_support and core_dofs must have the same length");

      for (size_t i = 0; i < core_dofs.size(); i++)
        {
          int local = core_in_support[i];
          if (local < 0 || local >= int(support_dofs.size()))
            throw Exception("AssembleNonlinearLocal: core_in_support index out of range");
          if (support_dofs[local] != core_dofs[i])
            throw Exception("AssembleNonlinearLocal: support_dofs[core_in_support[i]] must equal core_dofs[i]");
        }
    }

    void CheckElementDofsInSupport(const Array<int> & dnums,
                                   const std::vector<int> & global_to_support)
    {
      for (int i = 0; i < dnums.Size(); i++)
        if (dnums[i] >= 0 && global_to_support[dnums[i]] < 0)
          throw Exception("AssembleNonlinearLocal: element dof is missing from supplied support_dofs");
    }

    int LocalDof(int gdof, const std::vector<int> & global_to_local)
    {
      if (gdof < 0)
        return -1;
      if (size_t(gdof) >= global_to_local.size())
        throw Exception("AssembleNonlinearLocal: element dof is outside global_to_local map");
      return global_to_local[gdof];
    }
  }

  shared_ptr<BaseVector>
  MyAssembleNonlinearLocalResidual(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> core_elements,
                                   std::vector<int> support_elements,
                                   std::vector<int> core_dofs,
                                   std::vector<int> support_dofs,
                                   std::vector<int> core_in_support)
  {
    CheckSupportedProblem(fes, a);
    CheckSupportMetadata(core_dofs, support_dofs, core_in_support);
    auto global_to_support = BuildGlobalToLocal(fes, support_dofs);
    auto global_to_core = BuildGlobalToLocal(fes, core_dofs);

    auto ma = fes->GetMeshAccess();
    for (int elnr : core_elements)
      CheckElementNumber(ma, elnr);
    for (int elnr : support_elements)
      CheckElementNumber(ma, elnr);

    int dim = fes->GetDimension();
    auto local_res = make_shared<VVector<double>>(core_dofs.size() * dim);
    local_res->SetScalar(0.0);
    auto local_fv = local_res->FV<double>();

    LocalHeap lh(1000*1000);
    Array<int> dnums;

    for (int elnr : support_elements)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        if (!fes->DefinedOn(ei))
          continue;

        const FiniteElement & fel = fes->GetFE(ei, lh);
        ElementTransformation & eltrans = ma->GetTrafo(ei, lh);
        fes->GetDofNrs(ei, dnums);
        CheckElementDofsInSupport(dnums, global_to_support);

        int elsize = dnums.Size() * dim;
        FlatVector<double> elvecx(elsize, lh);
        FlatVector<double> elvecy(elsize, lh);

        // Matches S_BilinearForm::AddMatrix1: x.GetIndirect(dnums, elvecx).
        u.GetIndirect(dnums, elvecx);

        // Matches S_BilinearForm::AddMatrix1: TransformVec(..., TRANSFORM_SOL).
        fes->TransformVec(ei, elvecx, TRANSFORM_SOL);

        for (auto bfi : a->Integrators())
          {
            if (!bfi->DefinedOn(eltrans.GetElementIndex())) continue;
            if (!bfi->DefinedOnElement(ei.Nr())) continue;

            auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);

            // Matches S_BilinearForm::AddMatrix1: nonlinear element residual.
            bfi->ApplyElementMatrix(fel, mapped_trafo, elvecx, elvecy, 0, lh);

            // Matches S_BilinearForm::AddMatrix1: TransformVec(..., TRANSFORM_RHS).
            fes->TransformVec(ei, elvecy, TRANSFORM_RHS);

            for (int i = 0; i < dnums.Size(); i++)
              {
                int ldof = LocalDof(dnums[i], global_to_core);
                if (ldof < 0)
                  continue;
                for (int c = 0; c < dim; c++)
                  local_fv(dim*ldof+c) += elvecy(dim*i+c);
              }
          }
      }

    return local_res;
  }

  shared_ptr<BaseSparseMatrix>
  MyAssembleNonlinearLocalJacobian(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> core_elements,
                                   std::vector<int> support_elements,
                                   std::vector<int> core_dofs,
                                   std::vector<int> support_dofs,
                                   std::vector<int> core_in_support)
  {
    CheckSupportedProblem(fes, a);
    CheckSupportMetadata(core_dofs, support_dofs, core_in_support);
    auto global_to_support = BuildGlobalToLocal(fes, support_dofs);
    auto global_to_core = BuildGlobalToLocal(fes, core_dofs);

    auto ma = fes->GetMeshAccess();
    for (int elnr : core_elements)
      CheckElementNumber(ma, elnr);
    for (int elnr : support_elements)
      CheckElementNumber(ma, elnr);

    int dim = fes->GetDimension();
    int local_size = int(core_dofs.size()) * dim;

    Array<int> rows, cols;
    Array<double> vals;
    Array<int> dnums;
    LocalHeap lh(1000*1000);

    for (int elnr : support_elements)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        if (!fes->DefinedOn(ei))
          continue;

        const FiniteElement & fel = fes->GetFE(ei, lh);
        ElementTransformation & eltrans = ma->GetTrafo(ei, lh);
        fes->GetDofNrs(ei, dnums);
        CheckElementDofsInSupport(dnums, global_to_support);

        int elsize = dnums.Size() * dim;
        FlatVector<double> elveclin(elsize, lh);
        FlatMatrix<double> elmat(elsize, lh);
        FlatMatrix<double> sum_elmat(elsize, lh);
        sum_elmat = 0.0;

        // Matches S_BilinearForm::AssembleLinearization: lin.GetIndirect(...).
        u.GetIndirect(dnums, elveclin);

        // Matches S_BilinearForm::AssembleLinearization: TransformVec SOL.
        fes->TransformVec(ei, elveclin, TRANSFORM_SOL);

        for (auto bfi : a->Integrators())
          {
            HeapReset hr_bfi(lh);
            if (!bfi->DefinedOn(eltrans.GetElementIndex())) continue;
            if (!bfi->DefinedOnElement(ei.Nr())) continue;

            auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);

            // Matches S_BilinearForm::AssembleLinearization:
            // CalcLinearizedElementMatrix(..., elveclin, elmat, ...).
            bfi->CalcLinearizedElementMatrix(fel, mapped_trafo, elveclin, elmat, lh);
            sum_elmat += elmat;
          }

        // Matches S_BilinearForm::AssembleLinearization: element matrix DoF transform.
        fes->TransformMat(ei, sum_elmat, TRANSFORM_MAT_LEFT_RIGHT);

        for (int i = 0; i < dnums.Size(); i++)
          {
            int li = LocalDof(dnums[i], global_to_core);
            if (li < 0)
              continue;
            for (int j = 0; j < dnums.Size(); j++)
              {
                int lj = LocalDof(dnums[j], global_to_core);
                if (lj < 0)
                  continue;
                for (int ci = 0; ci < dim; ci++)
                  for (int cj = 0; cj < dim; cj++)
                    {
                      rows.Append(dim*li+ci);
                      cols.Append(dim*lj+cj);
                      vals.Append(sum_elmat(dim*i+ci, dim*j+cj));
                    }
              }
          }
      }

    return SparseMatrix<double>::CreateFromCOO(rows, cols, vals,
                                               local_size, local_size);
  }
}
