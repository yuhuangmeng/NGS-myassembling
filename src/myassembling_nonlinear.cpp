#include "myassembling_nonlinear.hpp"

namespace ngcomp
{
  namespace detail = ngcomp::myassembling_detail;

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

    void CheckElementDofsInSupport(const Array<int> & dnums,
                                   const std::vector<int> & global_to_support)
    {
      for (int i = 0; i < dnums.Size(); i++)
        if (dnums[i] >= 0 &&
            detail::LocalDof(dnums[i], global_to_support) < 0)
          throw Exception("AssembleNonlinearLocal: element dof is missing from supplied support_dofs");
    }
  }

  LocalNonlinearOperator ::
  LocalNonlinearOperator(shared_ptr<FESpace> fes_in,
                         shared_ptr<BilinearForm> a_in,
                         std::vector<int> core_elements_in,
                         std::vector<int> support_elements_in,
                         std::vector<int> core_dofs_in,
                         std::vector<int> support_dofs_in,
                         std::vector<int> core_in_support_in)
    : fes(std::move(fes_in)),
      a(std::move(a_in)),
      core_elements(std::move(core_elements_in)),
      support_elements(std::move(support_elements_in)),
      core_dofs(std::move(core_dofs_in)),
      support_dofs(std::move(support_dofs_in)),
      core_in_support(std::move(core_in_support_in))
  {
    CheckSupportedProblem(fes, a);
    detail::CheckSupportMetadata(core_dofs, support_dofs, core_in_support);
    ma = fes->GetMeshAccess();

    detail::CheckElementNumbers(ma, core_elements);
    detail::CheckElementNumbers(ma, support_elements);

    global_to_support = detail::BuildGlobalToLocal(fes, support_dofs);
    global_to_core = detail::BuildGlobalToLocal(fes, core_dofs);

    dim = fes->GetDimension();
    local_size = int(core_dofs.size()) * dim;
  }


  shared_ptr<LocalMatrix>
  LocalNonlinearOperator ::
  Jacobian(const BaseVector & u) const
  {

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
            int li = detail::LocalDof(dnums[i], global_to_core);
            if (li < 0)
              continue;
            for (int j = 0; j < dnums.Size(); j++)
              {
                int lj = detail::LocalDof(dnums[j], global_to_core);
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

    auto result = make_shared<LocalMatrix>();
    result->mat = SparseMatrix<double>::CreateFromCOO(rows, cols, vals,
                                                      local_size, local_size);
    detail::FillPatchMetadata(*result, core_elements, support_elements,
                              core_dofs, support_dofs, core_in_support);

    return result;
  }

  shared_ptr<LocalVector>
  LocalNonlinearOperator ::
  Residual(const BaseVector & u) const
  {
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
                int ldof = detail::LocalDof(dnums[i], global_to_core);
                if (ldof < 0)
                  continue;
                for (int c = 0; c < dim; c++)
                  local_fv(dim*ldof+c) += elvecy(dim*i+c);
              }
          }
      }

    auto result = make_shared<LocalVector>();
    result->vec = local_res;
    detail::FillPatchMetadata(*result, core_elements, support_elements,
                              core_dofs, support_dofs, core_in_support);

    return result;
  }

}
