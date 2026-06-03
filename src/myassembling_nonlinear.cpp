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
          throw Exception("AssembleNonlinearLocal: element dof is missing from computed support_dofs");
    }
  }

  LocalNonlinearOperator ::
  LocalNonlinearOperator(shared_ptr<FESpace> fes,
                         shared_ptr<BilinearForm> a,
                         std::vector<int> local_dofs)
    : fes_(std::move(fes)),
      a_(std::move(a))
  {
    CheckSupportedProblem(fes_, a_);

    ma_ = fes_->GetMeshAccess();

    auto info = BuildLocalSupportInfo(fes_, std::move(local_dofs));
    core_dofs_ = info->core_dofs;
    support_dofs_ = info->support_dofs;
    support_elements_ = info->support_elements;
    core_in_support_ = info->core_in_support;

    global_to_core_ = detail::BuildGlobalToLocal(fes_, core_dofs_);
    global_to_support_ = detail::BuildGlobalToLocal(fes_, support_dofs_);

    dim_ = fes_->GetDimension();
    local_size_ = int(core_dofs_.size()) * dim_;
  }


  shared_ptr<LocalVector>
  LocalNonlinearOperator ::
  Residual(const BaseVector & u) const
  {
    auto local_res = make_shared<VVector<double>>(local_size_);
    local_res->SetScalar(0.0);
    auto local_fv = local_res->FV<double>();

    LocalHeap lh(1000*1000);
    Array<int> dnums;

    for (int elnr : support_elements_)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        if (!fes_->DefinedOn(ei))
          continue;

        const FiniteElement & fel = fes_->GetFE(ei, lh);
        ElementTransformation & eltrans = ma_->GetTrafo(ei, lh);
        fes_->GetDofNrs(ei, dnums);
        CheckElementDofsInSupport(dnums, global_to_support_);

        int elsize = dnums.Size() * dim_;
        FlatVector<double> elvecx(elsize, lh);
        FlatVector<double> elvecy(elsize, lh);

        // Matches S_BilinearForm::AddMatrix1: x.GetIndirect(dnums, elvecx).
        u.GetIndirect(dnums, elvecx);

        // Matches S_BilinearForm::AddMatrix1: TransformVec(..., TRANSFORM_SOL).
        fes_->TransformVec(ei, elvecx, TRANSFORM_SOL);

        for (auto bfi : a_->Integrators())
          {
            if (!bfi->DefinedOn(eltrans.GetElementIndex())) continue;
            if (!bfi->DefinedOnElement(ei.Nr())) continue;

            auto & mapped_trafo = eltrans.AddDeformation(bfi->GetDeformation().get(), lh);

            // Matches S_BilinearForm::AddMatrix1: nonlinear element residual.
            bfi->ApplyElementMatrix(fel, mapped_trafo, elvecx, elvecy, 0, lh);

            // Matches S_BilinearForm::AddMatrix1: TransformVec(..., TRANSFORM_RHS).
            fes_->TransformVec(ei, elvecy, TRANSFORM_RHS);

            for (int i = 0; i < dnums.Size(); i++)
              {
                int ldof = detail::LocalDof(dnums[i], global_to_core_);
                if (ldof < 0)
                  continue;
                for (int c = 0; c < dim_; c++)
                  local_fv(dim_*ldof+c) += elvecy(dim_*i+c);
              }
          }
      }

    auto result = make_shared<LocalVector>();
    result->vec = local_res;
    detail::FillPatchMetadata(*result, support_elements_, support_elements_,
                              core_dofs_, support_dofs_, core_in_support_);

    return result;
  }


  shared_ptr<LocalMatrix>
  LocalNonlinearOperator ::
  Jacobian(const BaseVector & u) const
  {
    Array<int> rows, cols;
    Array<double> vals;
    Array<int> dnums;
    LocalHeap lh(1000*1000);

    for (int elnr : support_elements_)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        if (!fes_->DefinedOn(ei))
          continue;

        const FiniteElement & fel = fes_->GetFE(ei, lh);
        ElementTransformation & eltrans = ma_->GetTrafo(ei, lh);
        fes_->GetDofNrs(ei, dnums);
        CheckElementDofsInSupport(dnums, global_to_support_);

        int elsize = dnums.Size() * dim_;
        FlatVector<double> elveclin(elsize, lh);
        FlatMatrix<double> elmat(elsize, lh);
        FlatMatrix<double> sum_elmat(elsize, lh);
        sum_elmat = 0.0;

        // Matches S_BilinearForm::AssembleLinearization: lin.GetIndirect(...).
        u.GetIndirect(dnums, elveclin);

        // Matches S_BilinearForm::AssembleLinearization: TransformVec SOL.
        fes_->TransformVec(ei, elveclin, TRANSFORM_SOL);

        for (auto bfi : a_->Integrators())
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
        fes_->TransformMat(ei, sum_elmat, TRANSFORM_MAT_LEFT_RIGHT);

        for (int i = 0; i < dnums.Size(); i++)
          {
            int li = detail::LocalDof(dnums[i], global_to_core_);
            if (li < 0)
              continue;
            for (int j = 0; j < dnums.Size(); j++)
              {
                int lj = detail::LocalDof(dnums[j], global_to_core_);
                if (lj < 0)
                  continue;
                for (int ci = 0; ci < dim_; ci++)
                  for (int cj = 0; cj < dim_; cj++)
                    {
                      rows.Append(dim_*li+ci);
                      cols.Append(dim_*lj+cj);
                      vals.Append(sum_elmat(dim_*i+ci, dim_*j+cj));
                    }
              }
          }
      }

    auto result = make_shared<LocalMatrix>();
    result->mat = SparseMatrix<double>::CreateFromCOO(rows, cols, vals,
                                                      local_size_, local_size_);
    detail::FillPatchMetadata(*result, support_elements_, support_elements_,
                              core_dofs_, support_dofs_, core_in_support_);

    return result;
  }


  shared_ptr<LocalVector>
  MyAssembleLocalNonlinearResidual(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> local_dofs)
  {
    LocalNonlinearOperator op(fes, a, std::move(local_dofs));
    return op.Residual(u);
  }


  shared_ptr<LocalMatrix>
  MyAssembleLocalNonlinearJacobian(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> local_dofs)
  {
    LocalNonlinearOperator op(fes, a, std::move(local_dofs));
    return op.Jacobian(u);
  }
}
