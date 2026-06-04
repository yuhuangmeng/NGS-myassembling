#include "myassembling_nonlinear.hpp"

#include <algorithm>

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

    void PrepareBoundaryData(shared_ptr<FESpace> fes,
                             std::vector<int> boundary_dofs,
                             std::vector<double> boundary_values,
                             std::vector<int> & sorted_dofs,
                             std::vector<double> & sorted_values,
                             std::vector<int> & global_to_boundary)
    {
      if (boundary_dofs.size() != boundary_values.size())
        throw Exception("AssembleNonlinearLocal: boundary_dofs and boundary_values must have the same length");

      std::vector<std::pair<int, double>> pairs;
      pairs.reserve(boundary_dofs.size());
      for (size_t i = 0; i < boundary_dofs.size(); i++)
        {
          int gdof = boundary_dofs[i];
          if (gdof < 0 || gdof >= fes->GetNDof())
            throw Exception("AssembleNonlinearLocal: boundary dof is outside the FESpace dof range: " + ToString(gdof));
          pairs.emplace_back(gdof, boundary_values[i]);
        }

      std::sort(pairs.begin(), pairs.end(),
                [] (const auto & a, const auto & b) { return a.first < b.first; });

      sorted_dofs.clear();
      sorted_values.clear();
      sorted_dofs.reserve(pairs.size());
      sorted_values.reserve(pairs.size());
      global_to_boundary.assign(fes->GetNDof(), -1);

      for (size_t i = 0; i < pairs.size(); i++)
        {
          if (i > 0 && pairs[i].first == pairs[i-1].first)
            throw Exception("AssembleNonlinearLocal: boundary_dofs contains duplicate global dof: " + ToString(pairs[i].first));
          sorted_dofs.push_back(pairs[i].first);
          sorted_values.push_back(pairs[i].second);
          global_to_boundary[pairs[i].first] = int(i);
        }
    }
  }

  LocalNonlinearOperator ::
  LocalNonlinearOperator(shared_ptr<FESpace> fes,
                         shared_ptr<BilinearForm> a,
                         std::vector<int> local_dofs,
                         std::vector<int> boundary_dofs,
                         std::vector<double> boundary_values)
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
    PrepareBoundaryData(fes_, std::move(boundary_dofs), std::move(boundary_values),
                        boundary_dofs_, boundary_values_, global_to_boundary_);

    dim_ = fes_->GetDimension();
    if (dim_ != 1 && !boundary_dofs_.empty())
      throw Exception("AssembleNonlinearLocal: Dirichlet row treatment is currently implemented only for scalar FESpaces");

    local_size_ = int(core_dofs_.size()) * dim_;
  }


  int LocalNonlinearOperator ::
  BoundaryIndex(int gdof) const
  {
    if (gdof < 0 || size_t(gdof) >= global_to_boundary_.size())
      return -1;
    return global_to_boundary_[gdof];
  }


  bool LocalNonlinearOperator ::
  IsBoundaryDof(int gdof) const
  {
    return BoundaryIndex(gdof) >= 0;
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

    if (!boundary_dofs_.empty())
      {
        Array<int> bdnum(1);
        FlatVector<double> bdval(1, lh);
        for (size_t i = 0; i < core_dofs_.size(); i++)
          {
            int bidx = BoundaryIndex(core_dofs_[i]);
            if (bidx < 0)
              continue;
            bdnum[0] = core_dofs_[i];
            u.GetIndirect(bdnum, bdval);
            local_fv(int(i)) = bdval(0) - boundary_values_[bidx];
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

        // Row-only Dirichlet treatment:
        // boundary residual rows are F_i(u) = u_i - g_i, hence the
        // corresponding Jacobian rows are identity rows. We do NOT clear
        // boundary columns, because interior residual rows are left unchanged
        // and may still depend on boundary dofs. Clearing columns would
        // require additional residual/RHS corrections for nonzero g_i.
        for (int i = 0; i < dnums.Size(); i++)
          {
            int li = detail::LocalDof(dnums[i], global_to_core_);
            if (li < 0)
              continue;
            if (IsBoundaryDof(dnums[i]))
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

    for (size_t i = 0; i < core_dofs_.size(); i++)
      if (IsBoundaryDof(core_dofs_[i]))
        {
          rows.Append(int(i));
          cols.Append(int(i));
          vals.Append(1.0);
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
                                   std::vector<int> local_dofs,
                                   std::vector<int> boundary_dofs,
                                   std::vector<double> boundary_values)
  {
    LocalNonlinearOperator op(fes, a, std::move(local_dofs),
                              std::move(boundary_dofs),
                              std::move(boundary_values));
    return op.Residual(u);
  }


  shared_ptr<LocalMatrix>
  MyAssembleLocalNonlinearJacobian(shared_ptr<FESpace> fes,
                                   shared_ptr<BilinearForm> a,
                                   const BaseVector & u,
                                   std::vector<int> local_dofs,
                                   std::vector<int> boundary_dofs,
                                   std::vector<double> boundary_values)
  {
    LocalNonlinearOperator op(fes, a, std::move(local_dofs),
                              std::move(boundary_dofs),
                              std::move(boundary_values));
    return op.Jacobian(u);
  }
}
