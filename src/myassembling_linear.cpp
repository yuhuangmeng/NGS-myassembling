
/*

  Assemble the system matrix

  The input is
  - a finite element space, which provides the basis functions
  - an integrator, which computes the element matrices

  The result is a sparse matrix
*/

#include "myassembling_linear.hpp"

namespace ngcomp
{
  namespace detail = ngcomp::myassembling_detail;

  shared_ptr<BaseSparseMatrix> MyAssembleMatrix(shared_ptr<FESpace> fes,
                                                shared_ptr<BilinearFormIntegrator> bfi)
  {
    cout << "We assemble matrix" << endl;

    auto ma = fes->GetMeshAccess();

    int ndof = fes->GetNDof();
    int ne = ma->GetNE(VOL);

    // we build a sparse matrix
    // the non-zero pattern is given by the connectivity pattern provided by the FESpace
    // el2dof[i] stores dofs of element i
    Table<int> el2dof = fes->CreateDofTable(VOL);

    // generate sparse matrix of size ndof x ndof
    // from element-to-dof table for rows and columns
    auto mat = make_shared<SparseMatrix<double>> (ndof, ndof, el2dof, el2dof, false);
    mat -> SetZero();

    LocalHeap lh(1000*1000); // reserve 1MB
    Array<int> dnums;

    // loop over all volume elements
    for (int i = 0; i < ne; i++)
      {
        HeapReset hr(lh);  // cleanup heap at end of scope
        ElementId ei(VOL, i);

        // let FESpace generate the finite element
        FiniteElement & fel = fes->GetFE (ei, lh);

        // the global dof numbers of the element
        fes->GetDofNrs (ei, dnums);

        // the mesh knows the geometry of the element
        const ElementTransformation & trafo = ma->GetTrafo (ei, lh);

        // compute the element matrix
        FlatMatrix<> elmat (fel.GetNDof(), lh);
        bfi->CalcElementMatrix (fel, trafo, elmat, lh);

        mat->AddElementMatrix (dnums, dnums, elmat);
      }

    return mat;
  }


  shared_ptr<LocalMatrix>
  MyAssembleLocalMatrix(shared_ptr<FESpace> fes,
                        shared_ptr<BilinearFormIntegrator> bfi,
                        vector<int> core_elements,
                        vector<int> support_elements,
                        vector<int> core_dofs,
                        vector<int> support_dofs,
                        vector<int> core_in_support)
  {
    auto ma = fes->GetMeshAccess();

    detail::CheckElementNumbers(ma, core_elements);
    detail::CheckElementNumbers(ma, support_elements);

    detail::CheckSupportMetadata(core_dofs, support_dofs, core_in_support);
    auto global_to_support = detail::BuildGlobalToLocal(fes, support_dofs);
    auto global_to_core = detail::BuildGlobalToLocal(fes, core_dofs);

    Array<int> rows, cols;
    Array<double> vals;
    Array<int> dnums, core_ldnums;
    LocalHeap lh(1000*1000);

    /*
      Low-level support-patch assembly kernel. The support patch itself is
      constructed in Python; this loop assembles element matrices over the
      supplied support_elements but scatters only core rows and columns into
      the local numbering defined by core_dofs.
    */
    for (int elnr : support_elements)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        FiniteElement & fel = fes->GetFE(ei, lh);
        fes->GetDofNrs(ei, dnums);

        core_ldnums.SetSize(dnums.Size());
        for (int i = 0; i < dnums.Size(); i++)
          {
            if (dnums[i] < 0)
              {
                core_ldnums[i] = -1;
                continue;
              }

            if (detail::LocalDof(dnums[i], global_to_support) < 0)
              throw Exception("Element dof is missing from supplied support_dofs");

            core_ldnums[i] = detail::LocalDof(dnums[i], global_to_core);
          }

        const ElementTransformation & trafo = ma->GetTrafo(ei, lh);

        FlatMatrix<> elmat(fel.GetNDof(), lh);
        bfi->CalcElementMatrix(fel, trafo, elmat, lh);

        for (int i = 0; i < core_ldnums.Size(); i++)
          if (core_ldnums[i] >= 0)
            for (int j = 0; j < core_ldnums.Size(); j++)
              if (core_ldnums[j] >= 0)
                {
                  rows.Append(core_ldnums[i]);
                  cols.Append(core_ldnums[j]);
                  vals.Append(elmat(i, j));
                }
      }

    auto result = make_shared<LocalMatrix>();
    result->mat = SparseMatrix<double>::CreateFromCOO(rows, cols, vals,
                                                      core_dofs.size(),
                                                      core_dofs.size());
    detail::FillPatchMetadata(*result, core_elements, support_elements,
                              core_dofs, support_dofs, core_in_support);

    return result;
  }


  shared_ptr<LocalVector>
  MyAssembleLocalVector(shared_ptr<FESpace> fes,
                        shared_ptr<LinearFormIntegrator> lfi,
                        vector<int> core_elements,
                        vector<int> support_elements,
                        vector<int> core_dofs,
                        vector<int> support_dofs,
                        vector<int> core_in_support)
  {
    auto ma = fes->GetMeshAccess();

    detail::CheckElementNumbers(ma, core_elements);
    detail::CheckElementNumbers(ma, support_elements);

    detail::CheckSupportMetadata(core_dofs, support_dofs, core_in_support);
    auto global_to_support = detail::BuildGlobalToLocal(fes, support_dofs);
    auto global_to_core = detail::BuildGlobalToLocal(fes, core_dofs);

    auto vec = make_shared<VVector<double>> (core_dofs.size());
    vec->SetScalar(0.0);
    auto local_vec = vec->FV<double>();

    Array<int> dnums;
    LocalHeap lh(1000*1000);

    /*
      Low-level support-patch vector assembly kernel. Python supplies the
      support patch; C++ only assembles element vectors over support_elements
      and scatters only core entries into the local numbering defined by
      core_dofs.
    */
    for (int elnr : support_elements)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        FiniteElement & fel = fes->GetFE(ei, lh);
        fes->GetDofNrs(ei, dnums);

        const ElementTransformation & trafo = ma->GetTrafo(ei, lh);

        FlatVector<> elvec(fel.GetNDof(), lh);
        lfi->CalcElementVector(fel, trafo, elvec, lh);

        for (int i = 0; i < dnums.Size(); i++)
          {
            if (dnums[i] < 0)
              continue;
            if (detail::LocalDof(dnums[i], global_to_support) < 0)
              throw Exception("Element dof is missing from supplied support_dofs");

            int core_ldof = detail::LocalDof(dnums[i], global_to_core);
            if (core_ldof < 0)
              continue;
            local_vec(core_ldof) += elvec(i);
          }
      }

    auto result = make_shared<LocalVector>();
    result->vec = vec;
    detail::FillPatchMetadata(*result, core_elements, support_elements,
                              core_dofs, support_dofs, core_in_support);

    return result;
  }

  /*
    Exercise: implement a corresponding function for assembling the right hand side vector
   */
  shared_ptr<BaseVector> MyAssembleVector(shared_ptr<FESpace> fes,
                                          shared_ptr<LinearFormIntegrator> lfi)
  {
    shared_ptr<BaseVector> vec = make_shared<VVector<double>> (fes->GetNDof());
    vec->SetScalar(0.0);
    auto global_vec = vec->FV<double>();

    auto ma = fes->GetMeshAccess();
    int ne = ma->GetNE(VOL);
    Array<int> dnums;
    LocalHeap lh(1000*1000);

    for (int elnr = 0; elnr < ne; elnr++)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        FiniteElement & fel = fes->GetFE(ei, lh);
        fes->GetDofNrs(ei, dnums);

        const ElementTransformation & trafo = ma->GetTrafo(ei, lh);

        FlatVector<> elvec(fel.GetNDof(), lh);
        lfi->CalcElementVector(fel, trafo, elvec, lh);

        for (int i = 0; i < dnums.Size(); i++)
          if (dnums[i] >= 0)
            global_vec(dnums[i]) += elvec(i);
      }

    return vec;
  }
}
