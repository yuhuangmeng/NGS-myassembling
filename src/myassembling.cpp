
/*
  
  Assemble the system matrix

  The input is
  - a finite element space, which provides the basis functions
  - an integrator, which computes the element matrices
  
  The result is a sparse matrix
*/

#include "myassembling.hpp"
#include <map>

namespace ngcomp
{
  namespace
  {
    void CheckElementNumber(shared_ptr<MeshAccess> ma, int elnr)
    {
      if (elnr < 0 || elnr >= ma->GetNE(VOL))
        throw Exception("Volume element number out of range: " + ToString(elnr));
    }
  }

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


  shared_ptr<LocalSupportMatrix>
  MyAssembleGivenLocalSupportMatrix(shared_ptr<FESpace> fes,
                                    shared_ptr<BilinearFormIntegrator> bfi,
                                    vector<int> core_elements,
                                    vector<int> support_elements,
                                    vector<int> core_dofs,
                                    vector<int> support_dofs,
                                    vector<int> core_in_support)
  {
    auto ma = fes->GetMeshAccess();

    for (int elnr : core_elements)
      CheckElementNumber(ma, elnr);
    for (int elnr : support_elements)
      CheckElementNumber(ma, elnr);

    if (core_in_support.size() != core_dofs.size())
      throw Exception("core_in_support and core_dofs must have the same length");

    map<int, int> global_to_support;
    for (size_t i = 0; i < support_dofs.size(); i++)
      {
        int d = support_dofs[i];
        if (d < 0 || d >= fes->GetNDof())
          throw Exception("support_dofs contains invalid dof: " + ToString(d));
        global_to_support[d] = int(i);
      }

    for (size_t i = 0; i < core_dofs.size(); i++)
      {
        int local = core_in_support[i];
        if (local < 0 || local >= support_dofs.size())
          throw Exception("core_in_support index out of range");
        if (support_dofs[local] != core_dofs[i])
          throw Exception("support_dofs[core_in_support[i]] must equal core_dofs[i]");
      }

    Array<int> rows, cols;
    Array<double> vals;
    Array<int> dnums, ldnums;
    LocalHeap lh(1000*1000);

    /*
      Low-level support-patch assembly kernel. The support patch itself is
      constructed in Python; this loop only assembles element matrices over the
      supplied support_elements using the supplied support_dofs numbering.
    */
    for (int elnr : support_elements)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        FiniteElement & fel = fes->GetFE(ei, lh);
        fes->GetDofNrs(ei, dnums);

        ldnums.SetSize(dnums.Size());
        for (int i = 0; i < dnums.Size(); i++)
          {
            auto it = global_to_support.find(dnums[i]);
            if (dnums[i] < 0)
              ldnums[i] = -1;
            else if (it == global_to_support.end())
              throw Exception("Element dof is missing from supplied support_dofs");
            else
              ldnums[i] = it->second;
          }

        const ElementTransformation & trafo = ma->GetTrafo(ei, lh);

        FlatMatrix<> elmat(fel.GetNDof(), lh);
        bfi->CalcElementMatrix(fel, trafo, elmat, lh);

        for (int i = 0; i < ldnums.Size(); i++)
          if (ldnums[i] >= 0)
            for (int j = 0; j < ldnums.Size(); j++)
              if (ldnums[j] >= 0)
                {
                  rows.Append(ldnums[i]);
                  cols.Append(ldnums[j]);
                  vals.Append(elmat(i, j));
                }
      }

    auto result = make_shared<LocalSupportMatrix>();
    result->mat = SparseMatrix<double>::CreateFromCOO(rows, cols, vals,
                                                      support_dofs.size(),
                                                      support_dofs.size());
    result->core_elements = std::move(core_elements);
    result->support_elements = std::move(support_elements);
    result->core_dofs = std::move(core_dofs);
    result->support_dofs = std::move(support_dofs);
    result->core_in_support = std::move(core_in_support);

    return result;
  }

  /*
    Exercise: implement a corresponding function for assembling the right hand side vector
   */
  shared_ptr<BaseVector> MyAssembleVector(shared_ptr<FESpace> fes,
                                          shared_ptr<LinearFormIntegrator> lfi)
  {
    // A VVector (virtual vector) is derived from BaseVector
    shared_ptr<BaseVector> vec = make_shared<VVector<double>> (fes->GetNDof());

    // adding an element vector to the global vector
    FlatVector<double> elvec;
    Array<int> dnums;
    vec->AddIndirect (dnums, elvec);
    
    return vec;
  }
}
