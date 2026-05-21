
/*
  
  Assemble the system matrix

  The input is
  - a finite element space, which provides the basis functions
  - an integrator, which computes the element matrices
  
  The result is a sparse matrix
*/

#include "myassembling.hpp"
#include <map>
#include <set>

namespace ngcomp
{
  namespace
  {
    vector<int> SortedVector(const set<int> & values)
    {
      return vector<int> (values.begin(), values.end());
    }

    void CheckElementNumber(shared_ptr<MeshAccess> ma, int elnr)
    {
      if (elnr < 0 || elnr >= ma->GetNE(VOL))
        throw Exception("Volume element number out of range: " + ToString(elnr));
    }

    void AddDofsOfElements(shared_ptr<FESpace> fes,
                           const set<int> & elements,
                           set<int> & dofs)
    {
      Array<int> dnums;
      for (int elnr : elements)
        {
          fes->GetDofNrs(ElementId(VOL, elnr), dnums);
          for (int d : dnums)
            if (d >= 0)
              dofs.insert(d);
        }
    }

    void BuildDofToElementTable(shared_ptr<FESpace> fes,
                                vector<vector<int>> & dof_to_elements)
    {
      auto ma = fes->GetMeshAccess();
      dof_to_elements.assign(fes->GetNDof(), vector<int>());

      Array<int> dnums;
      for (int elnr = 0; elnr < ma->GetNE(VOL); elnr++)
        {
          fes->GetDofNrs(ElementId(VOL, elnr), dnums);
          for (int d : dnums)
            if (d >= 0)
              dof_to_elements[d].push_back(elnr);
        }
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
  MyBuildLocalSupportPatch(shared_ptr<FESpace> fes,
                           vector<int> core_elements)
  {
    auto ma = fes->GetMeshAccess();

    set<int> core_element_set;
    for (int elnr : core_elements)
      {
        CheckElementNumber(ma, elnr);
        core_element_set.insert(elnr);
      }

    set<int> core_dof_set;
    AddDofsOfElements(fes, core_element_set, core_dof_set);

    vector<vector<int>> dof_to_elements;
    BuildDofToElementTable(fes, dof_to_elements);

    set<int> support_element_set;
    for (int d : core_dof_set)
      for (int elnr : dof_to_elements[d])
        support_element_set.insert(elnr);

    set<int> support_dof_set;
    AddDofsOfElements(fes, support_element_set, support_dof_set);

    auto core_dofs = SortedVector(core_dof_set);
    auto support_dofs = SortedVector(support_dof_set);

    map<int, int> global_to_support;
    for (size_t i = 0; i < support_dofs.size(); i++)
      global_to_support[support_dofs[i]] = int(i);

    vector<int> core_in_support;
    core_in_support.reserve(core_dofs.size());
    for (int d : core_dofs)
      core_in_support.push_back(global_to_support.at(d));

    auto result = make_shared<LocalSupportMatrix>();
    result->core_elements = SortedVector(core_element_set);
    result->support_elements = SortedVector(support_element_set);
    result->core_dofs = std::move(core_dofs);
    result->support_dofs = std::move(support_dofs);
    result->core_in_support = std::move(core_in_support);

    return result;
  }


  shared_ptr<LocalSupportMatrix>
  MyAssembleLocalSupportMatrix(shared_ptr<FESpace> fes,
                               shared_ptr<BilinearFormIntegrator> bfi,
                               vector<int> core_elements)
  {
    auto support = MyBuildLocalSupportPatch(fes, std::move(core_elements));
    auto ma = fes->GetMeshAccess();

    map<int, int> global_to_support;
    for (size_t i = 0; i < support->support_dofs.size(); i++)
      global_to_support[support->support_dofs[i]] = int(i);

    Array<int> rows, cols;
    Array<double> vals;
    Array<int> dnums, ldnums;
    LocalHeap lh(1000*1000);

    /*
      Support-closure local assembly:

        A_support = sum_{K in support_elements} A_K

      where support_elements are exactly all elements touching core_dofs.
      This is designed so the core-core block of A_support matches the
      core-core block extracted from a globally assembled matrix.
    */
    for (int elnr : support->support_elements)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        FiniteElement & fel = fes->GetFE(ei, lh);
        fes->GetDofNrs(ei, dnums);

        ldnums.SetSize(dnums.Size());
        for (int i = 0; i < dnums.Size(); i++)
          {
            auto it = global_to_support.find(dnums[i]);
            ldnums[i] = (it == global_to_support.end()) ? -1 : it->second;
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

    support->mat = SparseMatrix<double>::CreateFromCOO(rows, cols, vals,
                                                       support->support_dofs.size(),
                                                       support->support_dofs.size());

    return support;
  }


  vector<shared_ptr<LocalSupportMatrix>>
  MyAssembleLocalSupportMatrices(shared_ptr<FESpace> fes,
                                 shared_ptr<BilinearFormIntegrator> bfi,
                                 vector<vector<int>> partition)
  {
    vector<shared_ptr<LocalSupportMatrix>> mats;
    mats.reserve(partition.size());

    for (auto & core_elements : partition)
      mats.push_back(MyAssembleLocalSupportMatrix(fes, bfi, core_elements));

    return mats;
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
