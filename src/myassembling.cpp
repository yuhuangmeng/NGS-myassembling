
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
    string OverlapModeName(OverlapMode overlap_mode)
    {
      switch (overlap_mode)
        {
        case OverlapMode::FACET:
          return "facet";
        case OverlapMode::VERTEX:
          return "vertex";
        default:
          throw Exception("Unknown overlap mode");
        }
    }

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

    /*
      Build an element patch from non-overlapping core cells.

      By default, overlap_layers counts layers in the cell adjacency graph
      induced by shared facets. In a 2D triangular mesh, facets are edges,
      so this is the usual edge-neighbor cell overlap used for element-based
      overlapping Schwarz subdomains from a METIS cell partition.

      The vertex mode is an optional vertex-patch expansion. It is wider in
      general and is not the default domain-decomposition overlap used here.
    */
    set<int> ExtendElements(shared_ptr<MeshAccess> ma,
                            const set<int> & core_elements,
                            int overlap_layers,
                            OverlapMode overlap_mode)
    {
      if (overlap_layers < 0)
        throw Exception("overlap_layers must be non-negative");

      set<int> patch_elements = core_elements;
      Array<int> facet_neighbors;

      for (int layer = 0; layer < overlap_layers; layer++)
        {
          set<int> new_patch_elements = patch_elements;
          for (int elnr : patch_elements)
            {
              ElementId ei(VOL, elnr);

              if (overlap_mode == OverlapMode::VERTEX)
                {
                  for (auto v : ma->GetElement(ei).Vertices())
                    {
                      auto vertex_neighbors = ma->GetVertexElements(size_t(v), VOL);
                      for (auto neighbor : vertex_neighbors)
                        new_patch_elements.insert(int(neighbor));
                    }
                }
              else
                {
                  for (auto facet : ma->GetElFacets(ei))
                    {
                      ma->GetFacetElements(size_t(facet), facet_neighbors);
                      for (int neighbor : facet_neighbors)
                        new_patch_elements.insert(neighbor);
                    }
                }
            }
          patch_elements.swap(new_patch_elements);
        }

      return patch_elements;
    }

    void PrintPatchSummary(const LocalPatchMatrix & patch)
    {
      cout << "LocalPatchMatrix(mode=" << patch.overlap_mode
           << ", layers=" << patch.overlap_layers
           << ", core_elements=" << patch.core_elements.size()
           << ", overlap_elements=" << patch.overlap_elements.size()
           << ", core_dofs=" << patch.core_dofs.size()
           << ", overlap_dofs=" << patch.overlap_dofs.size()
           << ")" << endl;
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


  shared_ptr<LocalPatchMatrix>
  MyAssembleLocalPatchMatrix(shared_ptr<FESpace> fes,
                             shared_ptr<BilinearFormIntegrator> bfi,
                             vector<int> core_elements,
                             int overlap_layers,
                             OverlapMode overlap_mode,
                             bool verbose)
  {
    auto ma = fes->GetMeshAccess();

    set<int> core_set;
    for (int elnr : core_elements)
      {
        CheckElementNumber(ma, elnr);
        core_set.insert(elnr);
      }

    auto overlap_set = ExtendElements(ma, core_set, overlap_layers, overlap_mode);
    auto result = MyAssembleGivenLocalPatchMatrix(fes, bfi,
                                                  SortedVector(core_set),
                                                  SortedVector(overlap_set),
                                                  false);
    result->overlap_mode = OverlapModeName(overlap_mode);
    result->overlap_layers = overlap_layers;
    if (verbose)
      PrintPatchSummary(*result);
    return result;
  }


  shared_ptr<LocalPatchMatrix>
  MyAssembleGivenLocalPatchMatrix(shared_ptr<FESpace> fes,
                                  shared_ptr<BilinearFormIntegrator> bfi,
                                  vector<int> core_elements,
                                  vector<int> overlap_elements,
                                  bool verbose)
  {
    auto ma = fes->GetMeshAccess();

    set<int> core_set;
    set<int> overlap_set;

    for (int elnr : core_elements)
      {
        CheckElementNumber(ma, elnr);
        core_set.insert(elnr);
        overlap_set.insert(elnr);
      }

    for (int elnr : overlap_elements)
      {
        CheckElementNumber(ma, elnr);
        overlap_set.insert(elnr);
      }

    set<int> core_dof_set;
    set<int> overlap_dof_set;
    AddDofsOfElements(fes, core_set, core_dof_set);
    AddDofsOfElements(fes, overlap_set, overlap_dof_set);

    auto core_dofs = SortedVector(core_dof_set);
    auto overlap_dofs = SortedVector(overlap_dof_set);

    map<int, int> glob2loc;
    for (size_t i = 0; i < overlap_dofs.size(); i++)
      glob2loc[overlap_dofs[i]] = int(i);

    vector<int> core_local_dofs;
    core_local_dofs.reserve(core_dofs.size());
    for (int d : core_dofs)
      core_local_dofs.push_back(glob2loc.at(d));

    Array<int> rows, cols;
    Array<double> vals;
    Array<int> dnums, ldnums;
    LocalHeap lh(1000*1000);

    /*
      This local matrix is assembled element-by-element on overlap_set:

        A_i = sum_{K in overlap_elements} A_K

      with global dof numbers compressed to the overlap-local numbering.
      It is not formed by taking an algebraic submatrix of an already
      assembled global matrix.
    */
    for (int elnr : overlap_set)
      {
        HeapReset hr(lh);
        ElementId ei(VOL, elnr);

        FiniteElement & fel = fes->GetFE(ei, lh);
        fes->GetDofNrs(ei, dnums);

        ldnums.SetSize(dnums.Size());
        for (int i = 0; i < dnums.Size(); i++)
          {
            auto it = glob2loc.find(dnums[i]);
            ldnums[i] = (it == glob2loc.end()) ? -1 : it->second;
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

    auto result = make_shared<LocalPatchMatrix>();
    result->mat = SparseMatrix<double>::CreateFromCOO(rows, cols, vals,
                                                      overlap_dofs.size(),
                                                      overlap_dofs.size());
    result->core_elements = SortedVector(core_set);
    result->overlap_elements = SortedVector(overlap_set);
    result->core_dofs = std::move(core_dofs);
    result->overlap_dofs = std::move(overlap_dofs);
    result->core_local_dofs = std::move(core_local_dofs);
    result->overlap_mode = "given";
    result->overlap_layers = 0;

    if (verbose)
      PrintPatchSummary(*result);

    return result;
  }


  vector<shared_ptr<LocalPatchMatrix>>
  MyAssembleLocalPatchMatrices(shared_ptr<FESpace> fes,
                               shared_ptr<BilinearFormIntegrator> bfi,
                               vector<vector<int>> partition,
                               int overlap_layers,
                               OverlapMode overlap_mode,
                               bool verbose)
  {
    vector<shared_ptr<LocalPatchMatrix>> mats;
    mats.reserve(partition.size());

    for (auto & core_elements : partition)
      mats.push_back(MyAssembleLocalPatchMatrix(fes, bfi, core_elements,
                                                overlap_layers,
                                                overlap_mode,
                                                verbose));

    return mats;
  }


  vector<shared_ptr<LocalPatchMatrix>>
  MyAssembleGivenLocalPatchMatrices(shared_ptr<FESpace> fes,
                                    shared_ptr<BilinearFormIntegrator> bfi,
                                    vector<vector<int>> core_partition,
                                    vector<vector<int>> overlap_partition,
                                    bool verbose)
  {
    if (core_partition.size() != overlap_partition.size())
      throw Exception("core_partition and overlap_partition must have the same length");

    vector<shared_ptr<LocalPatchMatrix>> mats;
    mats.reserve(core_partition.size());

    for (size_t i = 0; i < core_partition.size(); i++)
      mats.push_back(MyAssembleGivenLocalPatchMatrix(fes, bfi,
                                                     core_partition[i],
                                                     overlap_partition[i],
                                                     verbose));

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
