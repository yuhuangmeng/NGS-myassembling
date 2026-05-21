#ifndef FILE_MYASSEMBLING_HPP
#define FILE_MYASSEMBLING_HPP


#include <comp.hpp>
#include <vector>


/*
  Assembling matrix and vector
*/


namespace ngcomp
{
  enum class OverlapMode
  {
    FACET,
    VERTEX
  };

  struct LocalPatchMatrix
  {
    shared_ptr<BaseSparseMatrix> mat;
    vector<int> core_elements;
    vector<int> overlap_elements;
    vector<int> core_dofs;
    vector<int> overlap_dofs;
    vector<int> core_local_dofs;
    string overlap_mode;
    int overlap_layers = 0;
  };
  
  shared_ptr<BaseSparseMatrix>
  MyAssembleMatrix(shared_ptr<FESpace> fes,
                   shared_ptr<BilinearFormIntegrator> bfi);

  /*
    Assemble an element-based overlapping Schwarz local matrix.

    core_elements are the non-overlapping cell subdomain, for example from a
    METIS cell partition. overlap_layers counts layers in the cell adjacency
    graph induced by shared facets by default. In 2D triangular meshes, this is
    shared-edge expansion. The optional VERTEX mode builds a wider vertex-patch
    overlap and is not the recommended default.

    The matrix is reassembled over overlap_elements and compressed to
    overlap-local dof numbering. It is not an algebraic submatrix of a
    previously assembled global matrix.
  */
  shared_ptr<LocalPatchMatrix>
  MyAssembleLocalPatchMatrix(shared_ptr<FESpace> fes,
                             shared_ptr<BilinearFormIntegrator> bfi,
                             vector<int> core_elements,
                             int overlap_layers = 1,
                             OverlapMode overlap_mode = OverlapMode::FACET,
                             bool verbose = false);

  shared_ptr<LocalPatchMatrix>
  MyAssembleGivenLocalPatchMatrix(shared_ptr<FESpace> fes,
                                  shared_ptr<BilinearFormIntegrator> bfi,
                                  vector<int> core_elements,
                                  vector<int> overlap_elements,
                                  bool verbose = false);

  /*
    Assemble one local matrix per core subdomain. The default FACET mode is the
    standard cell-neighbor overlap for element-based domain decomposition.
  */
  vector<shared_ptr<LocalPatchMatrix>>
  MyAssembleLocalPatchMatrices(shared_ptr<FESpace> fes,
                               shared_ptr<BilinearFormIntegrator> bfi,
                               vector<vector<int>> partition,
                               int overlap_layers = 1,
                               OverlapMode overlap_mode = OverlapMode::FACET,
                               bool verbose = false);

  vector<shared_ptr<LocalPatchMatrix>>
  MyAssembleGivenLocalPatchMatrices(shared_ptr<FESpace> fes,
                                    shared_ptr<BilinearFormIntegrator> bfi,
                                    vector<vector<int>> core_partition,
                                    vector<vector<int>> overlap_partition,
                                    bool verbose = false);
    
  shared_ptr<BaseVector>
  MyAssembleVector(shared_ptr<FESpace> fes,
                   shared_ptr<LinearFormIntegrator> lfi);

  
  
}

#endif
