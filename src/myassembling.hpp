#ifndef FILE_MYASSEMBLING_HPP
#define FILE_MYASSEMBLING_HPP


#include <comp.hpp>
#include <vector>


/*
  Assembling matrix and vector
*/


namespace ngcomp
{
  struct LocalSupportMatrix
  {
    shared_ptr<BaseSparseMatrix> mat;
    vector<int> core_elements;
    vector<int> support_elements;
    vector<int> core_dofs;
    vector<int> support_dofs;
    vector<int> core_in_support;
  };
  
  shared_ptr<BaseSparseMatrix>
  MyAssembleMatrix(shared_ptr<FESpace> fes,
                   shared_ptr<BilinearFormIntegrator> bfi);

  /*
    Build the exact support closure of a local cell partition.

    core_dofs are all global dofs appearing on core_elements. support_elements
    are all volume elements whose dof list intersects core_dofs. Assembling on
    this support closure exactly reproduces the global core-core matrix block:

      A_support[core_in_support, core_in_support]
        == A_global[core_dofs, core_dofs]

    up to roundoff. The full support matrix is not expected to equal the
    algebraic global submatrix on support_dofs.
  */
  shared_ptr<LocalSupportMatrix>
  MyBuildLocalSupportPatch(shared_ptr<FESpace> fes,
                           vector<int> core_elements);

  shared_ptr<LocalSupportMatrix>
  MyAssembleLocalSupportMatrix(shared_ptr<FESpace> fes,
                               shared_ptr<BilinearFormIntegrator> bfi,
                               vector<int> core_elements);

  vector<shared_ptr<LocalSupportMatrix>>
  MyAssembleLocalSupportMatrices(shared_ptr<FESpace> fes,
                                 shared_ptr<BilinearFormIntegrator> bfi,
                                 vector<vector<int>> partition);
    
  shared_ptr<BaseVector>
  MyAssembleVector(shared_ptr<FESpace> fes,
                   shared_ptr<LinearFormIntegrator> lfi);

  
  
}

#endif
