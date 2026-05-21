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
    Low-level assembly kernel for a supplied local support patch.

    The caller provides the support closure and all index maps:
    core_elements, support_elements, core_dofs, support_dofs, and
    core_in_support. This C++ routine does not decide what the support patch is;
    it only assembles element matrices over support_elements and writes them
    into the local numbering defined by support_dofs.
  */
  shared_ptr<LocalSupportMatrix>
  MyAssembleGivenLocalSupportMatrix(shared_ptr<FESpace> fes,
                                    shared_ptr<BilinearFormIntegrator> bfi,
                                    vector<int> core_elements,
                                    vector<int> support_elements,
                                    vector<int> core_dofs,
                                    vector<int> support_dofs,
                                    vector<int> core_in_support);
    
  shared_ptr<BaseVector>
  MyAssembleVector(shared_ptr<FESpace> fes,
                   shared_ptr<LinearFormIntegrator> lfi);

  
  
}

#endif
