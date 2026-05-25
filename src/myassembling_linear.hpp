#ifndef FILE_MYASSEMBLING_LINEAR_HPP
#define FILE_MYASSEMBLING_LINEAR_HPP


#include "myassembling_common.hpp"


/*
  Assembling matrix and vector
*/


namespace ngcomp
{
  shared_ptr<BaseSparseMatrix>
  MyAssembleMatrix(shared_ptr<FESpace> fes,
                   shared_ptr<BilinearFormIntegrator> bfi);

  /*
    Low-level assembly kernel for a supplied local support patch.

    The caller provides the support closure and all index maps:
    core_elements, support_elements, core_dofs, support_dofs, and
    core_in_support. This C++ routine does not decide what the support patch is;
    it assembles element matrices over support_elements and returns the
    core-core matrix in the local numbering defined by core_dofs.
  */
  shared_ptr<LocalMatrix>
  MyAssembleLocalMatrix(shared_ptr<FESpace> fes,
                        shared_ptr<BilinearFormIntegrator> bfi,
                        vector<int> core_elements,
                        vector<int> support_elements,
                        vector<int> core_dofs,
                        vector<int> support_dofs,
                        vector<int> core_in_support);

  /*
    Low-level vector assembly over support_elements, returning only the core
    entries in the local numbering defined by core_dofs.
  */
  shared_ptr<LocalVector>
  MyAssembleLocalVector(shared_ptr<FESpace> fes,
                        shared_ptr<LinearFormIntegrator> lfi,
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
