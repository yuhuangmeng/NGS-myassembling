#ifndef FILE_MYASSEMBLING_NONLINEAR_HPP
#define FILE_MYASSEMBLING_NONLINEAR_HPP

#include <string>

namespace ngcomp
{
  /*
    Reserved for future true local-size nonlinear residual/Jacobian assembly.

    The current nonlinear workflow is intentionally implemented in the Python
    demo using NGSolve's automatic linearization:
      - BilinearForm::Apply
      - BilinearForm::AssembleLinearization
      - dx(definedonelements=...)

    That workflow verifies global-size residuals and Jacobians restricted to
    support elements. It does not yet implement a local support-dof nonlinear
    C++ kernel.
  */
  std::string NonlinearAssemblyStatus();
}

#endif
