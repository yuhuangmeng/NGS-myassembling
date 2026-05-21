#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/stl.h>

#include "myintegrator.hpp"
#include "myassembling.hpp"


namespace
{
  ngcomp::OverlapMode ParseOverlapMode(py::object overlap_mode,
                                       py::object vertex_overlap)
  {
    if (!overlap_mode.is_none())
      {
        if (py::isinstance<py::bool_>(overlap_mode))
          return overlap_mode.cast<bool>() ? ngcomp::OverlapMode::VERTEX
                                           : ngcomp::OverlapMode::FACET;

        string mode = overlap_mode.cast<string>();
        if (mode == "facet")
          return ngcomp::OverlapMode::FACET;
        if (mode == "vertex")
          return ngcomp::OverlapMode::VERTEX;

        throw Exception("overlap_mode must be 'facet' or 'vertex'");
      }

    if (!vertex_overlap.is_none() && vertex_overlap.cast<bool>())
      return ngcomp::OverlapMode::VERTEX;

    return ngcomp::OverlapMode::FACET;
  }
}


PYBIND11_MODULE(myassembling, m)
{
  cout << "Loading myassembling library" << endl;

  py::class_<ngfem::MyLaplaceIntegrator, shared_ptr<ngfem::MyLaplaceIntegrator>, ngfem::BilinearFormIntegrator>
    (m, "MyLaplace")
    .def(py::init<shared_ptr<ngfem::CoefficientFunction>>())
    ;
  
  py::class_<ngfem::MySourceIntegrator, shared_ptr<ngfem::MySourceIntegrator>, ngfem::LinearFormIntegrator>
    (m, "MySource")
    .def(py::init<shared_ptr<ngfem::CoefficientFunction>>())
    ;

  py::class_<ngcomp::LocalPatchMatrix, shared_ptr<ngcomp::LocalPatchMatrix>>
    (m, "LocalPatchMatrix")
    .def_readonly("mat", &ngcomp::LocalPatchMatrix::mat)
    .def_readonly("core_elements", &ngcomp::LocalPatchMatrix::core_elements)
    .def_readonly("overlap_elements", &ngcomp::LocalPatchMatrix::overlap_elements)
    .def_readonly("core_dofs", &ngcomp::LocalPatchMatrix::core_dofs)
    .def_readonly("overlap_dofs", &ngcomp::LocalPatchMatrix::overlap_dofs)
    .def_readonly("core_local_dofs", &ngcomp::LocalPatchMatrix::core_local_dofs)
    .def_readonly("overlap_mode", &ngcomp::LocalPatchMatrix::overlap_mode)
    .def_readonly("overlap_layers", &ngcomp::LocalPatchMatrix::overlap_layers)
    ;
  
  m.def ("MyAssembleMatrix",
         &ngcomp::MyAssembleMatrix,
         py::arg("fes"),py::arg("integrator"));

  m.def ("MyAssembleLocalPatchMatrix",
         [](shared_ptr<ngcomp::FESpace> fes,
            shared_ptr<ngfem::BilinearFormIntegrator> bfi,
            vector<int> core_elements,
            int overlap_layers,
            py::object overlap_mode,
            py::object vertex_overlap,
            bool verbose)
         {
           return ngcomp::MyAssembleLocalPatchMatrix(fes, bfi, core_elements,
                                                     overlap_layers,
                                                     ParseOverlapMode(overlap_mode,
                                                                      vertex_overlap),
                                                     verbose);
         },
         py::arg("fes"), py::arg("integrator"), py::arg("core_elements"),
         py::arg("overlap_layers") = 1,
         py::arg("overlap_mode") = py::none(),
         py::arg("vertex_overlap") = py::none(),
         py::arg("verbose") = false,
         R"raw_string(Assemble an element-based overlapping local matrix.

overlap_mode="facet" is the default and recommended mode. overlap_layers
counts layers in the cell adjacency graph induced by shared facets. In 2D
triangular meshes this is shared-edge expansion. overlap_mode="vertex" is an
optional wider vertex-patch mode. The legacy vertex_overlap keyword is still
accepted when overlap_mode is not provided.)raw_string");

  m.def ("MyAssembleGivenLocalPatchMatrix",
         &ngcomp::MyAssembleGivenLocalPatchMatrix,
         py::arg("fes"), py::arg("integrator"), py::arg("core_elements"),
         py::arg("overlap_elements"),
         py::arg("verbose") = false);

  m.def ("MyAssembleLocalPatchMatrices",
         [](shared_ptr<ngcomp::FESpace> fes,
            shared_ptr<ngfem::BilinearFormIntegrator> bfi,
            vector<vector<int>> partition,
            int overlap_layers,
            py::object overlap_mode,
            py::object vertex_overlap,
            bool verbose)
         {
           return ngcomp::MyAssembleLocalPatchMatrices(fes, bfi, partition,
                                                       overlap_layers,
                                                       ParseOverlapMode(overlap_mode,
                                                                        vertex_overlap),
                                                       verbose);
         },
         py::arg("fes"), py::arg("integrator"), py::arg("partition"),
         py::arg("overlap_layers") = 1,
         py::arg("overlap_mode") = py::none(),
         py::arg("vertex_overlap") = py::none(),
         py::arg("verbose") = false,
         R"raw_string(Assemble element-based overlapping local matrices for all core patches.

overlap_mode="facet" is the default and recommended mode. overlap_layers
counts layers in the cell adjacency graph induced by shared facets. In 2D
triangular meshes this is shared-edge expansion. overlap_mode="vertex" is an
optional wider vertex-patch mode. The legacy vertex_overlap keyword is still
accepted when overlap_mode is not provided.)raw_string");

  m.def ("MyAssembleGivenLocalPatchMatrices",
         &ngcomp::MyAssembleGivenLocalPatchMatrices,
         py::arg("fes"), py::arg("integrator"), py::arg("core_partition"),
         py::arg("overlap_partition"),
         py::arg("verbose") = false);
}    
