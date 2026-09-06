//#define CGAL_TETRAHEDRAL_REMESHING_DEBUG
//#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
//#define CGAL_DUMP_REMESHING_STEPS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <iostream>
#include <fstream>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/IO/File_medit.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;

  Remeshing_triangulation tr;

  const char* mesh = (argc > 2) ? argv[2] : "data/sphere.mesh";
  std::ifstream in(mesh);
  if (CGAL::SMDS_3::build_triangulation_from_file(in, tr))
    std::cout << "build triangulation ok" << std::endl;

  std::size_t n_complex = 0;
  for (auto e = tr.finite_edges_begin(); e != tr.finite_edges_end(); ++e) (void)e;
  std::cout << "input vertices: " << tr.number_of_vertices() << std::endl;

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length,
         CGAL::parameters::smooth_constrained_edges(true));
  (void)n_complex;

  std::ofstream ofs("dump.mesh");
  ofs.precision(17);
  CGAL::IO::write_MEDIT(ofs, tr);
  ofs.close();
  return EXIT_SUCCESS;
}
