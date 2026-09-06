#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <iostream>
#include <fstream>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/IO/File_medit.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, CGAL::Parallel_tag>
        Remeshing_triangulation;

int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;

  Remeshing_triangulation tr;

  std::ifstream in("data/sphere.mesh");
  if (CGAL::SMDS_3::build_triangulation_from_file(in, tr))
    std::cout << "build triangulation ok" << std::endl;

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length);

  std::ofstream ofs("dump_parallel.mesh");
  ofs.precision(17);
  CGAL::IO::write_MEDIT(ofs, tr);
  ofs.close();
  return EXIT_SUCCESS;
}
