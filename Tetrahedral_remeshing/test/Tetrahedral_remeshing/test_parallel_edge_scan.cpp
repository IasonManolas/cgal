// Does parallel_collect_finite_edges() visit exactly the finite edges, each
// once?  The helper underpins the parallel candidate-collection arms, so if it
// misses or duplicates an edge those arms do different work from the serial
// ones and the comparison between them means nothing.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, CGAL::Parallel_tag> Tr;

int main()
{
  Tr tr;
  std::ifstream in("data/sphere.mesh");
  if (!CGAL::SMDS_3::build_triangulation_from_file(in, tr))
  {
    std::cerr << "could not read data/sphere.mesh\n";
    return EXIT_FAILURE;
  }

  using Vh = Tr::Vertex_handle;
  using Vpair = std::pair<Vh, Vh>;

  std::vector<Vpair> serial;
  for (const Tr::Edge& e : tr.finite_edges())
    serial.push_back(CGAL::Tetrahedral_remeshing::make_vertex_pair(e));

  std::vector<Vpair> parallel =
    CGAL::Tetrahedral_remeshing::internal::parallel_collect_finite_edges<Vpair>(
      tr, [](const Tr::Edge& e, std::vector<Vpair>& out)
          { out.push_back(CGAL::Tetrahedral_remeshing::make_vertex_pair(e)); });

  std::sort(serial.begin(), serial.end());
  std::sort(parallel.begin(), parallel.end());

  const bool same_size = serial.size() == parallel.size();
  const bool no_dups = std::adjacent_find(parallel.begin(), parallel.end())
                       == parallel.end();
  const bool same_set = (serial == parallel);

  std::cout << "finite edges, serial walk : " << serial.size() << "\n"
            << "finite edges, parallel scan: " << parallel.size() << "\n"
            << "same count : " << (same_size ? "yes" : "NO") << "\n"
            << "no duplicates in parallel  : " << (no_dups ? "yes" : "NO") << "\n"
            << "identical sets             : " << (same_set ? "yes" : "NO") << "\n";

  return (same_size && no_dups && same_set) ? EXIT_SUCCESS : EXIT_FAILURE;
}
