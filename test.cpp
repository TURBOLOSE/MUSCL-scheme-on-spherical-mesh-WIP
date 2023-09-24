#include <iostream>
#include "pmp/surface_mesh.h"
#include "sph_gen.h"

using namespace pmp;

int main()
{

    double a=1;
    size_t N=3;

    SurfaceMesh mesh = sqare(a, N);
    auto points = mesh.get_vertex_property<Point>("v:point");

    std::cout << "vertices: " << mesh.n_vertices() << std::endl;
    std::cout << "edges: " << mesh.n_edges() << std::endl;
    std::cout << "faces: " << mesh.n_faces() << std::endl;


    for (auto e : mesh.vertices())
    {

        std::cout<<points[e]<<std::endl;
    }
    



    return 0;
}