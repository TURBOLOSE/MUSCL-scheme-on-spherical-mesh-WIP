#include <iostream>
#include "MUSCL_base.hpp"

using namespace pmp;

int main()
{

    double a = 1;
    size_t N = 4; // n of vertices in col/row (total = NxN)

    //SurfaceMesh mesh = square(a, N);
    SurfaceMesh mesh = quad_sphere(0);

    auto points = mesh.get_vertex_property<Point>("v:point");

    MUSCL_base test(mesh);
    test.print_vertices();
    //test.print_neighbors();


    std::vector<double> r {0.3,0.3,sqrt(3)/3};
    std::cout<<test.point_in_face(r)<<std::endl;


    return 0;
}