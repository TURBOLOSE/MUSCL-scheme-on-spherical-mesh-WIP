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
    //test.print_vertices();
    //test.print_neighbors();


    std::vector<double> r1 {sqrt(3)/3,0,0};
    std::vector<double> r2 {0,0,sqrt(3)/3};

    std::cout<<test.broken_distance(r1,r2)<<std::endl;

    //std::cout<<test.point_in_face(r)<<std::endl;


    return 0;
}