#include <iostream>
#include "MUSCL_base.hpp"

using namespace pmp;

int main()
{

    double a = 1;
    size_t N = 4; // n of vertices in col/row (total = NxN)


    //SurfaceMesh mesh = quad_sphere(0);
    SurfaceMesh mesh = icosphere(0);
    //SurfaceMesh mesh = icosphere_hex(0);


    MUSCL_base_geometry test(mesh);
    
    

    //test.print_vertices();
    //test.print_neighbors();


    //std::vector<double> r1 {0,1./3-1e-3,-0.00333333};
    vector3d<double> r1 {0.333333,0.666667-1e-3,-0.00745356};
    



    //std::cout<<test.broken_distance(r1,r2)<<std::endl;

    //std::cout<<test.point_in_face(r)<<std::endl;


    return 0;
}