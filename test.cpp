#include <iostream>
#include "MUSCL_geometry.hpp"
#include "MUSCL_base.hpp"

using namespace pmp;

int main()
{


    SurfaceMesh mesh = quad_sphere(0);
    //SurfaceMesh mesh = icosphere(0);
    //SurfaceMesh mesh = icosphere_hex(3);


    MUSCL_base_geometry test(mesh);
    //test.print_vertices();
    

    std::vector<std::array<double,4>> U_in;
    U_in.resize(mesh.n_faces());


    //MUSCL_base test2(mesh, U_in,1.,1.,1.,1.,1.);
    


    //test.print_neighbors();


    //std::vector<double> r1 {0,1./3-1e-3,-0.00333333};
    vector3d<double> r1 {0.333333,0.666667-1e-3,-0.00745356};
    



    //std::cout<<test.broken_distance(r1,r2)<<std::endl;

    //std::cout<<test.point_in_face(r)<<std::endl;


    return 0;
}