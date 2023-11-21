#include <iostream>
#include "MUSCL_geometry.hpp"
#include "MUSCL_base.hpp"

using namespace pmp;

int main()
{

    SurfaceMesh mesh = quad_sphere(1);
    // SurfaceMesh mesh = icosphere(2);
    //  SurfaceMesh mesh = icosphere_hex(1);

     MUSCL_base_geometry test(mesh);


    /*int dim = 4;
    std::vector<std::vector<double>> U_in;
    U_in.resize(mesh.n_faces());
    for (size_t i=0; i<mesh.n_faces(); i++)
    {
        U_in[i].push_back(1);
        U_in[i].push_back(1);
        U_in[i].push_back(1);
        U_in[i].push_back(1);
    }

    MUSCL_base test2(mesh, U_in, dim, 1, 5. / 3);*/


    // test.print_neighbors();

    // std::vector<double> r1 {0,1./3-1e-3,-0.00333333};

    return 0;
}