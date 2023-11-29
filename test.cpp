#include <iostream>
#include "HLLE.hpp"

using namespace pmp;

int main()
{

    SurfaceMesh mesh = quad_sphere(2);
    // SurfaceMesh mesh = icosphere(2);
    //  SurfaceMesh mesh = icosphere_hex(1);

    int dim = 4;
    double gam = 1;
    std::vector<std::vector<double>> U_in;
    U_in.resize(mesh.n_faces());
    vector3d<double> vel, L;

    for (size_t i = 0; i < mesh.n_faces(); i++)
    {
        U_in[i].resize(4);
        U_in[i][0] = 1;
        U_in[i][1] = 0;
        U_in[i][2] = 0;
        U_in[i][3] = 1;
        
    }


    U_in[1][0]=2;

    MUSCL_HLLE test2(mesh, U_in, dim, gam);

    


    test2.write_neighbors();
    test2.write_t_rho();


    

    std::vector<double> u0;
    u0.resize(4); u0[0]=2;u0[1]=0; u0[2]=0; u0[3]=1;

    size_t maxstep=5;
    double dt=0.01;

    for (size_t i = 0; i < maxstep; i++)
    {
        test2.do_step(dt);
        test2.write_t_rho();
    }

    // test.print_neighbors();

    // std::vector<double> r1 {0,1./3-1e-3,-0.00333333};

    return 0;
}