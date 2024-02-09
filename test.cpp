#include <iostream>
#include "HLLE.hpp"
#include "HLLE_p.hpp"
#include "HLLC.hpp"

using namespace pmp;

int main()
{

    //SurfaceMesh mesh = quad_sphere(0);
    SurfaceMesh mesh = icosphere(4);
    //SurfaceMesh mesh = icosphere_hex(4);

    double dt = 0.002;
    size_t maxstep = 10000;
    int dim = 4;
    double gam = 1.4;
    std::ifstream inData("input/input.dat");
    std::vector<std::vector<double>> U_in;
    U_in.resize(mesh.n_faces());
    vector3d<double> vel, L;

    for (size_t i = 0; i < mesh.n_faces(); i++)
    {
        U_in[i].resize(dim);
    }

    double element;
    std::vector<double> temp;

    while (!inData.eof() && inData >> element)
    {
        temp.push_back(element);
    }

    for (size_t i = 0; i < mesh.n_faces(); i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            U_in[i][j] = temp[i * dim + j];
        }
    }

    MUSCL_HLLE test2(mesh, U_in, dim, gam);
    
    //MUSCL_HLLE_p test2(mesh, U_in, dim, gam);

    // MUSCL_HLLC test2(mesh, U_in, dim, gam);

    // MUSCL_base_geometry test(mesh);

    test2.write_face_centers();
    test2.write_faces();
    test2.write_vertices();
    test2.write_t_rho();
    test2.write_t_curl();
    //test2.write_t_p();


   for (size_t i = 0; i < maxstep; i++)
    {
        test2.do_step(dt);
        test2.write_t_rho();
        test2.write_t_curl();
        //test2.write_t_p();

    }

    return 0;
}