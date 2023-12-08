#include <iostream>
#include "HLLE.hpp"

using namespace pmp;

int main()
{

     //SurfaceMesh mesh = quad_sphere(4);
     SurfaceMesh mesh = icosphere(3);
    //SurfaceMesh mesh = icosphere_hex(3);

    double dt = 0.002;
    size_t maxstep = 600;

    int dim = 4;
    double gam = 1;
    std::vector<std::vector<double>> U_in;
    U_in.resize(mesh.n_faces());
    vector3d<double> vel, L;

    // MUSCL_base_geometry test(mesh);

    for (size_t i = 0; i < mesh.n_faces(); i++)
    {
        U_in[i].resize(dim);
    }




    double element;
    std::vector<double> temp;
    std::ifstream inData("input/input.dat");
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

    test2.write_face_centers();
    test2.write_faces();
    test2.write_vertices();
    test2.write_t_rho();
    test2.write_t_curl();



    for (size_t i = 0; i < maxstep; i++)
    {
        test2.do_step(dt);
        test2.write_t_rho();
        test2.write_t_curl();
    }

    // test.print_neighbors();

    // std::vector<double> r1 {0,1./3-1e-3,-0.00333333};

    return 0;
}