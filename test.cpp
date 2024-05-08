#include <iostream>
#include "HLLE.hpp"
#include "HLLE_p.hpp"
#include "HLLC.hpp"
#include "HLLCplus.hpp"

using namespace pmp;

int main()
{
    //SurfaceMesh mesh = uv_sphere(40,40);
    //SurfaceMesh mesh = quad_sphere(0);
    //SurfaceMesh mesh = icosphere(3);
    SurfaceMesh mesh = icosphere_hex(4);

    //MUSCL_base_geometry test(mesh);

    double dt = 0.002;
    size_t maxstep = 10000+1;
    size_t skipstep = 1000;


    int dim = 5;
    double gam = 4./3;
    std::ifstream inData("input/input.dat");
    std::vector<std::vector<double>> U_in;
    U_in.resize(mesh.n_faces());



    for (size_t i = 0; i < mesh.n_faces(); i++)
    {
        U_in[i].resize(dim);
    }

    double element;
    std::vector<double> temp;

    int elements_read = 0;
    while (!inData.eof() && inData >> element)
    {
        temp.push_back(element);
        elements_read++;
    }

    if (elements_read < mesh.n_faces() * dim)
    {
        for (size_t i = elements_read; i < mesh.n_faces() * dim; i++)
        {
            temp.push_back(1);
        }

        std::cout << "input file does not have enough values!" << std::endl;
    }


    for (size_t i = 0; i < mesh.n_faces(); i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            U_in[i][j] = temp[i * dim + j];
        }

    }




    //MUSCL_HLLE test2(mesh, U_in, dim, gam);
    //MUSCL_HLLE_p test2(mesh, U_in, dim, gam);
    //MUSCL_HLLE_p_tracer test2(mesh, U_in, dim, gam);

    //MUSCL_HLLE_p test2(mesh, U_in, dim, gam);
    MUSCL_HLLC test2(mesh, U_in, dim, gam);
    //MUSCL_HLLCplus test2(mesh, U_in, dim, gam);
    //MUSCL_HLLC_tracer test2(mesh, U_in, dim, gam);


    // MUSCL_base_geometry test(mesh);

    test2.write_face_centers();
    test2.write_faces();
    test2.write_vertices();

    test2.write_t_rho();
    test2.write_t_p();
    test2.write_t_curl();
    test2.write_t_omega_z();
    //test2.write_t_tracer();

    for (size_t i = 0; i < maxstep; i++)
    {
        test2.do_step(dt);

        if (i % skipstep == 0)
        {
            test2.write_t_rho();
            test2.write_t_p();
            test2.write_t_curl();
            test2.write_t_omega_z();
            //test2.write_t_tracer();
        }

        if (test2.get_stop_check())
            break;

    }
    

    return 0;
}