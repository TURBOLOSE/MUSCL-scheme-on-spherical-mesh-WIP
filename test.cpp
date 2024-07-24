#include <iostream>
#include "src/Riemann_solvers/HLLE.hpp"
#include "src/Riemann_solvers/HLLE_p.hpp"
#include "src/Riemann_solvers/HLLC.hpp"
#include "src/Riemann_solvers/HLLCplus.hpp"
#include "src/json.hpp"

using json = nlohmann::json;

using namespace pmp;

int main()
{
    //SurfaceMesh mesh = uv_sphere(50,50);
    //SurfaceMesh mesh = quad_sphere(5);
    //SurfaceMesh mesh = icosphere(4);
    SurfaceMesh mesh = icosphere_hex(5);

    //MUSCL_base_geometry test(mesh);

    //double dt = 0.002;

    
    std::ifstream ifs("input/parameters.json");
    json parameters = json::parse(ifs);

    double dt = parameters["dt"];
    size_t maxstep = parameters["maxstep"];
    size_t skipstep = parameters["skipstep"];

    size_t dim = parameters["dim"];
    double gam3d = parameters["gam3d"];
    double gam=2-1/gam3d;
    double omega_ns=parameters["omega_ns"];
    size_t threads=parameters["threads"];

    bool accretion_on=parameters["accretion_on"];

   
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

        maxstep=0;
        std::cout << "input file does not have enough values!" << std::endl;
    }


    for (size_t i = 0; i < mesh.n_faces(); i++)
    {
        for (size_t j = 0; j < dim; j++)
        {
            U_in[i][j] = temp[i * dim + j];
        }

    }

    //MUSCL_HLLE test2(mesh, U_in, dim, gam, threads);
    //MUSCL_HLLE_p test2(mesh, U_in, dim, gam,omega_ns, threads);
    //MUSCL_HLLC test2(mesh, U_in, dim, gam,omega_ns, threads);
    MUSCL_HLLCplus test2(mesh, U_in, dim, gam, omega_ns,accretion_on, threads);



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