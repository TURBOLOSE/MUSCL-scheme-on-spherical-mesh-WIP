#pragma once

#include "MUSCL_geometry.hpp"

class MUSCL_base : MUSCL_base_geometry
{
protected:
    std::vector<std::array<double, 4>> U;
    std::vector<std::vector<double>> c, nu_plus, nu_minus;     // c_ij, nu^plus_ij nu^minus_ij
    std::vector<std::vector<std::vector<double>>> lamb, omega; // lambda_ijk, omega_ijk

    double dt, dx, C, gam, t;

public:
    MUSCL_base(pmp::SurfaceMesh mesh, std::vector<std::array<double, 4>> U_in, double dt, double dx, double C, double gam, double t) : MUSCL_base_geometry(mesh), U(U_in)
    {
        if (U_in.size() != mesh.n_faces())
        {
            std::cout << "U_in is wrong size, U_in.size() should be equal to mesh.n_faces()" << std::endl;
        }

        c.resize(this->n_faces());
        nu_plus.resize(this->n_faces());
        nu_minus.resize(this->n_faces());
        lamb.resize(this->n_faces());
        omega.resize(this->n_faces());

        for (size_t i; i < this->n_faces(); i++)//double check if needed
        {
            c[i].resize(neighbors[i].size());
            nu_plus[i].resize(neighbors[i].size());
            nu_minus[i].resize(neighbors[i].size());
            lamb[i].resize(neighbors[i].size());
            omega[i].resize(neighbors[i].size());
        }




    };
};