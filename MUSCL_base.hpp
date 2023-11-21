#pragma once

#include "MUSCL_geometry.hpp"

class MUSCL_base : MUSCL_base_geometry
{
protected:
    std::vector<std::vector<double>> U;
    std::vector<std::vector<std::vector<double>>> nu_plus, nu_minus, U_plus, U_minus; // nu^plus_ij nu^minus_ij
    double dt, C, gam, M, N, h0;
    int dim;

public:
    MUSCL_base(SurfaceMesh mesh, std::vector<std::vector<double>> U_in, int dim, double C, double gam) : MUSCL_base_geometry(mesh), U(U_in), gam(gam), C(C), dim(dim)
    {

        double h0_temp;

        if (U_in.size() != mesh.n_faces())
        {
            std::cout << "U_in is wrong size, U_in.size() should be equal to mesh.n_faces()" << std::endl;
        }

        // memory allocation
        nu_plus.resize(this->n_faces());
        nu_minus.resize(this->n_faces());
        U_minus.resize(this->n_faces());
        U_plus.resize(this->n_faces());

        for (size_t i = 0; i < this->n_faces(); i++)
        {

            nu_plus[i].resize(faces[i].size());
            nu_minus[i].resize(faces[i].size());
            U_minus[i].resize(faces[i].size());
            U_plus[i].resize(faces[i].size());

            for (size_t j = 0; j < faces[i].size(); j++)
            {
                nu_plus[i][j].resize(dim);
                nu_minus[i][j].resize(dim);
                U_minus[i][j].resize(dim);
                U_plus[i][j].resize(dim);
            }
        }

        N = 1;
        h0 = 10;
        for (size_t n_face = 0; n_face < this->n_faces(); n_face++)
        {
            if (faces[n_face].size() > N)
            {
                N = faces[n_face].size();
            }

            for (size_t face_element = 0; face_element < faces[n_face].size(); face_element++)
            {

                if (face_element < faces[n_face].size() - 1)
                {
                    h0_temp = (vertices[faces[n_face][face_element]] - vertices[faces[n_face][face_element + 1]]).norm() /
                              surface_area[n_face];
                }
                else
                {
                    h0_temp = (vertices[faces[n_face][face_element]] - vertices[faces[n_face][0]]).norm() /
                              surface_area[n_face];
                }

                if (h0_temp < h0)
                {
                    h0 = h0_temp;
                }
            }
        }

        M = 1; // temp
        find_M();
        dt = h0 / (2 * M * N); // first dt set



    };

    void do_step()
    {
    }

protected:
    double limiter(double r, double etha_plus, double etha_minus)
    { // classical Superbee limiter for irregular grids
        // CFL independent
        return std::max(0., std::max(std::min(1., etha_minus * r), std::min(r, etha_plus)));
    };

private:
    void find_M() // computes value M = max (d phi/ d U1, - d phi/ d U2 )
    {

        // std::cout<<U_H_plus(1,1)[0]<<" "<<U_H_plus(1,1)[1]<<" "<<U_H_plus(1,1)[2]<<" "<<U_H_plus(1,1)[3]<<std::endl;
    }

    void find_U_edges()
    {
        for (size_t i = 0; i < this->n_faces(); i++)
        {

            for (size_t j = 0; j < faces[i].size(); j++)
            {
                for (size_t k = 0; k < dim; k++)
                {
                }
            }
        }
    };

    std::vector<double> U_H_plus(int n_face, int face_edge)
    {

        std::vector<double> res;
        res.resize(dim);

        for (size_t U_element = 0; U_element < dim; U_element++)
        {
            res[U_element] = betas_plus[n_face][face_edge][0] * U[flux_faces_plus[n_face][face_edge][0]][U_element] +
                             betas_plus[n_face][face_edge][1] * U[flux_faces_plus[n_face][face_edge][1]][U_element];
        }

        return res;
    };

    std::vector<double> U_H_minus(int n_face, int face_edge)
    {
        std::vector<double> res;
        res.resize(dim);

        for (size_t U_element = 0; U_element < dim; U_element++)
        {
            res[U_element] = betas_minus[n_face][face_edge][0] * U[flux_faces_minus[n_face][face_edge][0]][U_element] +
                             betas_minus[n_face][face_edge][1] * U[flux_faces_minus[n_face][face_edge][1]][U_element];
        }

        return res;
    };

    std::vector<double> p_plus(int n_face, int face_edge)
    {
        std::vector<double> res, Up;
        res.resize(dim);
        Up = U_H_plus(n_face, face_edge);
        for (size_t U_element = 0; U_element < dim; U_element++)
        {
            res[U_element] = (Up[U_element] - U[n_face][U_element]) / H_plus[n_face][face_edge];
        }
        return res;
    }

    std::vector<double> p_minus(int n_face, int face_edge)
    {
        std::vector<double> res, Um;
        res.resize(dim);
        Um = U_H_minus(n_face, face_edge);
        for (size_t U_element = 0; U_element < dim; U_element++)
        {
            res[U_element] = (-Um[U_element] + U[n_face][U_element]) / H_minus[n_face][face_edge];
        }
        return res;
    }

    std::vector<double> flux_temp(std::vector<double> u0)
    {
        // placeholder, remove later
        return u0;
    }

    std::vector<double> HLLC_temp(std::vector<double> u0)
    {
        // placeholder, remove later
        return flux_temp(u0);
    }
};