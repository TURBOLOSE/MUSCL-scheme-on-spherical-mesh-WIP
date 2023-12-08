#pragma once
#include "MUSCL_base.hpp"
#include "MUSCL_geometry.hpp"

class MUSCL_HLLE : public MUSCL_base
{

private:
    std::ofstream outfile;

public:
    MUSCL_HLLE(SurfaceMesh mesh, std::vector<std::vector<double>> U_in, int dim, double gam)
        : MUSCL_base(mesh, U_in, dim, gam)
    {
        outfile.open("results/rho.dat", std::ios::out | std::ios::trunc);
        outfile.close();
        outfile.open("results/rho.dat", std::ios::out | std::ios::app);
    }

    void print_rho()
    {
        for (auto U_i : U)
        {
            std::cout << U_i[0] << std::endl;
        }
    };

    void write_t_rho()
    {
        outfile << this->time() << "  ";
        for (auto U_i : U)
        {

            outfile << U_i[0] << " ";
        }
        outfile << "\n";
    };

protected:
    const double a = 1;

    // U = {rho, l1, l2, l3}
    std::vector<double> flux(std::vector<double> u_in, int n_face, int n_edge)
    {

        std::vector<double> res;
        res.resize(4);
        double PI, ndv, L, A, R;
        vector3d<double> vel, l_vec, nxR;

        R = face_centers[n_face].norm();
        PI = a * a * u_in[0];

        l_vec[0] = u_in[1];
        l_vec[1] = u_in[2];
        l_vec[2] = u_in[3];

        vel = cross_product(face_centers[n_face], l_vec);
        vel /= (-u_in[0] * R * R);

        ndv = dot_product(edge_normals[n_face][n_edge], vel);
        nxR = cross_product(edge_normals[n_face][n_edge], face_centers[n_face]);

        if (n_edge < faces[n_face].size() - 1)
        {

            L = (vertices[faces[n_face][n_edge]] - vertices[faces[n_face][n_edge + 1]]).norm();
        }
        else
        {
            L = (vertices[faces[n_face][n_edge]] - vertices[faces[n_face][0]]).norm();
        }

        A = surface_area[n_face];

        res[0] = u_in[0] * dot_product(vel, edge_normals[n_face][n_edge]);
        res[1] = L / A * (u_in[1] * ndv + nxR[0] * PI);
        res[2] = L / A * (u_in[2] * ndv + nxR[1] * PI);
        res[3] = L / A * (u_in[3] * ndv + nxR[2] * PI);

        /*if (n_face == 1)
        {
            std::cout<<face_centers[n_face].norm()<<std::endl;
            face_centers[n_face].print();
            edge_normals[n_face][n_edge].print();
            vel.print();

            std::cout<<"=>"<<L / A<<" "<<ndv<<" "<<nxR[2]<<" "<<PI<<std::endl;


            std::cout << res[0] << " " << res[1] << " " << res[2] << " " << res[3] << std::endl;
            std::cout << std::endl;
        }*/

        /*if(n_face==1){

            std::cout<<res[0]<<" "<<res[1]<<" "<<res[2]<<" "<<res[3]<<std::endl;
        }*/

        return res;
    }
    
    std::vector<double> flux_star(std::vector<double> ul, std::vector<double> ur, int n_face, int n_edge)
    {
        std::vector<double> FL, FR, F, c_vel;
        double S_R, S_L;
        F.resize(4);

        c_vel = char_vel(ul, ur, n_face, n_edge);
        S_L = c_vel[0];
        S_R = c_vel[1];

        FL = flux(ul, n_face, n_edge);
        FR = flux(ur, n_face, n_edge);

        if (S_L >= 0)
        {
            F = FL;
        }
        else if (S_L < 0 && S_R > 0)
        {

            for (size_t i = 0; i < 4; i++)
            {
                F[i] = (S_R * FL[i] - S_L * FR[i] + S_R * S_L * (ur[i] - ul[i])) / (S_R - S_L);
            }
        }
        else if (S_R <= 0)
        {
            F = FR;
        }

        return F;
    }
    std::vector<double> char_vel(std::vector<double> u_L, std::vector<double> u_R, int n_face, int n_edge)
    {
        // returns vector {S_L, S_R}
        std::vector<double> res;
        double a_L, a_R, S_L, S_R, p_L, p_R;
        vector3d<double> vel_r, vec_r, vel_l, vec_l;

        vec_l[1] = u_L[1];
        vec_l[2] = u_L[2];
        vec_l[3] = u_L[3];

        vec_r[1] = u_R[1];
        vec_r[2] = u_R[2];
        vec_r[3] = u_R[3];

        vel_r = cross_product(face_centers[n_face], vec_r) / (-u_R[0] * (face_centers[n_face].norm()) * (face_centers[n_face].norm()));
        vel_l = cross_product(face_centers[n_face], vec_l) / (-u_L[0] * (face_centers[n_face].norm()) * (face_centers[n_face].norm()));

        p_L = a * a * u_L[0];
        p_R = a * a * u_R[0];

        a_L = sqrt(gam * p_L / u_L[0]);
        a_R = sqrt(gam * p_R / u_R[0]);

        int shift = 0;

        S_L = std::min(dot_product(vel_l, edge_normals[n_face][n_edge]), dot_product(vel_r, edge_normals[n_face][n_edge])) - std::max(a_L, a_R);
        S_R = std::max(dot_product(vel_l, edge_normals[n_face][n_edge]), dot_product(vel_r, edge_normals[n_face][n_edge])) + std::max(a_L, a_R);

        res.resize(2);
        res[0] = S_L;
        res[1] = S_R;

        return res;
    }
};
