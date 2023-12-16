#pragma once
#include "MUSCL_base.hpp"
#include "MUSCL_geometry.hpp"

class MUSCL_HLLE : public MUSCL_base
{

private:
    std::ofstream outfile, outfile_curl;

public:
    MUSCL_HLLE(SurfaceMesh mesh, std::vector<std::vector<double>> U_in, int dim, double gam)
        : MUSCL_base(mesh, U_in, dim, gam)
    {
        outfile.open("results/rho.dat", std::ios::out | std::ios::trunc);
        outfile.close();
        outfile.open("results/rho.dat", std::ios::out | std::ios::app);

        outfile_curl.open("results/curl.dat", std::ios::out | std::ios::trunc);
        outfile_curl.close();
        outfile_curl.open("results/curl.dat", std::ios::out | std::ios::app);
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

    void write_t_curl()
    {
         vector3d<double> vel, l_vec, rxV;
        outfile_curl << this->time() << "  ";
        for (size_t n_face = 0; n_face < this->n_faces(); n_face++)
        {
            
            l_vec[0] = U[n_face][1];
            l_vec[1] = U[n_face][2];
            l_vec[2] = U[n_face][3];


            vel = cross_product(face_centers[n_face], l_vec);
            vel /= (-U[n_face][0] * (face_centers[n_face].norm() * face_centers[n_face].norm()));

            rxV=cross_product(face_centers[n_face],vel);
            outfile_curl << rxV.norm()<< " ";
        }
        outfile_curl << "\n";
    };

public:
    const double a = 1;

    // U = {rho, l1, l2, l3}
    std::vector<double> flux(std::vector<double> u_in, int n_face, int n_edge)
    {
        std::vector<double> res;
        res.resize(4);
        double PI, ndv, L, A, R;
        vector3d<double> R_vec,vel, l_vec, nxR,edge_center;

        int n_edge_1=n_edge+1;
        if((n_edge_1)==faces[n_face].size()){
            n_edge_1=0;
        }


        edge_center=(vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]])/2.;

        PI = a * a * u_in[0];


        R_vec=face_centers[n_face];
        //R_vec=edge_center;
        R=R_vec.norm();


        l_vec[0] = u_in[1];
        l_vec[1] = u_in[2];
        l_vec[2] = u_in[3];

        vel = cross_product(R_vec, l_vec);
        vel /= (-u_in[0] * R * R);


        ndv = dot_product(edge_normals[n_face][n_edge], vel);
        nxR = cross_product(edge_normals[n_face][n_edge], R_vec);

        if (n_edge < faces[n_face].size() - 1)
        {

            L = (vertices[faces[n_face][n_edge]] - vertices[faces[n_face][n_edge + 1]]).norm();
        }
        else
        {
            L = (vertices[faces[n_face][n_edge]] - vertices[faces[n_face][0]]).norm();
        }

        A = surface_area[n_face];


        //std::cout<<dot_product(vel, edge_normals[n_face][n_edge])<<std::endl;


        res[0] = u_in[0] * dot_product(vel, edge_normals[n_face][n_edge]);
        res[1] = L / A * (u_in[1] * ndv + nxR[0] * PI);
        res[2] = L / A * (u_in[2] * ndv + nxR[1] * PI);
        res[3] = L / A * (u_in[3] * ndv + nxR[2] * PI);

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
        vector3d<double> vel_r, vec_r, vel_l, vec_l,edge_center;

        int n_edge_1=n_edge+1;
        if((n_edge_1)==faces[n_face].size()){
            n_edge_1=0;
        }


        edge_center=face_centers[n_face];
        //edge_center=(vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]])/2.;



        vec_l[1] = u_L[1];
        vec_l[2] = u_L[2];
        vec_l[3] = u_L[3];

        vec_r[1] = u_R[1];
        vec_r[2] = u_R[2];
        vec_r[3] = u_R[3];

        vel_r = cross_product(edge_center, vec_r) / (-u_R[0] * (edge_center.norm()) * (edge_center.norm()));
        vel_l = cross_product(edge_center, vec_l) / (-u_L[0] * (edge_center.norm()) * (edge_center.norm()));

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
