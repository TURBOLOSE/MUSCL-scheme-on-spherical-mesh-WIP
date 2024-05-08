#pragma once
#include "MUSCL_base.hpp"
#include "MUSCL_geometry.hpp"

class MUSCL_HLLE_p : public MUSCL_base
{

private:
    std::ofstream outfile, outfile_curl, outfile_p, outfile_omega;
    double omega_ns;

public:
    MUSCL_HLLE_p(SurfaceMesh mesh, std::vector<std::vector<double>> U_in, int dim, double gam)
        : MUSCL_base(mesh, U_in, dim, gam) //U_in should be n_faces * dim(=5)
    {
         if(dim!=5){
            std::cout<<"check dim \n";
            stop_check=true;
        }
        omega_ns=5;

        outfile.open("results/rho.dat", std::ios::out | std::ios::trunc);
        outfile.close();
        outfile.open("results/rho.dat", std::ios::out | std::ios::app);

        outfile_curl.open("results/curl.dat", std::ios::out | std::ios::trunc);
        outfile_curl.close();
        outfile_curl.open("results/curl.dat", std::ios::out | std::ios::app);

        outfile_p.open("results/p.dat", std::ios::out | std::ios::trunc);
        outfile_p.close();
        outfile_p.open("results/p.dat", std::ios::out | std::ios::app);

        outfile_omega.open("results/omega.dat", std::ios::out | std::ios::trunc);
        outfile_omega.close();
        outfile_omega.open("results/omega.dat", std::ios::out | std::ios::app);
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
    void write_t_p()
    {
        vector3d<double> vel, l_vec;
        double pres;
        outfile_p << this->time() << "  ";
        for (size_t n_face = 0; n_face < faces.size(); n_face++)
        {

            l_vec[0] = U[n_face][1];
            l_vec[1] = U[n_face][2];
            l_vec[2] = U[n_face][3];

            vel = cross_product(face_centers[n_face] / face_centers[n_face].norm(), l_vec);
            vel /= (-U[n_face][0]);

            pres=pressure(U[n_face],vel, face_centers[n_face]);
            outfile_p << pres << " ";
        }
        outfile_p << "\n";
    };

    void write_t_curl()
    {
        vector3d<double> vel, l_vec, rxV, r, edge_center;
        double vort;
        outfile_curl << this->time() << "  ";
        size_t n_edge_1;
        for (size_t n_face = 0; n_face < this->n_faces(); n_face++)
        {

            // vel = cross_product(face_centers[n_face]/face_centers[n_face].norm(), l_vec);
            // vel /= (-U[n_face][0]);

            vort = 0;
            for (size_t n_edge = 0; n_edge < faces[n_face].size(); n_edge++)
            {
                l_vec[0] = U_plus[n_face][n_edge][1];
                l_vec[1] = U_plus[n_face][n_edge][2];
                l_vec[2] = U_plus[n_face][n_edge][3];

                n_edge_1 = n_edge + 1;
                if (n_edge == faces[n_face].size() - 1)
                    n_edge_1 = 0;

                edge_center = (vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]]) / 2.;
                edge_center /= edge_center.norm();

                vel = cross_product(edge_center, l_vec);
                vel /= (-U[n_face][0]);

                r = (vertices[faces[n_face][n_edge]] - vertices[faces[n_face][n_edge_1]]);
                vort += dot_product(vel, r);
            }

            // rxV = cross_product(face_centers[n_face], vel);
            // outfile_curl << rxV.norm() << " ";

            outfile_curl << vort / surface_area[n_face] << " ";
        }
        outfile_curl << "\n";
    };


    void write_t_omega_z()
    {
        vector3d<double> vel, l_vec, rxV;
        outfile_omega << this->time() << "  ";
        size_t n_edge_1;
        double theta;
        for (size_t n_face = 0; n_face < this->n_faces(); n_face++)
        {
            theta=std::acos(face_centers[n_face][2]/face_centers[n_face].norm());
            l_vec[0]=U[n_face][1]; l_vec[1]=U[n_face][2]; l_vec[2]=U[n_face][3];
            vel = cross_product(face_centers[n_face]/face_centers[n_face].norm(), l_vec);
            vel /= (-U[n_face][0]);


            rxV = cross_product(face_centers[n_face], vel);
            //outfile_omega << rxV[2] << " ";
            outfile_omega<< (vel.norm() * vel.norm()-omega_ns*omega_ns*std::sin(theta)*std::sin(theta)) / 2<<" ";
        }
        outfile_omega << "\n";
    };

protected:
    // U = {rho, l1, l2, l3, E}
    std::vector<double> flux(std::vector<double> u_in, int n_face, int n_edge)
    {

        std::vector<double> res;
        res.resize(dim);
        double PI, ndv, L, A, R;
        vector3d<double> vel, l_vec, nxR, edge_center;

        // R = face_centers[n_face].norm();

        int n_edge_1 = n_edge + 1;
        if ((n_edge_1) == faces[n_face].size())
        {
            n_edge_1 = 0;
        }

        edge_center = (vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]]) / 2.;
        edge_center /= edge_center.norm();

        l_vec[0] = u_in[1];
        l_vec[1] = u_in[2];
        l_vec[2] = u_in[3];

        vel = cross_product(edge_center, l_vec);
        vel /= (-u_in[0]) * edge_center.norm();

        //PI = (u_in[4] - u_in[0] * vel.norm() * vel.norm() / 2.) * (gam - 1);

        PI=pressure(u_in,vel,edge_center);

        ndv = dot_product(edge_normals[n_face][n_edge], vel);
        nxR = cross_product(edge_normals[n_face][n_edge], (edge_center / edge_center.norm()));

        res[0] = u_in[0] * dot_product(vel, edge_normals[n_face][n_edge]);
        res[1] = (u_in[1] * ndv - nxR[0] * PI);
        res[2] = (u_in[2] * ndv - nxR[1] * PI);
        res[3] = (u_in[3] * ndv - nxR[2] * PI);
        res[4] = (u_in[4] + PI) * dot_product(vel, edge_normals[n_face][n_edge]);

        return res;
    }

    std::vector<double> flux_star(std::vector<double> ul, std::vector<double> ur, int n_face, int n_edge)
    {

        std::vector<double> FL, FR, F, c_vel;
        double S_R, S_L;

        F.resize(dim);

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

            for (size_t i = 0; i < dim; i++)
            {
                F[i] = (S_R * FL[i] - S_L * FR[i] + S_R * S_L * (ur[i] - ul[i])) / (S_R - S_L);
            }
        }
        else if (S_R <= 0)
        {
            F = FR;
        }
        else
        {
            std::cout << "flux_star: check char vel, S_R=  " << S_R << " S_L= " << S_L << std::endl;
            stop_check = true;
        }

        return F;
    };

    std::vector<double> source(std::vector<double> u, int n_face){ // du/dt


        std::vector<double>res;
        res.resize(dim);
        vector3d<double> edge_center, l_vec, vel, RxV;

        for (size_t i = 0; i < dim; i++)
        res[i]=0;


        l_vec[0] = u[1];
        l_vec[1] = u[2];
        l_vec[2] = u[3];

        vel = cross_product(face_centers[n_face]/face_centers[n_face].norm(), l_vec);
        vel /= (-u[0]);

        //compressed star test
       double theta=std::acos(face_centers[n_face][2]/face_centers[n_face].norm());
        double phi=std::atan2(face_centers[n_face][1]/face_centers[n_face].norm(),
        face_centers[n_face][0]/face_centers[n_face].norm());


        res[1]=-omega_ns*omega_ns*u[0] * std::cos(theta) * std::sin(theta) *(-std::sin(phi));
        res[2]=-omega_ns*omega_ns*u[0] * std::cos(theta) * std::sin(theta) * std::cos(phi);




    return res;

    };

    double pressure(std::vector<double> u, vector3d<double> vel, vector3d<double> r){


        double theta=std::acos(r[2]/r.norm());
        //return  (u[4] - u[0] * vel.norm() * vel.norm() / 2) * (gam - 1); //v1 = uncompressed
        return  (u[4] - u[0] * (vel.norm() * vel.norm()-omega_ns*omega_ns*std::sin(theta)*std::sin(theta)) / 2) * (gam - 1)/gam; //v3 = compressed star + sin 

    }

    std::vector<double> char_vel(std::vector<double> u_L, std::vector<double> u_R, int n_face, int n_edge)
    {
        // returns vector {S_L, S_R}
        std::vector<double> res;
        double a_l, a_r, S_L, S_R, p_L, p_R;
        vector3d<double> vel_r, vec_r, vel_l, vec_l,edge_center_l,edge_center_r;


        int n_edge_1 = n_edge + 1;
        if ((n_edge_1) == faces[n_face].size())
        {
            n_edge_1 = 0;
        }

        edge_center_r = (vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]]) / 2.;
        edge_center_l = (vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]]) / 2.;

        edge_center_r/=edge_center_r.norm();
        edge_center_l/=edge_center_l.norm();

        vec_l[0] = u_L[1];
        vec_l[1] = u_L[2];
        vec_l[2] = u_L[3];

        vec_r[0] = u_R[1];
        vec_r[1] = u_R[2];
        vec_r[2] = u_R[3];

        vel_l=cross_product(edge_center_l, vec_l);
        vel_l /= (-u_L[0])*edge_center_l.norm();

        vel_r=cross_product(edge_center_r, vec_r);
        vel_r /= (-u_R[0])*edge_center_r.norm();


        //p_L = (u_L[4] - u_L[0] * vel_l.norm() * vel_l.norm() / 2.) * (gam - 1);
        //p_R = (u_R[4] - u_R[0] * vel_r.norm() * vel_r.norm() / 2.) * (gam - 1);

        p_L=pressure(u_L, vel_l, edge_center_l);
        p_R=pressure(u_R, vel_r, edge_center_r);

        a_l = std::sqrt(gam * p_L / u_L[0]);
        a_r = std::sqrt(gam * p_R / u_R[0]);

        int shift = 0;

        S_L = std::min(dot_product(vel_l, edge_normals[n_face][n_edge]), dot_product(vel_r, edge_normals[n_face][n_edge])) - std::max(a_l, a_r);
        S_R = std::max(dot_product(vel_l, edge_normals[n_face][n_edge]), dot_product(vel_r, edge_normals[n_face][n_edge])) + std::max(a_l, a_r);

        res.resize(2);
        res[0] = S_L;
        res[1] = S_R;
        
        /*if(abs(p_R-1)>0.1||abs(p_L-1)>0.1){
            std::cout<<"err p_r: "<<p_R<<" p_l: "<<p_L<<std::endl;
        }*/
        

         if(std::isnan(S_L)||std::isnan(S_R)){
            std::cout<<"p_r: "<<p_R<<" p_l: "<<p_L<<std::endl;
            std::cout<<"rho_l: "<<u_L[0]<<" rho_r: "<<u_R[0]<<std::endl;
            stop_check=true;
        }


        return res;
    }


    std::vector<double> limiter(std::vector<double> u_r, int n_face, int n_edge)
    { // classical Superbee limiter for irregular grids
        // CFL independent
        double etha_minus, etha_plus;
        vector3d<double> R_vec, l_vec, vel, edge_center;
        double R, c, nu_plus;
        std::vector<double> res;
        res.resize(4);

        int n_edge_1 = n_edge + 1;
        if ((n_edge_1) == faces[n_face].size())
        {
            n_edge_1 = 0;
        }
        edge_center = (vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]]) / 2.;
        edge_center /= edge_center.norm();

        l_vec[0] = U[n_face][1];
        l_vec[1] = U[n_face][2];
        l_vec[2] = U[n_face][3];

        vel = cross_product(edge_center, l_vec);
        vel /= (-U[n_face][0]) * edge_center.norm();

        //double p = (U[n_face][4] - U[n_face][0] * vel.norm() * vel.norm() / 2) * (gam - 1);
        double p = pressure(U[n_face],vel, edge_center);
        c = std::sqrt(gam * p / U[n_face][0]);

        nu_plus = (c + dot_product(vel, edge_normals[n_face][n_edge])) * dt *
                  (distance(vertices[faces[n_face][n_edge]], vertices[faces[n_face][n_edge_1]]) / surface_area[n_face]);

        etha_plus = H_plus[n_face][n_edge] / BM_dist[n_face][n_edge];
        etha_plus = H_minus[n_face][n_edge] / BM_dist[n_face][n_edge];

        for (size_t i = 0; i < dim; i++)
        {
            res[i] = std::max(0.,
                              std::max(std::min(1., etha_minus * u_r[i] * 2 / (2 * faces[n_face].size() * nu_plus)),
                                       std::min(u_r[i], etha_plus)));

            if (std::isnan(u_r[i]))
            {
                res[i] = 0;
            }
        }

        return res;
    };
};
