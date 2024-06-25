#pragma once
#include "../MUSCL_base/MUSCL_base.hpp"
#include "../geometry/MUSCL_geometry.hpp"

class adiabatic : public MUSCL_base
{

protected:
    std::ofstream outfile, outfile_curl, outfile_p, outfile_omega;
    double omega_ns;

public:
    adiabatic(SurfaceMesh mesh, std::vector<std::vector<double>> U_in, int dim, double gam, double omega_ns_i)
        : MUSCL_base(mesh, U_in, dim, gam)
    {

        omega_ns=omega_ns_i;

        set_analytical_solution();
        if (dim != 5)
        {
            std::cout << "check dim \n";
            stop_check = true;
        }
        
       

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

        /*for (size_t i = 0; i < faces.size(); i++)
        {
            outfile << surface_area[i] << " ";
        }*/
        

        outfile << "\n";
    };
    void write_t_p()
    {
        vector3d<double> vel, l_vec, edge_center;
        double pres;
        outfile_p << this->time() << "  ";
        for (size_t n_face = 0; n_face < faces.size(); n_face++)
        {

            l_vec[0] = U[n_face][1];
            l_vec[1] = U[n_face][2];
            l_vec[2] = U[n_face][3];

            vel = cross_product(face_centers[n_face] / face_centers[n_face].norm(), l_vec);
            vel /= (-U[n_face][0]);

            pres = pressure(U[n_face], vel, face_centers[n_face]);
            outfile_p << pres << " ";
            //outfile_p <<U[n_face][4] << " ";
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
            theta = std::acos(face_centers[n_face][2] / face_centers[n_face].norm());
            l_vec[0] = U[n_face][1];
            l_vec[1] = U[n_face][2];
            l_vec[2] = U[n_face][3];
            vel = cross_product(face_centers[n_face] / face_centers[n_face].norm(), l_vec);
            vel /= (-U[n_face][0]);

            rxV = cross_product(face_centers[n_face], vel);
            outfile_omega << rxV[2] << " ";

            // std::cout<<std::setprecision(9)<<face_centers[n_face][2]<<"\n";
            //outfile_omega << (vel.norm() * vel.norm() - omega_ns * omega_ns * std::sin(theta) * std::sin(theta)) / 2 << " ";

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

        // PI = (u_in[4] - u_in[0] * vel.norm() * vel.norm() / 2.) * (gam - 1);
        PI = pressure(u_in, vel, edge_center);

        ndv = dot_product(edge_normals[n_face][n_edge], vel);
        nxR = cross_product(edge_normals[n_face][n_edge], (edge_center / edge_center.norm()));

        double theta = std::acos(edge_center[2] / edge_center.norm());

        res[0] = u_in[0] * dot_product(vel, edge_normals[n_face][n_edge]);
        res[1] = (u_in[1] * ndv - nxR[0] * PI);
        res[2] = (u_in[2] * ndv - nxR[1] * PI);
        res[3] = (u_in[3] * ndv - nxR[2] * PI);
        res[4] = (u_in[4] + PI) * dot_product(vel, edge_normals[n_face][n_edge]);


        return res;
    }

    virtual std::vector<double> flux_star(std::vector<double> ul, std::vector<double> ur, int n_face, int n_edge) = 0;
    

    std::vector<double> source(std::vector<double> u, int n_face)
    { // du/dt

        std::vector<double> res;
        res.resize(dim);
        vector3d<double> l_vec, vel, vel_dot, omega_acc, omxr, rxv, fc_normed, edge_center;

        double tilt_angle = M_PI / 30;

        for (size_t i = 0; i < dim; i++)
            res[i] = 0;

        l_vec[0] = u[1];
        l_vec[1] = u[2];
        l_vec[2] = u[3];

        fc_normed = face_centers[n_face] / face_centers[n_face].norm();
        vel = cross_product(fc_normed, l_vec);
        vel /= (-u[0]);

        // compressed star test
        double theta = std::acos(face_centers[n_face][2] / face_centers[n_face].norm());

        double phi = std::atan2(face_centers[n_face][1] / face_centers[n_face].norm(), face_centers[n_face][0] / face_centers[n_face].norm());
        //res[1] = -omega_ns * omega_ns * u[0] * std::cos(theta) * std::sin(theta) * (-std::sin(phi)); // x
        //res[2] = -omega_ns * omega_ns * u[0] * std::cos(theta) * std::sin(theta) * std::cos(phi);    // y

        int n_edge_1;
        res[1]=0;
        res[2]=0;
        for (size_t n_edge = 0; n_edge < faces[n_face].size(); n_edge++)
        {   
            n_edge_1=n_edge+1;
            if(n_edge==faces[n_face].size()-1)
            n_edge_1=0;
            edge_center = (vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]]);
            edge_center/=edge_center.norm();
            theta = std::acos(edge_center[2]);
            res[1] += -omega_ns * omega_ns * (U_plus[n_face][n_edge][0]+U_minus[n_face][n_edge][0])/2. 
            * std::cos(theta) * std::sin(theta) * (-std::sin(phi)); // x
            res[2] += -omega_ns * omega_ns * (U_plus[n_face][n_edge][0]+U_minus[n_face][n_edge][0])/2. 
            * std::cos(theta) * std::sin(theta) * std::cos(phi);    // y

        };
        
        res[1]/=faces[n_face].size();
        res[2]/=faces[n_face].size();

        //=============================================
        //p_an[n_face]=3.2e-06+(vel.norm()*vel.norm()-omega_ns*omega_ns*std::sin(theta)*std::sin(theta))/2.;
        //=============================================


         /*if(std::abs(face_centers[n_face][2]*std::cos(tilt_angle) +face_centers[n_face][1]*std::sin(tilt_angle))  <0.1){
         //res[0]=1.6e-6; //= 10^-8 M_sun/yr
         res[0]=1.6e-2; //= 10^-4 M_sun/yr
         //res[0]=0.16; //= 10^-3 M_sun/yr
         omega_acc[0]=0; omega_acc[1]=std::sin(tilt_angle)*0.1; omega_acc[2]=std::cos(tilt_angle)*0.1;

         omxr=cross_product(omega_acc, fc_normed);
         rxv=cross_product(fc_normed, omxr);

         res[1]+=res[0]*rxv[0];
         res[2]+=res[0]*rxv[1];
         res[3]+=res[0]*rxv[2];


         //res[4]=(omxr.norm()*omxr.norm())/2. * res[0];
         res[4]=(omxr.norm()*omxr.norm()-omega_ns*omega_ns*std::sin(theta)*std::sin(theta))/2. * res[0];
         }*/
         

        // res[4] = -omega_ns * omega_ns * u[0] * std::sin(theta) * std::cos(theta) * (-std::sin(phi)*l_vec[0]+std::cos(phi)*l_vec[1]); //v2

        return res;
    };

    double pressure(std::vector<double> u, vector3d<double> vel, vector3d<double> r)
    {
        double theta = std::acos(r[2] / r.norm());

        // return  (u[4] - u[0] * vel.norm() * vel.norm() / 2) * (gam - 1); //v1 = uncompressed
        // return (u[4] - u[0] * vel.norm() * vel.norm() / 2) * (gam - 1) / gam; // v2 = different P
        //return (u[4] - u[0] * (vel.norm() * vel.norm() - omega_ns * omega_ns * std::sin(theta) * std::sin(theta)) / 2) * (gam - 1) / gam; // v3 = compressed star + sin
        return (u[4] - u[0] * (vel.norm() * vel.norm() - omega_ns * omega_ns * std::sin(theta) * std::sin(theta)) / 2) * (gam - 1); // v4 = compressed star new gamma

    }

    std::vector<double> char_vel(std::vector<double> u_L, std::vector<double> u_R, int n_face, int n_edge)
    {
        // returns vector {S_L, S_R}
        std::vector<double> res;
        double a_L, a_R, S_L, S_R, p_L, p_R, q_R, q_L;
        vector3d<double> vel_r, vec_r, vel_l, vec_l, edge_center_l, edge_center_r;

        int n_edge_1 = n_edge + 1;
        if ((n_edge_1) == faces[n_face].size())
        {
            n_edge_1 = 0;
        }

        edge_center_r = (vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]]) / 2.;
        edge_center_l = (vertices[faces[n_face][n_edge]] + vertices[faces[n_face][n_edge_1]]) / 2.;

        edge_center_r /= edge_center_r.norm();
        edge_center_l /= edge_center_l.norm();

        vec_l[0] = u_L[1];
        vec_l[1] = u_L[2];
        vec_l[2] = u_L[3];

        vec_r[0] = u_R[1];
        vec_r[1] = u_R[2];
        vec_r[2] = u_R[3];

        vel_l = cross_product(edge_center_l, vec_l);
        vel_l /= (-u_L[0]) * edge_center_l.norm();

        vel_r = cross_product(edge_center_r, vec_r);
        vel_r /= (-u_R[0]) * edge_center_r.norm();

        // p_L = (u_L[4] - u_L[0] * vel_l.norm() * vel_l.norm() / 2) * (gam - 1);
        // p_R = (u_R[4] - u_R[0] * vel_r.norm() * vel_r.norm() / 2) * (gam - 1);

        p_L = pressure(u_L, vel_l, edge_center_l);
        p_R = pressure(u_R, vel_r, edge_center_r);

        a_L = std::sqrt(gam * p_L / u_L[0]);
        a_R = std::sqrt(gam * p_R / u_R[0]);

        double z = (gam - 1) / (2 * gam);
        // double p_star=std::pow((a_L+a_R-(gam-1)/2. * (dot_product(edge_normals[n_face][n_edge], vel_r)-dot_product(edge_normals[n_face][n_edge], vel_l)))
        //  /( a_L/std::pow(p_L,z) + a_R/std::pow(p_R,z) ),1./z);

        double p_pvrs = 1 / 2. * (p_L + p_R) - 1 / 8. * (dot_product(edge_normals[n_face][n_edge], vel_r) - dot_product(edge_normals[n_face][n_edge], vel_l)) * (u_L[0] + u_R[0]) * (a_L + a_R);
        double p_star = std::max(0., p_pvrs);

        if (p_star <= p_R)
        {
            q_R = 1;
        }
        else
        {
            q_R = std::sqrt(1 + (gam + 1) / (2 * gam) * (p_star / p_R - 1));
        }

        if (p_star <= p_L)
        {
            q_L = 1;
        }
        else
        {
            q_L = std::sqrt(1 + (gam + 1) / (2 * gam) * (p_star / p_L - 1));
        }

        S_L = dot_product(vel_l, edge_normals[n_face][n_edge]) - a_L * q_L;
        S_R = dot_product(vel_r, edge_normals[n_face][n_edge]) + a_R * q_R;

        // S_L = std::min(dot_product(vel_l, edge_normals[n_face][n_edge]), dot_product(vel_r, edge_normals[n_face][n_edge])) - std::max(a_L, a_R);
        // S_R = std::max(dot_product(vel_l, edge_normals[n_face][n_edge]), dot_product(vel_r, edge_normals[n_face][n_edge])) + std::max(a_L, a_R);

        if (std::isnan(S_L) || std::isnan(S_R))
        {
            std::cout<<"S_L or S_R is NaN"<<"\n";
            stop_check = true;
        }

        res.resize(2);
        res[0] = S_L;
        res[1] = S_R;

        return res;
    }

    std::vector<double> limiter(std::vector<double> u_r, int n_face, int n_edge)
    {
        //std::cout<<u_r[0]<<" "<<u_r[1]<<" "<<u_r[2]<<" "<<u_r[3]<<" "<<u_r[4]<<"\n";

        std::vector<double> res;
        res.resize(dim);

        double a = 4, b = 2, c = 0.1, d = 10, e = 3, f = 6; // switch function parameters
        auto h = [a, b, c, d, e, f](double r)
        {
            double res = 0;
            if (r < 1 && r > 0)
                res = (1 - std::tanh(a * std::pow(r, b) * std::pow(1 - r, c)));
            if (r >= 1)
                res = std::pow(std::tanh(d * std::pow(r - 1, e)), f);

            return res;
        };

        std::vector<double> to = limiter_third_order(u_r, n_face, n_edge);
        std::vector<double> sb = limiter_superbee(u_r, n_face, n_edge);

        for (size_t i = 0; i < dim; i++)
        {
            res[i] = ((1 - h(u_r[i])) * to[i] + h(u_r[i]) * sb[i]);
            //res[i] = ((1 - h(u_r[0])) * to[0] + h(u_r[0]) * sb[0]);
            // res[i] = 0;
            if (std::isnan(u_r[i]))
            {
                res[i] = 0;
            }
            //res[i]=1;
        }


        //std::cout<<to[4]<<" "<<sb[4]<<" "<<res[4]<<"\n";
         //return limiter_third_order(u_r, n_face, n_edge);
        //return limiter_superbee(u_r, n_face, n_edge);
        return res;
    }

    std::vector<double> limiter_third_order(std::vector<double> u_r, int n_face, int n_edge)
    {
        std::vector<double> supb = limiter_superbee(u_r, n_face, n_edge);
        vector3d<double> R_vec, l_vec, vel, edge_center;
        double R, c, nu_plus;
        std::vector<double> res;
        res.resize(dim);

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

        //tag1
        //vel /= (-U[n_face][0]) * edge_center.norm();

        vel /= (-U[n_face][0]-rho_an[n_face]) * edge_center.norm();
        //double p = pressure(U[n_face], vel, edge_center);
        double p=U[n_face][4]+p_an[n_face];

        c = std::sqrt(gam * p / (U[n_face][0]+rho_an[n_face]));
        //c = std::sqrt(gam * p / U[n_face][0]);

        nu_plus = (c + dot_product(vel, edge_normals[n_face][n_edge])) * dt *
                  (distance(vertices[faces[n_face][n_edge]], vertices[faces[n_face][n_edge_1]]) / surface_area[n_face]);


        for (size_t i = 0; i < dim; i++)
        {

            res[i] = std::max(0., std::min(supb[i], 1 + (1 + nu_plus) / 3 * (u_r[i] - 1)));
            

            if (std::isnan(u_r[i]))
            {
                res[i] = 0;
            }
        }

        return res;
    }

    std::vector<double> limiter_superbee(std::vector<double> u_r, int n_face, int n_edge)
    { // classical Superbee limiter for irregular grids
        // CFL independent
        double etha_minus, etha_plus;
        vector3d<double> R_vec, l_vec, vel, edge_center;
        double R, c, nu_plus;
        std::vector<double> res;
        res.resize(dim);

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
        //tag1
        //vel /= (-U[n_face][0]) * edge_center.norm();

        vel /= (-U[n_face][0]-rho_an[n_face]) * edge_center.norm();
        //double p = pressure(U[n_face], vel, edge_center);

        double p=U[n_face][4]+p_an[n_face];

        //c = std::sqrt(gam * p / U[n_face][0]);

        c = std::sqrt(gam * p / (U[n_face][0]+rho_an[n_face]));


        nu_plus = (c + dot_product(vel, edge_normals[n_face][n_edge])) * dt *
                  (distance(vertices[faces[n_face][n_edge]], vertices[faces[n_face][n_edge_1]]) / surface_area[n_face]);

        etha_plus = H_plus[n_face][n_edge] / BM_dist[n_face][n_edge];
        etha_minus = H_minus[n_face][n_edge] / BM_dist[n_face][n_edge];



        for (size_t i = 0; i < dim; i++)
        {

            res[i] = std::max(0.,
                              std::max(std::min(1., etha_minus * u_r[i] / (2 * faces[n_face].size() * nu_plus)),
                                       std::min(u_r[i], etha_plus)));
            if (std::isnan(u_r[i]))
            {
                res[i] = 0;
            }
        }

        return res;
    };

    void set_analytical_solution()// analytical solution to be preserved                               
    {                             // if no AS is required, thish should set rho_an and p_an to 0
        vector3d<double> vec_l, vel,r;
        for (size_t i = 0; i < faces.size(); i++)
        {
            vec_l[0] = U[i][1];
            vec_l[1] = U[i][2];
            vec_l[2] = U[i][3];


            vel = cross_product(face_centers[i]/face_centers[i].norm(), vec_l);
            vel /= -U[i][0];

            //p_an[i] = pressure(U[i], vel, face_centers[i]);
            //rho_an[i] = U[i][0];   //will try to conserve current profile

            rho_an[i] = 0;   //no profile to be conserved
            p_an[i] = 0;
        }
    }

};
