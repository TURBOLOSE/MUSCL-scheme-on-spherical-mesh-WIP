#pragma once
#include "../MUSCL_base/MUSCL_base.hpp"
#include "../geometry/MUSCL_geometry.hpp"
#include "../json.hpp"

using json = nlohmann::json;


class adiabatic : public MUSCL_base
{

protected:
    std::ofstream outfile, outfile_p, outfile_omega,outfile_curl;
    std::ofstream outfile_l[3];
    bool accretion_on;
    double total_mass, acc_rate, e_acc, omega_acc_abs;

public:
    adiabatic(SurfaceMesh mesh, std::vector<std::vector<double>> U_in, 
    int dim, double gam, double omega_ns_i, bool accretion_on_i, size_t threads)
        : MUSCL_base(mesh, U_in, dim, gam, omega_ns_i, threads), accretion_on(accretion_on_i)
    {

        omega_ns=omega_ns_i;

        set_analytical_solution();
        if (dim != 5)
        {
            std::cout << "check dim \n";
            stop_check = true;
        }
        
        //clean file and append next line
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

        std::string adrs[] = {"results/Lx.dat", "results/Ly.dat","results/Lz.dat"};

        for(size_t i; i<3; i++){
        
        outfile_l[i].open(adrs[i], std::ios::out | std::ios::trunc);
        outfile_l[i].close();
        outfile_l[i].open(adrs[i], std::ios::out | std::ios::app);

        }

       
        std::ifstream ifs("input/parameters.json");
        json parameters = json::parse(ifs);

        acc_rate=parameters["accretion_rate"];
        e_acc=parameters["e_acc"];
        omega_acc_abs=parameters["omega_acc_abs"];

        total_mass=0;
        for (size_t n_face = 0; n_face < faces.size(); n_face++)
        {
            total_mass+=U[n_face][0]*surface_area[n_face];
        }
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
            //out_lc<< U_i[0] << " ";
        }
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

    void write_t_L()
    {

        for(size_t i; i<3; i++){

            outfile_l[i] << this->time() << "  ";
            for (auto U_j : U)
            {

                outfile_l[i] << U_j[i+1] << " ";
                //out_lc<< U_i[0] << " ";
            }
            outfile_l[i] << "\n";
        }
    };

    void write_t_omega_z()
    {
        vector3d<double> vel, l_vec, rxV;
        outfile_omega << this->time() << "  ";
        size_t n_edge_1;
        double theta,om_z;
        for (size_t n_face = 0; n_face < this->n_faces(); n_face++)
        {
            theta = std::acos(face_centers[n_face][2] / face_centers[n_face].norm());
            l_vec[0] = U[n_face][1];
            l_vec[1] = U[n_face][2];
            l_vec[2] = U[n_face][3];
            vel = cross_product(face_centers[n_face] / face_centers[n_face].norm(), l_vec);
            vel /= (-U[n_face][0]);

            rxV = cross_product(face_centers[n_face]/ face_centers[n_face].norm(), vel);
            om_z=rxV[2]/(std::sin(theta)*std::sin(theta));

            if(std::isnan(om_z)||std::isinf(om_z)||std::abs(om_z)>1e4){
                //std::cout<<rxV[2]<<" "<<std::sin(theta)<<"\n";
                //(face_centers[n_face]/ face_centers[n_face].norm()).print();
                om_z=0;
            }

            outfile_omega << om_z  << " ";

            // std::cout<<std::setprecision(9)<<face_centers[n_face][2]<<"\n";
            //outfile_omega << (vel.norm() * vel.norm() - omega_ns * omega_ns * std::sin(theta) * std::sin(theta)) / 2 << " ";

        }
        outfile_omega << "\n";
    };

    std::vector<double> get_light_curves(){
        //roation around z axis is implied
        std::vector<double> result;

        vector3d<double> obs_vector_0, obs_vector_45, obs_vector_90, l_vec, vel;
        obs_vector_0[0]=9.4*3*1e19/15000; obs_vector_0[1]=0; obs_vector_0[2]=0; // dist = 9400pc (in R_ns)
        obs_vector_45[0]=std::sqrt(9.4*3*1e19/15000); obs_vector_45[1]=0; obs_vector_45[2]=std::sqrt(9.4*3*1e19/15000); // dist = 9400pc (in R_ns)
        obs_vector_90[0]=0; obs_vector_90[1]=0; obs_vector_90[2]=9.4*3*1e19/15000; // dist = 9400pc (in R_ns)

        double flux_tot_0=0,flux_tot_45=0,flux_tot_90=0; 
        double phi_fc, theta_fc, d_vec, cos_alpha,PI;

        double GM=0.217909; //grav parameter in R_unit^3/t_unit^2
        double g_eff,C, beta,E;
        double c_sigma=4.85e36; //c/sigma_SB in R_unit*t_unit^2*K^4/M_unit
        double k_m=1.6e-13; // k/m in V_unit(speed of light)^2/K
        double kappa=3.4e6; //scattering opacity in 1/Sigma_unit (R_unit^2/M_unit)
        double beta_switch=0.736194670678821; //switch point for function
        double C_switch = -1-1/(beta_switch-1);
        double beta_ceil=1-1e-8, beta_floor=1e-8;

        for (size_t n_face = 0; n_face < this->n_faces(); n_face++)
        {
            phi_fc=std::atan2(face_centers[n_face][1]/face_centers[n_face].norm(), 
            face_centers[n_face][0]/face_centers[n_face].norm());
            theta_fc = std::acos(face_centers[n_face][2] / face_centers[n_face].norm());

            l_vec[0] = U[n_face][1];
            l_vec[1] = U[n_face][2];
            l_vec[2] = U[n_face][3];

            vel = cross_product(face_centers[n_face]/face_centers[n_face].norm(), l_vec);
            vel /= (-U[n_face][0]);
            PI = pressure(U[n_face], vel, face_centers[n_face]/face_centers[n_face].norm());

            g_eff=GM-vel.norm()*vel.norm();

            C=12./5*k_m*U[n_face][0]/PI*pow(3./4 *c_sigma*g_eff*U[n_face][0],1./4);

            if(C<=C_switch){
                
                beta=1-1/(1+C);
            }else{
                beta=1-pow(2/C,4);
            }

            if(beta<beta_floor|| std::isnan(beta)) //beta limitations
            beta=beta_floor;

            if(beta>beta_ceil)
            beta=beta_ceil;

        
            E=g_eff/kappa*(1-beta);

            if(phi_fc<M_PI/2 && phi_fc >-M_PI/2){
            d_vec=dot_product(obs_vector_0, face_centers[n_face]/face_centers[n_face].norm());
            cos_alpha=std::abs(d_vec)/obs_vector_0.norm();
            //flux_tot_0+=PI*cos_alpha*surface_area[n_face];
            flux_tot_0+=E*cos_alpha*surface_area[n_face];
            }

            if(theta_fc<M_PI/2){
            d_vec=dot_product(obs_vector_90, face_centers[n_face]/face_centers[n_face].norm());
            cos_alpha=std::abs(d_vec)/obs_vector_90.norm();
            //flux_tot_90+=PI*cos_alpha*surface_area[n_face];
            flux_tot_90+=E*cos_alpha*surface_area[n_face];
            }

            if(dot_product(obs_vector_45, face_centers[n_face]/face_centers[n_face].norm())>0){
            d_vec=dot_product(obs_vector_45, face_centers[n_face]/face_centers[n_face].norm());
            cos_alpha=std::abs(d_vec)/obs_vector_45.norm();
            //flux_tot_45+=PI*cos_alpha*surface_area[n_face];
            flux_tot_45+=E*cos_alpha*surface_area[n_face];
            }
        }

        result.push_back(flux_tot_0); result.push_back(flux_tot_45); result.push_back(flux_tot_90);
        return result;
    };

    void write_final_state(){
        std::ofstream out_final;
        out_final.open("results/final_state.dat");

        for (size_t n_face = 0; n_face < this->n_faces(); n_face++)
        {
            for(size_t j = 0; j < dim; j++)
            {
                out_final<<U[n_face][j]<<" ";
            }
            out_final<<"\n";
        }

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
      //due to the weird reconstruction algorithm, u[4] here is the pressure \Pi


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
        //old stuff for compressed star

    
        //res[1] = -omega_ns * omega_ns * u[0] * std::cos(theta) * std::sin(theta) * (-std::sin(phi)); // x
        //res[2] = -omega_ns * omega_ns * u[0] * std::cos(theta) * std::sin(theta) * std::cos(phi);    // y

        


        //accretion terms
        //double acc_rate=3.3e-5;//= 10^-7 M_sun/yr (d \Sigma/dt in local units)
        //double acc_rate; //= 10^-8 M_sun/yr
        //double acc_rate=3.3e-2; //= 10^-4 M_sun/yr
        //double acc_rate=0.33; //= 10^-3 M_sun/yr

        //double e_acc=1e-6; //local E units
        //double e_acc=6e-3; //local E units

        if(accretion_on && std::abs(face_centers[n_face][2]*std::cos(tilt_angle) +face_centers[n_face][1]*std::sin(tilt_angle))  <0.1)
        {
            res[0]=acc_rate;
            omega_acc[0]=0; omega_acc[1]=std::sin(tilt_angle)*omega_acc_abs; omega_acc[2]=std::cos(tilt_angle)*omega_acc_abs;

            omxr=cross_product(omega_acc, fc_normed);
            rxv=cross_product(fc_normed, omxr);

            res[1]+=acc_rate*rxv[0];
            res[2]+=acc_rate*rxv[1];
            res[3]+=acc_rate*rxv[2];
            //res[4]=(omxr.norm()*omxr.norm())/2. * res[0];
            res[4]= acc_rate *( e_acc+ ((omxr-vel).norm()*(omxr-vel).norm()-omega_ns*omega_ns*std::sin(theta)*std::sin(theta))/2. );
        }



       
        if(accretion_on)
        {
             //sink term for mass into crust (+ energy and angular momentum)
            rxv=cross_product(fc_normed, vel);
            //double dmdt=-1./512100000*u[0]; //time constant in T_unit^{-1}
            double dmdt=-(acc_rate*0.4*M_PI)/total_mass * u[0]; //time constant in T_unit^{-1} times surf density
            //double dmdt=0;
            res[0]+=dmdt;
            res[1]+=dmdt*rxv[0];
            res[2]+=dmdt*rxv[1];
            res[3]+=dmdt*rxv[2];
            res[4]+=dmdt*(1/(gam-1)*u[4]/u[0]+ (vel.norm()*vel.norm()-omega_ns*omega_ns*std::sin(theta)*std::sin(theta))/2.) ;


            //energy sink term (radiation energy diffusion)
            double GM=0.217909; //grav parameter in R_unit^3/t_unit^2
            double g_eff=GM-vel.norm()*vel.norm();
            //double g_eff=GM;
            double c_sigma=4.85e36; //c/sigma_SB in R_unit*t_unit^2*K^4/M_unit
            double k_m=1.6e-13; // k/m in V_unit(speed of light)^2/K
            double kappa=3.4e6; //scattering opacity in 1/Sigma_unit (R_unit^2/M_unit)

            double C=12./5*k_m*u[0]/u[4]*pow(3./4 *c_sigma*g_eff*u[0],1./4);

            double beta_switch=0.736194670678821; //switch point for function
            double C_switch = -1-1/(beta_switch-1);
            double beta;

            if(C<=C_switch){
                
                beta=1-1/(1+C);
            }else{
                beta=1-pow(2/C,4);
            }


            double beta_ceil=1-1e-8, beta_floor=1e-8;

            if(beta<beta_floor|| std::isnan(beta)) //beta limitations
            beta=beta_floor;

            if(beta>beta_ceil)
            beta=beta_ceil;

          
            
            res[4]-=g_eff/kappa*(1-beta);




            //if(n_face==2720){
            // std::cout<<acc_rate *( e_acc+ ((omxr-vel).norm()*(omxr-vel).norm()-omega_ns*omega_ns*std::sin(theta)*std::sin(theta))/2. )<<" "
            // <<dmdt*(gam/(gam-1)*u[4]/u[0]+ (vel.norm()*vel.norm()-omega_ns*omega_ns*std::sin(theta)*std::sin(theta))/2.) <<" "
            // <<g_eff/kappa*(1-beta)<<" beta= "<<beta<< "\n";

            //std::cout<<beta<< "\n";
            //}

            // static double max_E=0;
            // if(std::abs(res[4]/(1/(gam-1)*u[4]+vel.norm()*vel.norm()/2))>max_E){
            // max_E=std::abs(res[4]/(1/(gam-1)*u[4]+vel.norm()*vel.norm()/2));
            // std::cout<<max_E<<"\n";
            // }

            /*static double max_Mach=0;
            if(vel.norm()/std::sqrt(gam*u[4]/u[0])>max_Mach){
            max_Mach=vel.norm()/std::sqrt(gam*u[4]/u[0]);
            std::cout<<max_Mach<<"\n";
            }*/

            //if(this->time()>590)
            //std::cout<<max_Mach<<"\n";

            
        }

            

        return res;
    };

    double pressure(std::vector<double> u, vector3d<double> vel, vector3d<double> r)
    {
        double theta = std::acos(r[2] / r.norm());

        double pressure_floor=1e-16;
        // return  (u[4] - u[0] * vel.norm() * vel.norm() / 2) * (gam - 1); //v1 = uncompressed
        // return (u[4] - u[0] * vel.norm() * vel.norm() / 2) * (gam - 1) / gam; // v2 = different P
        //return (u[4] - u[0] * (vel.norm() * vel.norm() - omega_ns * omega_ns * std::sin(theta) * std::sin(theta)) / 2) * (gam - 1) / gam; // v3 = compressed star + sin
        return std::max(pressure_floor,(u[4] - u[0] * (vel.norm() * vel.norm() - omega_ns * omega_ns * std::sin(theta) * std::sin(theta)) / 2) * (gam - 1)); // v4 = compressed star new gamma

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

        p_L =  pressure(u_L, vel_l, edge_center_l);
        p_R =  pressure(u_R, vel_r, edge_center_r);




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
            if(std::isnan(q_R)||std::isinf(q_R))
            q_R = 1;
        }

        if (p_star <= p_L)
        {
            q_L = 1;
        }
        else
        {
            q_L = std::sqrt(1 + (gam + 1) / (2 * gam) * (p_star / p_L - 1));
            if(std::isnan(q_L)||std::isinf(q_L))
            q_L = 1;
        }

        S_L = dot_product(vel_l, edge_normals[n_face][n_edge]) - a_L * q_L;
        S_R = dot_product(vel_r, edge_normals[n_face][n_edge]) + a_R * q_R;

        // S_L = std::min(dot_product(vel_l, edge_normals[n_face][n_edge]), dot_product(vel_r, edge_normals[n_face][n_edge])) - std::max(a_L, a_R);
        // S_R = std::max(dot_product(vel_l, edge_normals[n_face][n_edge]), dot_product(vel_r, edge_normals[n_face][n_edge])) + std::max(a_L, a_R);

        if (std::isnan(S_L) || std::isnan(S_R))
        {
            std::cout<<"S_L or S_R is NaN in face "<<n_face<<"\n";
            std::cout<<"p_L/R = "<<p_L<<" "<<p_R<<" rho_L/R = "<<u_L[0]<<" "<<u_R[0]<<"\n";
            stop_check = true;
        }

        res.resize(2);
        res[0] = S_L;
        res[1] = S_R;

        return res;
    }

    std::vector<double> limiter(std::vector<double> u_r, int n_face, int n_edge)
    {   // here U[4] is also pressure
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
    { // here U[4] is also pressure
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
    {   // here U[4] is also pressure
        // classical Superbee limiter for irregular grids
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
