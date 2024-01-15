#pragma once

#include "MUSCL_geometry.hpp"

class MUSCL_base : public MUSCL_base_geometry
{
protected:
    std::vector<std::vector<double>> U, U_temp;
    std::vector<std::vector<std::vector<double>>> flux_var_plus, flux_var_minus, U_plus, U_minus;
    // flux_var^plus_ij flux_var^minus_ij, U_ij (short), U_ji(short)
    double dt, gam, M, N, h0, t, max_vel, rho_full;
    int dim;
    size_t steps;

public:
    MUSCL_base(SurfaceMesh mesh, std::vector<std::vector<double>> U_in, int dim, double gam) : MUSCL_base_geometry(mesh), U(U_in), gam(gam), dim(dim)
    {

        double h0_temp;

        if (U_in.size() != mesh.n_faces())
        {
            std::cout << "U_in is wrong size, U_in.size() should be equal to mesh.n_faces()" << std::endl;
        }

        // memory allocation
        flux_var_plus.resize(this->n_faces());
        flux_var_minus.resize(this->n_faces());
        U_plus.resize(this->n_faces());
        U_minus.resize(this->n_faces());
        U_temp.resize(this->n_faces());

        for (size_t i = 0; i < this->n_faces(); i++)
        {

            flux_var_plus[i].resize(faces[i].size());
            flux_var_minus[i].resize(faces[i].size());
            U_plus[i].resize(faces[i].size());
            U_minus[i].resize(faces[i].size());
            U_temp[i].resize(dim);

            for (size_t j = 0; j < faces[i].size(); j++)
            {
                flux_var_plus[i][j].resize(dim);
                flux_var_minus[i][j].resize(dim);
                U_plus[i][j].resize(dim);
                U_minus[i][j].resize(dim);
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

        t = 0;
        steps = 0;
    };

    void do_step(double dt0)
    { // RK2

        dt = dt0;
        U_temp = U;
        find_U_edges();
        find_flux_var();
        // find_M();
        find_v_max();

        /*if(dt > h0 / (2 * M * N)){
            dt=h0 / (2 * M * N);
        }*/

        if (dt > h0 * 0.1 / max_vel)
        {
            dt = h0 * 0.1 / max_vel;
        }

        // std::cout<<h0 / (2 * M * N)<<std::endl;

        res2d(dt / 2.); // res2d makes U = dt/2*phi(U)


        //std::cout<<std::endl;
        for (size_t i = 0; i < this->n_faces(); i++)
        {

            for (size_t k = 0; k < dim; k++)
            {
                U[i][k] += U_temp[i][k]; // results in U=U+dt/2*phi(U)
            }

            //std::cout<<U[i][0]<<" "<<U[i][1]<<" "<<U[i][2]<<" "<<U[i][3]<<std::endl;
        }
        //std::cout<<std::endl;

        find_U_edges();
        find_flux_var();
        res2d(dt); // U=dt*phi(U+dt/2*phi(U))

        double temp = 0;
        double l1, l2, l3 = 0;

        l1=0; l2=0;
        for (size_t i = 0; i < this->n_faces(); i++)
        {

            for (size_t k = 0; k < dim; k++)
            {
                U[i][k] += U_temp[i][k]; // U=U+dt*phi(U+dt/2*phi(U))
            }
            temp += U[i][0];
            l1 += U[i][1];
            l2 += U[i][2];
            l3 += U[i][3];

            //std::cout<<U[i][0]<<" "<<U[i][1]<<" "<<U[i][2]<<" "<<U[i][3]<<std::endl;
        }


        if (steps == 0)
            rho_full = temp;

        std::cout << "t= " << t << " rho_total_err=" << 1. - temp / rho_full << " l_total_norm= " << sqrt(l1 * l1 + l2 * l2 + l3 * l3)
                  << " l_total = (" << l1 << "," << l2 << "," << l3 << ")" << std::endl;


        //std::cout<<std::endl;
        // std::cout <<t<<" "<<1.-temp/faces.size() << std::endl;

        t += dt;
        steps++;
    }

    double time()
    {

        return t;
    }

protected:
    /*double limiter(double r, double etha_plus, double etha_minus)
    { // classical Superbee limiter for irregular grids
        // CFL independent

        if (std::isnan(r))
        {
            r = 0;
        }

        return std::max(0., std::max(std::min(1., etha_minus * r), std::min(r, etha_plus)));
    };*/

    void res2d(double dt_here) // space step
    // takes data from U matrix
    {

        for (size_t i = 0; i < this->n_faces(); i++)
        {
            for (size_t k = 0; k < dim; k++)
                U[i][k] = 0;

            for (size_t j = 0; j < faces[i].size(); j++)
            {

                int neighboor_num = neighbors_edge[i][j];
                int j0 = std::find(neighbors_edge[neighboor_num].begin(), neighbors_edge[neighboor_num].end(), i) - neighbors_edge[neighboor_num].begin();
                int j01 = j0 + 1;
                int j1 = j + 1;

                if (j == (faces[i].size() - 1))
                    j1 = 0;

                if (j0 == (faces[neighboor_num].size() - 1))
                    j01 = 0;

                for (size_t k = 0; k < dim; k++)
                {

                    U[i][k] -= dt_here * ((vertices[faces[i][j]] - vertices[faces[i][j1]]).norm() / surface_area[i]) *
                               (flux_var_plus[i][j][k] + flux_var_minus[i][j][k]);

                    if (std::isnan((flux_var_plus[i][j][k] + flux_var_minus[i][j][k])))
                    {
                        std::cout <<"time: "<<t<<" face: "<< i << " edge: " << j << " NaN in flux detected!" << std::endl;
                    }

                }

                //std::cout<<" i= "<< i << " j= " << j << " flux: " << flux_var_minus[i][j][0] << " " << flux_var_minus[i][j][1] << " " << flux_var_minus[i][j][2] << " " << flux_var_minus[i][j][3] << std::endl;

                //int component =2;
                //std::cout << i << " " << j << " "<< neighboor_num<<" " << flux_var_plus[i][j][component] + flux_var_minus[i][j][component] << std::endl;
                //std::cout << i << " " << j << " "<< neighboor_num<<" " << flux_var_plus[neighboor_num][j0][component] + flux_var_minus[neighboor_num][j0][component] << std::endl;
                //std::cout << std::endl;
            }
        }
    };

    virtual std::vector<double> flux_star(std::vector<double> ul, std::vector<double> ur, int n_face, int n_edge) = 0;
    virtual std::vector<double> limiter(std::vector<double> u_r, int n_face, int n_edge) = 0;

private:
    void find_M()
    // computes value M = max (d phi/ d U1, - d phi/ d U2 )
    // requires recently updated flux_var_plus, flux_var_minus
    {
        double max, f1, f2;
        max = 0.001;
        for (size_t i = 0; i < this->n_faces(); i++)
        {
            for (size_t j = 0; j < faces[i].size(); j++)
            {
                for (size_t k = 0; k < dim; k++)
                {
                    f1 = flux_var_plus[i][j][k] / (U_plus[i][j][k] - U[i][k]);
                    f2 = -flux_var_minus[i][j][k] / (U_minus[i][j][k] - U[i][k]);

                    if (!std::isnan(f1) && !std::isinf(f1) && f1 > max)
                        max = f1;

                    if (!std::isnan(f2) && !std::isinf(f2) && f2 > max)
                        max = f2;
                }
            }
        }

        M = max;
    }
    void find_v_max()
    {
        double max;
        max = 0.001;
        vector3d<double> R_vec, vel, l_vec;
        for (size_t i = 0; i < this->n_faces(); i++)
        {

            l_vec[0] = U[i][1];
            l_vec[1] = U[i][2];
            l_vec[2] = U[i][3];
            vel = cross_product(face_centers[i], l_vec);
            if (vel.norm() > max)
                max = vel.norm();
        }

        max_vel = max;
    }

    void find_flux_var()
    {

        std::vector<double> phi_ii, phi_iji, phi_ijji, r1,r2;
        phi_ii.resize(dim);
        phi_iji.resize(dim);
        phi_ijji.resize(dim);

        for (int i = 0; i < this->n_faces(); ++i)
        {
            for (int j = 0; j < faces[i].size(); ++j)
            {
                
                phi_ijji = flux_star(U_plus[i][j], U_minus[i][j], i, j);

                //std::cout << i << " " << j << std::endl;
                //std::cout << U_plus[i][j][0] << " " << U_plus[i][j][1] << " " << U_plus[i][j][2] << " " << U_plus[i][j][3] << std::endl;

                //std::cout << U_minus[i][j][0] << " " << U_minus[i][j][1] << " " << U_minus[i][j][2] << " " << U_minus[i][j][3] << std::endl;
                //std::cout<<"HLLE_ijji: " << phi_ijji[0] << " " << phi_ijji[1] << " " << phi_ijji[2] << " " << phi_ijji[3] << std::endl;
                //std::cout<<std::endl;
                //phi_ii = flux_star(U[i], U[i], i, j);
                //phi_iji = flux_star(U_plus[i][j], U[i], i, j);

                flux_var_minus[i][j] = phi_ijji;



                /*if ((i == 4 && j == 1) || (i == 2 && j == 1))
                {

                    std::cout <<" i= "<< i << " j= " << j << std::endl;
                    std::cout<<"U_L:" << U_plus[i][j][0] << " " << U_plus[i][j][1] << " " << U_plus[i][j][2] << " " << U_plus[i][j][3] << std::endl;
                    std::cout<<"U_R:"  << U_minus[i][j][0] << " " << U_minus[i][j][1] << " " << U_minus[i][j][2] << " " << U_minus[i][j][3] << std::endl;

                    std::cout<<"HLLE_ijji: " << phi_ijji[0] << " " << phi_ijji[1] << " " << phi_ijji[2] << " " << phi_ijji[3] << std::endl;
                    
                    phi_ijji=flux_star(U_plus[i][j], U_minus[i][j], i, j);
                    std::cout <<"HLLE_ijji_manual: "<< phi_ijji[0] << " " << phi_ijji[1]
                    << " " << phi_ijji[2] << " " << phi_ijji[3] << std::endl;

                    std::cout <<flux_var_plus[i][j][0]+flux_var_minus[i][j][0]<<std::endl;
                    //std::cout << phi_ii[0] << " " << phi_iji[0] << " " << phi_ijji[0] << std::endl;
                    std::cout << std::endl;
                }*/


            }
        }

        // std::cout << std::endl;
    }

    void find_U_edges() // finding U_ij and U_ji
    {

        std::vector<double> pp, pm, lim;
        int j0, neighboor_num;
        for (size_t i = 0; i < this->n_faces(); ++i)
        {
            for (size_t j = 0; j < faces[i].size(); ++j)
            {

                neighboor_num = neighbors_edge[i][j];
                j0 = std::find(neighbors_edge[neighboor_num].begin(), neighbors_edge[neighboor_num].end(), i) - neighbors_edge[neighboor_num].begin();

                pp = p_plus(i, j);
                pm = p_minus(i, j);

                for (size_t k = 0; k < dim; k++)
                    pm[k] /= pp[k];

                lim = limiter(pm, i, j);
                for (size_t k = 0; k < dim; k++)
                {
                    // U_plus[i][j][k] = U[i][k] + pp[k] * limiter(pm[k] / pp[k], H_plus[i][j] / BM_dist[i][j], H_minus[i][j] / BM_dist[i][j]) * BM_dist[i][j];

                    U_plus[i][j][k] = U[i][k] + pp[k] * lim[k] * BM_dist[i][j];
                    U_minus[neighboor_num][j0][k] = U_plus[i][j][k];
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
};