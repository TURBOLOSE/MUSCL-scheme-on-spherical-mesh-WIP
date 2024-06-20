#pragma once
#include "../MUSCL_base/MUSCL_base.hpp"
#include "../geometry/MUSCL_geometry.hpp"
#include "../physics/isothermal.hpp"

class MUSCL_HLLE : public isothermal
{

private:
    std::ofstream outfile, outfile_curl, outfile_b;

public:
    MUSCL_HLLE(SurfaceMesh mesh, std::vector<std::vector<double>> U_in, int dim, double gam)
        : isothermal(mesh, U_in, dim, gam){}


    std::vector<double> flux_star(std::vector<double> ul, std::vector<double> ur, int n_face, int n_edge)
    {

        std::vector<double> FL, FR, F, c_vel;
        F.resize(dim);

        double S_R, S_L;
        
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
        {   F = FR;
            std::cout<<"flux_star: check char vel"<<std::endl;
            stop_check=true;
        }


        return F;
    };

};
