#include <iostream>
#include "MUSCL_base.hpp"


using namespace pmp;

int main()
{

    double a=1;
    size_t N=3;

    SurfaceMesh mesh = square(a, N);
    auto points = mesh.get_vertex_property<Point>("v:point");

    MUSCL_base test(mesh);


    /*for (auto e : mesh.vertices())
    {
        std::cout<<points[e]<<std::endl;
    }*/
    
    

    return 0;
}