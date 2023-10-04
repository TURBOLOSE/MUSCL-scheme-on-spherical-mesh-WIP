#pragma once

#include "pmp/surface_mesh.h"
#include "sph_gen.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <cmath>

using namespace pmp;

class MUSCL_base
{

protected:
    SurfaceMesh mesh;
    std::vector<std::vector<double>> vertices; // coordinates of points
    std::vector<std::vector<int>> faces;       // indexes of vertices that make a face
    std::vector<std::vector<int>> neighbors;   // indexes of neighbor faces of a given face
    std::vector<std::vector<double>> face_centers;
    std::vector<std::vector<double>> normals;
    std::vector<std::vector<std::vector<int>>> faces_for_flux;
    // 4 indeces(2+ and 2-) of faces for + and - flux through edge for every edge of every face

public:
    MUSCL_base(SurfaceMesh mesh_in) : mesh(mesh_in)
    {
        // std::cout << "vertices: " << mesh.n_vertices() << std::endl;
        // std::cout << "edges: " << mesh.n_edges() << std::endl;
        // std::cout << "faces: " << mesh.n_faces() << std::endl;
        std::vector<double> r1, r2;

        vertices.resize(mesh.n_vertices());
        faces.resize(mesh.n_faces());
        face_centers.resize(mesh.n_faces());
        neighbors.resize(mesh.n_faces());
        normals.resize(mesh.n_faces());

        size_t i = 0; // getting vertices from SurfaceMesh
        auto points = mesh.get_vertex_property<Point>("v:point");
        for (auto vertice : mesh.vertices())
        {
            vertices[i].resize(points[vertice].size());

            for (size_t j = 0; j < points[vertice].size(); j++)
            {
                vertices[i][j] = points[vertice][j];
            }
            i++;
        }

        i = 0;

        for (auto face : mesh.faces()) // getting vertice numbers for faces
        {
            auto face_vertices = mesh.vertices(face);

            for (auto vertice = face_vertices.begin(); vertice != face_vertices.end(); ++vertice)
            {
                faces[i].push_back((*vertice).idx());
            }
            i++;
        }

        i = 0;
        for (auto face : mesh.faces()) // finding all neighbors for faces
        {
            auto face_vertices = mesh.vertices(face);

            for (auto vertice = face_vertices.begin(); vertice != face_vertices.end(); ++vertice)
            {
                auto faces_around = mesh.faces(*vertice);

                for (auto face_in : faces_around)
                {
                    neighbors[i].push_back((face_in).idx());
                }
            }
            i++;
        }

        i = 0;
        for (auto &&neighbor : neighbors)
        { // removing duplicates from neighbors

            std::sort(neighbor.begin(), neighbor.end());
            auto last = std::unique(neighbor.begin(), neighbor.end());
            neighbor.erase(last, neighbor.end());
            neighbor.erase(std::remove_if(neighbor.begin(), neighbor.end(), [&i](int x)
                                          { return x == i; }),
                           neighbor.end());
            i++;
        }

        for (size_t i = 0; i < mesh.n_faces(); i++)
        {
            face_centers[i].resize(3);
            for (size_t j = 0; j < 3; j++) // 3d points
            {
                for (size_t k = 0; k < faces[i].size(); k++)
                {
                    face_centers[i][j] += vertices[faces[i][k]][j];
                }

                face_centers[i][j] /= faces[i].size();
            }
        }

        r1.resize(3);
        r2.resize(3);

        for (size_t i = 0; i < mesh.n_faces(); i++)
        {

            normals[i].resize(3);

            std::transform(vertices[faces[i][0]].begin(), vertices[faces[i][0]].end(),
                           vertices[faces[i][1]].begin(), r1.begin(), std::minus<double>()); // r1=vertice_1-vertice_0

            std::transform(vertices[faces[i][1]].begin(), vertices[faces[i][1]].end(),
                           vertices[faces[i][2]].begin(), r2.begin(), std::minus<double>()); // r2=vertice_2-vertice_1

            normals[i].resize(3);

            normals[i][0] = r1[1] * r2[2] - r1[2] * r2[1]; // cross product
            normals[i][1] = -(r1[0] * r2[2] - r1[2] * r2[0]);
            normals[i][2] = r1[0] * r2[1] - r1[1] * r2[0];

            double norm = sqrt(normals[i][0] * normals[i][0] + normals[i][1] * normals[i][1] + normals[i][2] * normals[i][2]);

            std::transform(normals[i].begin(), normals[i].end(), normals[i].begin(),
                           std::bind(std::multiplies<double>(), std::placeholders::_1, 1 / norm));
        }

        process_mesh();
    };

    void process_mesh(){

    };

    double broken_distance(std::vector<double> a, std::vector<double> b)
    {
        // dim(a)=3; dim(b)=3
        // distance on a spherical mesh between points a and b

        double dist = 0;
        int start_face = point_in_face(a);
        int end_face = point_in_face(b);
        std::vector<double> bma, bma_face; // b-a and (b-a) projected to face
        std::vector<double> len;
        size_t n_edges=faces[start_face].size();


        bma.resize(3);
        bma_face.resize(3);

        for (size_t i = 0; i < 3; i++) // line vector
        {
            bma[i] = b[i] - a[i];
        }

        if (start_face == end_face)
        {
            dist = sqrt((b[0] - a[0]) * (b[0] - a[0]) + (b[1] - a[1]) * (b[1] - a[1]) + (b[2] - a[2]) * (b[2] - a[2]));
        }
        else
        {
            for (size_t i = 0; i < 3; i++) // line vector
            {
                bma_face[i] = normals[start_face][i] - bma[i];
            }


        
            for (size_t i = 0; i < n_edges; i++)
            {
                //WIP
            }
            



        }

        return dist;
    };

    int point_in_face(std::vector<double> r)
    { // dim(r)=3;
        // finds face in which point r is in
        int count = 0;
        int res;
        int flag = 0;
        double dotpr;
        std::vector<double> r_face;
        double eps = 1e-5; // abs precision

        for (auto face : faces)
        {
            dotpr = 0;
            r_face = vertices[face[0]];
            for (size_t i = 0; i < 3; i++)
            {
                dotpr += (r_face[i] - r[i]) * normals[count][i];
            }

            if (abs(dotpr) < eps) // compare with 0
            {
                res = count;
                flag += 1;
            }

            count++;
        }

        if (flag == 0)
            std::cout << "point_in_face: no face has been found" << std::endl;

        if (flag > 1)
            std::cout << "point_in_face: more than 1 face has been found" << std::endl;

        return res;
    };

    void print_vertices()
    {
        for (auto face : faces)
        {
            for (auto face_element : face)
            {
                std::cout << "(" << vertices[face_element][0] << "|" << vertices[face_element][1] << "|" << vertices[face_element][2] << ")";
            }
            std::cout << std::endl;
        }
    };

    void print_neighbors()
    {
        for (auto neighbor : neighbors)
        {
            std::cout << "|";
            for (auto neighbor_element : neighbor)
            {
                std::cout << neighbor_element << "|";
            }
            std::cout << std::endl;
        }
    };

    void print_normals()
    {
        for (auto normal : normals)
        {
            std::cout << "(" << normal[0] << "|" << normal[1] << "|" << normal[2] << ")";
            std::cout << std::endl;
        }
    };
};