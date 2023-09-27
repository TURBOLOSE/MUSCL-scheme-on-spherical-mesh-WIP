#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "pmp/surface_mesh.h"
#include "sph_gen.h"

using namespace pmp;

class MUSCL_base
{

protected:
    SurfaceMesh mesh;
    std::vector<std::vector<double>> vertices; // coordinates of points
    std::vector<std::vector<int>> faces;       // indexes of vertices that make a face
    std::vector<std::vector<int>> neighboors;  // indexes of neighboor faces of a given face

public:
    MUSCL_base(SurfaceMesh mesh_in) : mesh(mesh_in)
    {
        std::cout << "vertices: " << mesh.n_vertices() << std::endl;
        std::cout << "edges: " << mesh.n_edges() << std::endl;
        std::cout << "faces: " << mesh.n_faces() << std::endl;

        vertices.resize(mesh.n_vertices());
        faces.resize(mesh.n_faces());
        neighboors.resize(mesh.n_faces());

        size_t i = 0; // getting vertices from SurfaceMesh
        auto points = mesh.get_vertex_property<Point>("v:point");
        for (auto vertice : mesh.vertices())
        {
            vertices[i].resize(points[vertice].size());

            for (size_t j = 0; j < points[vertice].size(); j++)
            {
                vertices[i][j] = points[vertice][j];
            }
        }

        i = 0;

        for (auto face : mesh.faces())
        {
            auto face_vertices = mesh.vertices(face);

            for (auto vertice = face_vertices.begin(); vertice != face_vertices.end(); ++vertice)
            {
                faces[i].push_back((*vertice).idx());
            }
            i++;
        }
  

        //============================================================================
        /*
        Point face_center(0, 0, 0);
        int points_in_face = 0;

        for (auto i = v.begin(); i != v.end(); ++i)
        {
            face_center += points[*i];
            points_in_face++;
            std::cout << points[*i] << std::endl;
        }

        face_center /= points_in_face;
        std::cout << face_center << std::endl;

        auto edge0 = mesh.from_vertex(mesh.halfedge(*face0));
        auto edge1 = mesh.from_vertex(mesh.halfedge(*face0));

        std::cout<<points[edge0]<<std::endl;
        std::cout<<points[edge1]<<std::endl;*/
    };
};