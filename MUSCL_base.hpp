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
    { // broken distance between 2 points if broken line made out of 2 parts

        // dim(a)=3; dim(b)=3
        // distance on a spherical mesh between points a and b

        double dist = 0;
        double min_etha;
        int start_face = point_in_face(a);
        int end_face = point_in_face(b);
        int edge1, edge2;
        signed int sign;
        signed int sign_prev = 0;
        std::vector<double> bma, bma_face, r, r_edge; // b-a and (b-a) projected to face, r-- edge vector
        std::vector<double> intersection;             // rhs of linear system and intersection points
        Eigen::Matrix<double, 3, 3> dist_matrix;
        Eigen::Matrix<double, 3, 1> rhs;

        size_t n_edges = faces[start_face].size();

        bma.resize(3);
        // rhs.resize(3);
        intersection.resize(3);
        bma_face.resize(3);
        r.resize(3);
        r_edge.resize(3);

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
            for (size_t i = 0; i < 3; i++) // projected line vector
            {
                bma_face[i] = bma[i] - normals[start_face][i] *
                                           (bma[0] * normals[start_face][0] +
                                            bma[1] * normals[start_face][1] +
                                            bma[2] * normals[start_face][2]);
            }

            double bma_face_norm = sqrt(bma_face[0] * bma_face[0] + bma_face[1] * bma_face[1] + bma_face[2] * bma_face[2]);

            for (size_t i = 0; i < 3; i++) // projected line vector
            {
                bma_face[i] /= bma_face_norm;
            }

            for (size_t i = 0; i < n_edges; i++) // finding edge that our projected vector crosses
            {

                for (size_t j = 0; j < 3; j++) // vector from a to edge centers
                {

                    if (i < n_edges - 1)
                    {
                        r[j] = -a[j] + (vertices[faces[start_face][i]][j] + vertices[faces[start_face][i + 1]][j]) / 2;
                    }
                    else
                    {
                        r[j] = -a[j] + (vertices[faces[start_face][i]][j] + vertices[faces[start_face][0]][j]) / 2;
                    }
                }

                double sinetha = sqrt((bma_face[1] * r[2] - bma_face[2] * r[1]) * (bma_face[1] * r[2] - bma_face[2] * r[1]) +
                                      (bma_face[2] * r[0] - bma_face[0] * r[2]) * (bma_face[2] * r[0] - bma_face[0] * r[2]) +
                                      (bma_face[0] * r[1] - bma_face[1] * r[0]) * (bma_face[0] * r[1] - bma_face[1] * r[0])) /
                                 (sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));

                double cosetha = (bma_face[0] * r[0] + bma_face[1] * r[1] + bma_face[2] * r[2]) /
                                 (sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));

                if (i == 0)
                {
                    min_etha = abs(std::atan2(sinetha, cosetha));
                    edge1 = 0;
                    edge2 = 1;
                }
                else if (abs(std::atan2(sinetha, cosetha)) < min_etha)
                {
                    min_etha = abs(std::atan2(sinetha, cosetha));
                    edge1 = i;
                    edge2 = i + 1;
                    if (edge2 == n_edges)
                    {
                        edge2 = 0;
                    }
                }
            }

            for (size_t i = 0; i < 3; i++)
            // edge vector
            {
                r_edge[i] = vertices[faces[start_face][edge2]][i] - vertices[faces[start_face][edge1]][i];
            }

            double r_edge_norm = sqrt(r_edge[0] * r_edge[0] + r_edge[1] * r_edge[1] + r_edge[2] * r_edge[2]);

            for (size_t i = 0; i < 3; i++)
            {
                r_edge[i] /= r_edge_norm;
            }

            intersection = find_lines_intersection(a, vertices[faces[start_face][edge1]], bma_face, r_edge);

            dist += sqrt((intersection[0] - a[0]) * (intersection[0] - a[0]) +
                         (intersection[1] - a[1]) * (intersection[1] - a[1]) +
                         (intersection[2] - a[2]) * (intersection[2] - a[2]));

            // dist+=broken_distance(intersection,b);
        }

        return dist;
    };

    std::vector<double> find_lines_intersection(std::vector<double> x1, std::vector<double> x2,
                                                std::vector<double> n1, std::vector<double> n2)
    {
        // x1 and x2 -- starting points of lines
        // n1 and n2 -- vectors of lines
        // method returns coordinates of point of intersection
        std::vector<double> res;
        res.resize(3);
        double t0;
        double prec = 1e-7; // compare with 0

        if (abs(n2[0]) > prec)
        {
            if (abs(n1[1] - n2[1] * n1[0] / n2[0]) > prec)
            {
                t0 = (x2[1] - x1[1] + n2[1] * (x1[0] - x2[0]) / n2[0]) / (n1[1] - n2[1] * n1[0] / n2[0]);
            }
            else if (abs(n1[2] - n2[2] * n1[0] / n2[0]) > prec)
            {
                t0 = (x2[2] - x1[2] + n2[2] * (x1[0] - x2[0]) / n2[0]) / (n1[2] - n2[2] * n1[0] / n2[0]);
            }
            else
            {
                std::cout << "find_lines_intersection error: lines dont intersect" << std::endl;
                t0 = 0;
            }
        }
        else if (abs(n2[1]) > prec)
        {
            if (abs(n1[0] - n2[0] * n1[1] / n2[1]) > prec)
            {
                t0 = (x2[0] - x1[0] + n2[0] * (x1[1] - x2[1]) / n2[1]) / (n1[0] - n2[0] * n1[1] / n2[1]);
            }
            else if (abs(n1[2] - n2[2] * n1[1] / n2[1]) > prec)
            {
                t0 = (x2[2] - x1[2] + n2[2] * (x1[1] - x2[1]) / n2[1]) / (n1[2] - n2[2] * n1[1] / n2[1]);
            }
            else
            {
                std::cout << "find_lines_intersection error: lines dont intersect" << std::endl;
                t0 = 0;
            }
        }
        else if (abs(n2[2]) > prec)
        {
            if (abs(n1[0] - n2[0] * n1[2] / n2[2]) > prec)
            {
                t0 = (x2[0] - x1[0] + n2[0] * (x1[2] - x2[2]) / n2[2]) / (n1[0] - n2[0] * n1[2] / n2[2]);
            }
            else if (abs(n1[1] - n2[1] * n1[2] / n2[2]) > prec)
            {
                t0 = (x2[1] - x1[1] + n2[1] * (x1[2] - x2[2]) / n2[2]) / (n1[1] - n2[1] * n1[2] / n2[2]);
            }
            else
            {
                std::cout << "find_lines_intersection error: lines dont intersect" << std::endl;
                t0 = 0;
            }
        }
        else
        {
            std::cout << "find_lines_intersection error: lines dont intersect" << std::endl;
            t0 = 0;
        }

        for (size_t i = 0; i < 3; i++)
        {
            res[i] = x1[i] + t0 * n1[i];
        }

        return res;
    }

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