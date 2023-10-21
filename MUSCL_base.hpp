#pragma once

#include "pmp/surface_mesh.h"
#include "sph_gen.h"
#include "vec3d.hpp"

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
    std::vector<vector3d<double>> vertices;  // coordinates of points
    std::vector<std::vector<int>> faces;     // indexes of vertices that make a face
    std::vector<std::vector<int>> neighbors; // indexes of neighbor faces of a given face
    std::vector<vector3d<double>> face_centers;
    std::vector<vector3d<double>> normals;
    std::vector<std::vector<std::vector<int>>> faces_for_flux;
    // 4 indeces(2+ and 2-) of faces for + and - flux through edge for every edge of every face

public:
    MUSCL_base(SurfaceMesh mesh_in) : mesh(mesh_in)
    {
        // std::cout << "vertices: " << mesh.n_vertices() << std::endl;
        // std::cout << "edges: " << mesh.n_edges() << std::endl;
        // std::cout << "faces: " << mesh.n_faces() << std::endl;
        vector3d<double> r1, r2;

        vertices.resize(mesh.n_vertices());
        faces.resize(mesh.n_faces());
        face_centers.resize(mesh.n_faces());
        neighbors.resize(mesh.n_faces());
        normals.resize(mesh.n_faces());

        size_t i = 0; // getting vertices from SurfaceMesh
        auto points = mesh.get_vertex_property<Point>("v:point");
        for (auto vertice : mesh.vertices())
        {
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
            // for (size_t j = 0; j < 3; j++) // 3d points
            //{
            for (size_t k = 0; k < faces[i].size(); k++)
            {
                face_centers[i] += vertices[faces[i][k]];
            }

            face_centers[i] /= faces[i].size();
            //}
        }

        for (size_t i = 0; i < mesh.n_faces(); i++)
        {
            std::transform(vertices[faces[i][0]].begin(), vertices[faces[i][0]].end(),
                           vertices[faces[i][1]].begin(), r1.begin(), std::minus<double>()); // r1=vertice_1-vertice_0

            std::transform(vertices[faces[i][1]].begin(), vertices[faces[i][1]].end(),
                           vertices[faces[i][2]].begin(), r2.begin(), std::minus<double>()); // r2=vertice_2-vertice_1

            normals[i] = cross_product(r1, r2);
            normals[i] /= normals[i].norm();
        }

        process_mesh();
    };

    void process_mesh()
    {
        vector3d<double> BM, B_nB, B_nB_face, BM_normal, B_left, B1B2, r, Hm; // vectors from B to M and from B_n to B
        double max_cos = -1;
        int left_face1, left_face2, right_face1, right_face2;

        // for (auto face : faces)
        //{
        int n_face = 1;
        std::vector<int> face = faces[n_face];

        for (size_t i = 0; i < face.size(); i++)
        {
            max_cos = -1;

            if (i == face.size() - 1)
            {
                BM = (vertices[face[i]] + vertices[face[0]]) / 2 - face_centers[n_face];
                r = (vertices[face[0]] - vertices[face[i]]);

                // vertices[face[i]].print();
                // vertices[face[0]].print();
            }
            else
            {
                BM = (vertices[face[i]] + vertices[face[i + 1]]) / 2 - face_centers[n_face];
                r = (vertices[face[i + 1]] - vertices[face[i]]);

                // vertices[face[i]].print();
                // vertices[face[i+1]].print();
            }

            for (auto neighboor : neighbors[n_face]) // find 1st most left element
            {

                B_nB = face_centers[n_face] - face_centers[neighboor];                   // not projected yet
                B_nB_face = B_nB - normals[n_face] * dot_product(B_nB, normals[n_face]); // projection

                double cos_etha = dot_product(BM, B_nB_face) / (BM.norm() * B_nB_face.norm());
                if (cos_etha > max_cos)
                {
                    max_cos = cos_etha;
                    left_face1 = neighboor;
                }
            }

            BM_normal = cross_product(BM, normals[n_face]);
            max_cos = -1;
            left_face2 = left_face1; // in case "if" doesent go off

            B_left = face_centers[left_face1] - face_centers[n_face]; // vector in direction of B_left_1

            vector3d<double> B_left_face;
            B_left_face = B_left - normals[n_face] * dot_product(B_left, normals[n_face]); // projection

            for (auto neighboor : neighbors[n_face]) // find second most left element
            {

                B_nB = face_centers[n_face] - face_centers[neighboor];
                B_nB_face = B_nB - normals[n_face] * dot_product(B_nB, normals[n_face]); // projection

                double cos_etha = dot_product(BM, B_nB_face) / (BM.norm() * B_nB_face.norm());

                double line_cross_check = (cross_product(B_nB_face, BM_normal)).norm() / (B_nB_face.norm() * BM_normal.norm()) *
                                          (cross_product(B_left_face, BM_normal).norm()) / (B_left_face.norm() * BM_normal.norm()) *
                                          dot_product(normals[n_face], cross_product(B_nB_face, BM_normal)) * dot_product(normals[n_face], cross_product(B_left_face, BM_normal));
                if ((cos_etha > max_cos) && (neighboor != left_face1) && (line_cross_check <= 0))
                {
                    max_cos = cos_etha;
                    left_face2 = neighboor;
                }
            }

            // std::cout << left_face1 << " " << left_face2 << std::endl;

            // finding H -- point of intersection between lines BiMij and line between face centers of lf1 and lf2

            // vector3d<double> H_0 = find_lines_intersection(face_centers[n_face], vertices[face[i]], BM, r);
        }
       


        vector3d<double> B1B2_proj1 = (face_centers[left_face2] - face_centers[left_face1]) -
                                      normals[left_face1] * dot_product((face_centers[left_face2] - face_centers[left_face1]), normals[left_face1]);

        vector3d<double> B1B2_proj2 = (face_centers[left_face1] - face_centers[left_face2]) -
                                      normals[left_face2] * dot_product((face_centers[left_face1] - face_centers[left_face2]), normals[left_face2]);

        vector3d<double> p1 = find_line_surf_intersection(face_centers[left_face1], B1B2_proj1,
                                                          face_centers[n_face], BM, normals[n_face]);

        vector3d<double> p2 = find_line_surf_intersection(face_centers[left_face1], B1B2_proj2,
                                                          face_centers[n_face], BM, normals[n_face]);


        if (is_on_surface(p1))
        {
            Hm = p1;
        }
        else if (is_on_surface(p2))
        {
            Hm = p2;
        }
        else
        {
            std::cout << "process_mesh: could not find H_minus" << std::endl;
        }

        double Hm_dist=broken_distance(face_centers[n_face], Hm);
        double B1B2_d = broken_distance(face_centers[left_face1], face_centers[left_face2]);


    };

    double broken_distance(vector3d<double> a, vector3d<double> b)
    {
        // broken distance between 2 points if broken line made out of 2 parts
        // dim(a)=3; dim(b)=3
        // distance on a spherical mesh between points a and b

        double dist = 0;
        int maxiter = 4;
        int iter = 0;

        int start_face = point_in_face(a)[0];
        int current_face = point_in_face(a)[0];
        int end_face = point_in_face(b)[0];
        vector3d<double> bma, intersection, intersection_prev;
        intersection_prev = a;

        while ((current_face != end_face) && (iter < maxiter))
        {

            bma = b - a; // line vector

            intersection = broken_distance_base(intersection_prev, b, bma, current_face);

            dist += (intersection - intersection_prev).norm();

            if (point_in_face(intersection)[0] == current_face)
            {
                current_face = point_in_face(intersection)[1];
            }
            else
            {
                current_face = point_in_face(intersection)[0];
            }

            intersection_prev = intersection;
            iter++;
        }

        intersection = intersection_prev;

        dist += (b - intersection).norm();
        return dist;
    }

    vector3d<double> broken_distance_base(vector3d<double> a, vector3d<double> b, vector3d<double> bma, int face_num)
    {
        // returns intersection point from a to b among projection of bma inside  face with index face_num

        double dist = 0;
        double prec = 1e-9;
        double min_etha, sinetha, cosetha;
        int edge1, edge2;
        signed int sign;
        signed int sign_prev = 0;
        vector3d<double> bma_face, r, r_edge; // b-a and (b-a) projected to face, r-- edge vector
        vector3d<double> intersection;        // rhs of linear system and intersection points
        Eigen::Matrix<double, 3, 3> dist_matrix;
        Eigen::Matrix<double, 3, 1> rhs;

        size_t n_edges = faces[face_num].size();

        bma_face = bma - normals[face_num] * dot_product(bma, normals[face_num]);
        bma_face /= bma_face.norm();

        for (size_t i = 0; i < n_edges; i++) // finding edge that our projected vector crosses
        {

            if (i < n_edges - 1)
            {
                r = (vertices[faces[face_num][i]] + vertices[faces[face_num][i + 1]]) / 2 - a;
            }
            else
            {
                r = (vertices[faces[face_num][i]] + vertices[faces[face_num][0]]) / 2 - a;
            }

            if (abs(dot_product(cross_product(bma_face, r), normals[face_num])) < prec)
            {

                sinetha = 0;
            }
            else
            {

                sinetha = (cross_product(bma_face, r)).norm() / (bma_face.norm() * r.norm()) *
                          dot_product(cross_product(bma_face, r), normals[face_num]) / abs(dot_product(cross_product(bma_face, r), normals[face_num]));
            }

            cosetha = dot_product(bma_face, r) / (bma_face.norm() * r.norm());

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

            // std::cout << dot_product(cross_product(bma_face, r), normals[face_num]) << std::endl;

            // cross_product(bma_face, r).print();
            // normals[face_num].print();
        }

        r_edge = vertices[faces[face_num][edge2]] - vertices[faces[face_num][edge1]];

        intersection = find_lines_intersection(a, vertices[faces[face_num][edge1]], bma_face, r_edge);

        return intersection;
    };

    vector3d<double> find_lines_intersection(vector3d<double> x1, vector3d<double> x2,
                                             vector3d<double> n1, vector3d<double> n2)
    {
        // x1 and x2 -- starting points of lines
        // n1 and n2 -- vectors of lines
        // method returns coordinates of point of intersection
        vector3d<double> res;
        double t0, t1, t, tau;
        double prec = 1e-7; // compare with 0

        // std::cout<<"here->"<<abs(n2[0]) <<" "<<abs(n2[1])<<" "<< abs(n2[2])<<std::endl;

        if (abs(n2[0]) > prec)
        {

            // std::cout<<abs(n1[1] - n2[1] * n1[0] / n2[0])<<" "<<abs(n1[2] - n2[2] * n1[0] / n2[0])<<std::endl;

            if (abs(n1[1] - n2[1] * n1[0] / n2[0]) > prec)
            {
                t0 = (x2[1] - x1[1] + n2[1] * (x1[0] - x2[0]) / n2[0]) / (n1[1] - n2[1] * n1[0] / n2[0]);
            }
            if (abs(n1[2] - n2[2] * n1[0] / n2[0]) > prec)
            {
                t1 = (x2[2] - x1[2] + n2[2] * (x1[0] - x2[0]) / n2[0]) / (n1[2] - n2[2] * n1[0] / n2[0]);
            }
            if ((abs(n1[1] - n2[1] * n1[0] / n2[0]) < prec) && (abs(n1[2] - n2[2] * n1[0] / n2[0]) < prec))
            {
                std::cout << "find_lines_intersection error: lines dont intersect" << std::endl;
                t0 = 0;
            }

            double err0 = x1[2] + n1[2] * t0 - x2[2] - n2[2] * (x1[0] - x2[0] + n1[0] * t0) / n2[0];
            double err1 = x1[1] + n1[1] * t0 - x2[1] - n2[1] * (x1[0] - x2[0] + n1[0] * t1) / n2[0];

            if (abs(err0) < abs(err1))
            {
                t = t0;
            }
            else
            {
                t = t1;
            }

            tau = (x1[0] - x2[0] + n1[0] * t) / n2[0];
        }
        else if (abs(n2[1]) > prec)
        {
            if (abs(n1[0] - n2[0] * n1[1] / n2[1]) > prec)
            {
                t0 = (x2[0] - x1[0] + n2[0] * (x1[1] - x2[1]) / n2[1]) / (n1[0] - n2[0] * n1[1] / n2[1]);
            }
            if (abs(n1[2] - n2[2] * n1[1] / n2[1]) > prec)
            {
                t0 = (x2[2] - x1[2] + n2[2] * (x1[1] - x2[1]) / n2[1]) / (n1[2] - n2[2] * n1[1] / n2[1]);
            }
            if (abs(n1[2] - n2[2] * n1[1] / n2[1]) < prec && abs(n1[0] - n2[0] * n1[1] / n2[1]) < prec)
            {
                std::cout << "find_lines_intersection error: lines dont intersect" << std::endl;
                t0 = 0;
                tau = 0;
            }

            double err0 = x1[2] + n1[2] * t0 - x2[2] - n2[2] * (x1[1] - x2[1] + n1[1] * t0) / n2[1];
            double err1 = x1[0] + n1[0] * t0 - x2[0] - n2[0] * (x1[1] - x2[1] + n1[1] * t1) / n2[1];

            if (abs(err0) < abs(err1))
            {
                t = t0;
            }
            else
            {
                t = t1;
            }
            tau = (x1[1] - x2[1] + n1[1] * t) / n2[1];
        }
        else if (abs(n2[2]) > prec)
        {
            if (abs(n1[0] - n2[0] * n1[2] / n2[2]) > prec)
            {
                t0 = (x2[0] - x1[0] + n2[0] * (x1[2] - x2[2]) / n2[2]) / (n1[0] - n2[0] * n1[2] / n2[2]);
            }
            if (abs(n1[1] - n2[1] * n1[2] / n2[2]) > prec)
            {
                t0 = (x2[1] - x1[1] + n2[1] * (x1[2] - x2[2]) / n2[2]) / (n1[1] - n2[1] * n1[2] / n2[2]);
            }
            if (abs(n1[1] - n2[1] * n1[2] / n2[2]) < prec && abs(n1[0] - n2[0] * n1[2] / n2[2]) < prec)
            {
                std::cout << "find_lines_intersection error: lines dont intersect" << std::endl;
                t = 0;
                tau = 0;
            }

            double err0 = x1[1] + n1[1] * t0 - x2[1] - n2[1] * (x1[2] - x2[2] + n1[2] * t0) / n2[2];
            double err1 = x1[0] + n1[0] * t0 - x2[0] - n2[0] * (x1[2] - x2[2] + n1[2] * t1) / n2[2];

            if (abs(err0) < abs(err1))
            {
                t = t0;
            }
            else
            {
                t = t1;
            }
            tau = (x1[2] - x2[2] + n1[2] * t) / n2[2];
        }
        else
        {
            std::cout << "find_lines_intersection error: lines dont intersect" << std::endl;
            t = 0;
            tau = 0;
        }

        for (size_t i = 0; i < 3; i++)
        {
            res[i] = x2[i] + tau * n2[i];
        }

        return res;
    }

    vector3d<double> find_line_surf_intersection(vector3d<double> xl, vector3d<double> nl,
                                                 vector3d<double> xs, vector3d<double> vs1, vector3d<double> vs2)
    {
        // xl -- line starting point, nl -- line vector
        // xs -- surface starting point,  vs1,2 -- any vector in surface(non parallel)

        vector3d<double> res;
        double A, B, C, D, t;
        double eps = 1e-7;
        if (cross_product(vs1, vs2).norm() < eps)
        {
            std::cout << "find_line_surf_intersection:parallel vectors cannot define surface" << std::endl;
        }

        A = vs1[1] * vs2[2] - vs1[2] * vs2[1];
        B = vs1[2] * vs2[0] - vs1[0] * vs2[2];
        C = vs1[0] * vs2[1] - vs1[1] * vs2[0];
        D = -xs[0] * A - xs[1] * B - xs[2] * C;

        t = -(D + A * xl[0] + B * xl[1] + C * xl[2]) / (A * nl[0] + B * nl[1] + C * nl[2]);

        res = xl + nl * t;

        return res;
    }

    std::vector<int> point_in_face(vector3d<double> r)
    { // dim(r)=3;
        // finds face in which point r is in
        int count = 0;
        std::vector<int> res;
        int flag = 0;
        double dotpr;
        vector3d<double> r_face;
        double eps = 1e-5; // abs precision

        for (auto face : faces)
        {
            r_face = face_centers[count];
            dotpr = dot_product((r_face - r), normals[count]);
            if (abs(dotpr) < eps) // compare with 0
            {
                res.push_back(count);
                flag += 1;
            }

            count++;
        }
        if (flag == 0)
            std::cout << "point_in_face: no face has been found" << std::endl;

        return res;
    };

    bool is_on_surface(vector3d<double> r)
    {
        int count = 0;
        bool res = false;
        double dotpr;
        vector3d<double> r_face;
        double eps = 1e-5; // abs precision

        for (auto face : faces)
        {
            r_face = face_centers[count];
            dotpr = dot_product((r_face - r), normals[count]);
            if (abs(dotpr) < eps) // compare with 0
            {
                res = true;
            }

            count++;
        }

        return res;
    }

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

protected:
    vector3d<double> cross_product(vector3d<double> a, vector3d<double> b)
    {
        vector3d<double> res;
        res[0] = a[1] * b[2] - a[2] * b[1]; // cross prod
        res[1] = a[2] * b[0] - a[0] * b[2];
        res[2] = a[0] * b[1] - a[1] * b[0];

        return res;
    }

    double dot_product(vector3d<double> a, vector3d<double> b)
    {
        double res;
        res = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

        return res;
    }
};