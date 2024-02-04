//
// Created by egor on 2/4/24.
//

#include <algorithm>
#include <tuple>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

#include "anchor_plane.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Kernel::Point_3 Point_3;
typedef std::pair<Point_3, size_t> Point_with_index;
typedef CGAL::First_of_pair_property_map<Point_with_index> Property_map;
typedef CGAL::Extreme_points_traits_adapter_3<Property_map, CGAL::Convex_hull_traits_3<Kernel>> CHT;
typedef CGAL::Surface_mesh<Point_with_index> Surface_mesh;

std::tuple<size_t, size_t, size_t> find_anchor_plane_points(std::vector<Eigen::Vector3d> const & x, size_t face_num) {
    std::vector<Point_with_index> points;
    points.resize(x.size());

    for (size_t i = 0; i < x.size(); i ++) {
        points.emplace_back(Point_3{x[i][0], x[i][1], x[i][2]}, i);
    }

    Surface_mesh sm;
    CGAL::convex_hull_3(points.begin(), points.end(), sm, CHT(Property_map()));

    Surface_mesh::Face_index face_index(face_num);
    CGAL::Vertex_around_face_circulator<Surface_mesh> vcirc(sm.halfedge(face_index), sm);
    std::array<size_t, 3> vertex_indices {0};

    for (size_t i = 0; i < 3; i ++) {
        vertex_indices[i] = size_t(*vcirc++);
    }

    size_t a = sm.point(*(sm.vertices_begin() + vertex_indices[0])).second;
    size_t b = sm.point(*(sm.vertices_begin() + vertex_indices[1])).second;
    size_t c = sm.point(*(sm.vertices_begin() + vertex_indices[2])).second;

    return {a, b, c};
}
