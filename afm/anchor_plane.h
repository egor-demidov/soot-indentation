//
// Created by egor on 2/4/24.
//

#ifndef SOOT_INDENTATION_ANCHOR_PLANE_H
#define SOOT_INDENTATION_ANCHOR_PLANE_H

#include <vector>
#include <tuple>
#include <cstddef>

#include <Eigen/Eigen>

std::tuple<size_t, size_t, size_t> find_anchor_plane_points(std::vector<Eigen::Vector3d> const & x, size_t face_num);

#endif //SOOT_INDENTATION_ANCHOR_PLANE_H
