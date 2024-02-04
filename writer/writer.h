//
// Created by egor on 1/23/24.
//

#ifndef LIBGRAN_WRITER_H
#define LIBGRAN_WRITER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <Eigen/Eigen>

void write_particles(const std::string & dir, std::vector<Eigen::Vector3d> const & x, std::vector<Eigen::Vector3d> const & theta, double r_part);

void write_spring_connectors(std::string const & dir,
                             double r_part,
                             std::vector<std::array<Eigen::Vector3d, 5>> const & spring_connectors,
                             std::vector<bool> const & enabled_contacts);

void write_necks(std::string const & dir, std::vector<Eigen::Vector3d> const & x, double r_part,
                 std::vector<bool> const & bonded_contacts);

void write_boundary_particles(std::string const & dir, std::vector<Eigen::Vector3d> const & x, double r_part,
    std::tuple<size_t, size_t, size_t> anchored_particles, size_t tip_particle);

#endif //LIBGRAN_WRITER_H
