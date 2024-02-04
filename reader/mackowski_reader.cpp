//
// Created by egor on 2/3/24.
//

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

#include "mackowski_reader.h"

// Function that loads an aggregate from a text file
// the aggregate needs to be generated by the program from Mackowski, 1995
std::vector<Eigen::Vector3d> load_mackowski_aggregate(std::string const & file_name) {
    std::ifstream ifs(file_name);

    if (!ifs.good()) {
        std::cerr << "Unable to read the aggregate file: " << file_name << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<Eigen::Vector3d> particles;
    std::string line;

    while (getline(ifs, line)) {
        if (!line.empty()) {
            std::istringstream oss(line);
            double _, x, y, z;
            oss >> _ >> x >> y >> z;
            particles.emplace_back(x, y, z);
        }
    }

    return particles;
}

// Iteratively removes overlap between adjacent particles
void remove_mackowski_overlap(std::vector<Eigen::Vector3d> & x, size_t num_iter, double critical_separation, double r_part) {
    // Build a neighbor list
    std::vector<std::pair<size_t, size_t>> neighbors;
    for (size_t i = 0; i < x.size() - 1; i ++) {
        for (size_t j = i + 1; j < x.size(); j ++) {
            if ((x[i] - x[j]).norm() - 2.0 * r_part <= critical_separation)
                neighbors.emplace_back(i, j);
        }
    }

    for (size_t i = 0; i < num_iter; i ++) {
        for (auto const & pair : neighbors) {
            auto & particleA = x[pair.first];
            auto & particleB = x[pair.second];

            Eigen::Vector3d n = (particleB - particleA).normalized();
            Eigen::Vector3d dist = ((particleB - particleA).norm() - 2.0) * n;

            // Move each particle by an equal amount so that they are in point-touch contact
            particleA += dist / 2.0;
            particleB -= dist / 2.0;
        }
    }
}