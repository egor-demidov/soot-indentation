//
// Created by egor on 2/3/24.
//

#ifndef SOOT_INDENTATION_MACKOWSKI_READER_H
#define SOOT_INDENTATION_MACKOWSKI_READER_H

#include <vector>

#include <Eigen/Eigen>

std::vector<Eigen::Vector3d> load_mackowski_aggregate(std::string const & file_name);
void remove_mackowski_overlap(std::vector<Eigen::Vector3d> & x, size_t num_iter,
                              double critical_separation, double r_part);

#endif //SOOT_INDENTATION_MACKOWSKI_READER_H
