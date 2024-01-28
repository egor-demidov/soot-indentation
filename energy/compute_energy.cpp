//
// Created by egor on 1/19/24.
//

#include "compute_energy.h"

double compute_rotational_kinetic_energy(std::vector<Eigen::Vector3d> const & omegas, double inertia) {
    double val = 0.0;
    for (auto const & omega : omegas) {
        val += inertia * omega.dot(omega);
    }
    return val / 2.0;
}

double compute_translational_kinetic_energy(std::vector<Eigen::Vector3d> const & vs, double m) {
    double val = 0.0;
    for (auto const & v : vs) {
        val += m * v.dot(v);
    }
    return val / 2.0;
}

double compute_linear_momentum(std::vector<Eigen::Vector3d> const & vs, double m) {
    Eigen::Vector3d val = Eigen::Vector3d::Zero();
    for (auto const & v : vs) {
        val += v;
    }
    return m * val.norm();
}

double compute_angular_momentum(std::vector<Eigen::Vector3d> const & xs, std::vector<Eigen::Vector3d> const & vs, double m, std::vector<Eigen::Vector3d> const & omegas, double inertia) {
    Eigen::Vector3d val = Eigen::Vector3d::Zero();
    for (auto const & omega : omegas) {
        val += inertia * omega;
    }
    for (size_t i = 0; i < xs.size(); i ++) {
        auto const & x = xs[i];
        auto const & v = vs[i];
        val += x.cross(m*v);
    }
    return val.norm();
}
