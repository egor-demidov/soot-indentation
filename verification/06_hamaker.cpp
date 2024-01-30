//
// Created by egor on 1/30/24.
//

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

// Using Eigen for linear algebra
#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/granular_system/granular_system.h>

#include "../energy/compute_energy.h"
#include "../writer.h"

// "Assemble" the force models used in this simulation
using contact_force_functor_t = contact_force_functor<Eigen::Vector3d, double>; // Contact force
using vdw_force_dunctor_t = hamaker_functor<Eigen::Vector3d, double>; // Van der Waals force
using binary_force_container_t
    = binary_force_functor_container<Eigen::Vector3d, double, contact_force_functor_t, vdw_force_dunctor_t>; // Binary force container

using unary_force_container_t = unary_force_functor_container<Eigen::Vector3d, double>; // Unary force container (empty)

using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    rotational_step_handler, binary_force_container_t, unary_force_container_t>; // Granular system representation

int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 5.0e-9;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 300;
    const size_t dump_period = n_steps / n_dumps;
    const size_t n_thermo_dumps = 10000;
    const size_t thermo_dump_period = n_steps / n_thermo_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 1.4e-8;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

    // Parameters for the contact model
    const double k = 10000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the Van der Waals model
    const double A = 1.0e-20;
    const double h0 = 1.0e-9;

    // Initialize the particles
    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
    // A sphere
    const double r_sphere = 20.0 * r_part;
    const Eigen::Vector3d sphere_1_center = Eigen::Vector3d::Zero();
    const Eigen::Vector3d sphere_2_center = Eigen::Vector3d::Zero() + 2.05 * Eigen::Vector3d::UnitX() * r_sphere;

    const auto max_num_parts_per_dim = size_t(2.0 * r_sphere / r_part);

    // Create sphere 1
    for (size_t i = 0; i < max_num_parts_per_dim; i ++) {
        double x = sphere_1_center[0] - r_sphere + r_part + 2.0 * r_part * double(i);
        for (size_t j = 0; j < max_num_parts_per_dim; j ++) {
            double y = sphere_1_center[1] - r_sphere + r_part + 2.0 * r_part * double(j);
            for (size_t m = 0; m < max_num_parts_per_dim; m ++) {
                double z = sphere_1_center[2] - r_sphere + r_part + 2.0 * r_part * double(m);
                if ((Eigen::Vector3d{x, y, z} - sphere_1_center).norm() < r_sphere) {
                    x0.emplace_back(x, y, z);
                    v0.emplace_back(2, 0, 0);
                }
            }
        }
    }

    // Create sphere 2
    for (size_t i = 0; i < max_num_parts_per_dim; i ++) {
        double x = sphere_2_center[0] - r_sphere + r_part + 2.0 * r_part * double(i);
        for (size_t j = 0; j < max_num_parts_per_dim; j ++) {
            double y = sphere_2_center[1] - r_sphere + r_part + 2.0 * r_part * double(j);
            for (size_t m = 0; m < max_num_parts_per_dim; m ++) {
                double z = sphere_2_center[2] - r_sphere + r_part + 2.0 * r_part * double(m);
                if ((Eigen::Vector3d{x, y, z} - sphere_2_center).norm() < r_sphere) {
                    x0.emplace_back(x, y, z);
                    v0.emplace_back(-2, 0, 0);
                }
            }
        }
    }

    std::cout << "Number of particles: " << x0.size() << std::endl;

    // Initialize the remaining buffers
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    // Set the initial velocity of the above-plane particle

    // Create an instance of step_handler
    // Using field type Eigen::Vector3d with container std::vector
    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;

    // Create an instance of contact force model
    // Using field type Eigen::Vector3d with real type double
    contact_force_functor_t contact_force_model(x0.size(),
                                       k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                       r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of Hamaker model
    vdw_force_dunctor_t hamaker_model(A, h0,
                       r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    binary_force_container_t
            binary_force_functors{contact_force_model, hamaker_model};

    unary_force_container_t
            unary_force_functors;

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system_t system(x0,
                           v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
                           step_handler_instance, binary_force_functors, unary_force_functors);

    // Buffers for thermo data
    std::vector<double> t_span, ke_trs_span, ke_rot_span, ke_tot_span, lm_span, am_span;

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump #" << n / dump_period << std::endl;
            write_particles("run", system.get_x(), system.get_theta(), r_part);
        }
        if (n % thermo_dump_period == 0) {
            auto ke_trs = compute_translational_kinetic_energy(system.get_v(), mass);
            auto ke_rot = compute_rotational_kinetic_energy(system.get_omega(), inertia);
            auto ke_tot = ke_trs + ke_rot;
            auto lm = compute_linear_momentum(system.get_v(), mass);
            auto am = compute_angular_momentum(system.get_x(), system.get_v(), mass, system.get_omega(), inertia);
            t_span.emplace_back(double(n) * dt);
            ke_trs_span.emplace_back(ke_trs);
            ke_rot_span.emplace_back(ke_rot);
            ke_tot_span.emplace_back(ke_tot);
            lm_span.emplace_back(lm);
            am_span.emplace_back(am);
        }
        system.do_step(dt);
    }

    // Output stream for data
    std::ofstream ofs("../plots/verification/06_hamaker.dat");

    if (!ofs.good()) {
        std::cerr << "Unable to create a data file" << std::endl;
        return EXIT_FAILURE;
    }

    double ke_trs_max = *std::max_element(ke_trs_span.begin(), ke_trs_span.end());
    double ke_rot_max = *std::max_element(ke_rot_span.begin(), ke_rot_span.end());
    double ke_tot_max = *std::max_element(ke_tot_span.begin(), ke_tot_span.end());
    double lm_max = *std::max_element(lm_span.begin(), lm_span.end());
    double am_max = *std::max_element(am_span.begin(), am_span.end());

    // Normalize the buffers
    std::transform(ke_trs_span.begin(), ke_trs_span.end(), ke_trs_span.begin(), [ke_trs_max, ke_tot_max] (auto ke) {
        return ke / ke_tot_max;
    });
    std::transform(ke_rot_span.begin(), ke_rot_span.end(), ke_rot_span.begin(), [ke_rot_max, ke_tot_max] (auto ke) {
        return ke / ke_tot_max;
    });
    std::transform(lm_span.begin(), lm_span.end(), lm_span.begin(), [lm_max] (auto lm) {
        return lm / lm_max;
    });
    std::transform(am_span.begin(), am_span.end(), am_span.begin(), [am_max] (auto am) {
        return am / am_max;
    });

    // Write the data
    ofs << "t\tE_trs\tE_rot\tE_tot\tP\tL\n";
    for (size_t i = 0; i < ke_trs_span.size(); i ++) {
        ofs << t_span[i] << "\t"
            << ke_trs_span[i] << "\t"
            << ke_rot_span[i] << "\t"
            << ke_tot_span[i] << "\t"
            << lm_span[i] << "\t"
            << am_span[i] << "\n";
    }

    return 0;
}
