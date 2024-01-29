//
// Created by egor on 1/27/24.
//

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/no_unary_force/no_unary_force.h>
#include <libgran/granular_system/granular_system.h>

#include "../energy/compute_energy.h"
#include "../afm/afm_step_handler.h"
#include "../writer/writer.h"

template <typename field_container_t, typename field_value_t>
using afm_rotational_step_handler = afm_step_handler<field_container_t, field_value_t, rotational_step_handler>;

int main() {
    // General simulation parameters
    const double dt = 1e-13;
    const double t_tot = 5.0e-8;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 300;
    const size_t dump_period = n_steps / n_dumps;

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

    std::vector<Eigen::Vector3d> x0, v0, theta0, omega0;
    x0.emplace_back(Eigen::Vector3d::Zero());
    x0.emplace_back(0, 2.5*r_part, 0);

    v0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());
    theta0.resize(x0.size());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    omega0.resize(x0.size());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    // Create an instance of step_handler
    // Using field type Eigen::Vector3d with container std::vector
    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler_instance;

    // Fix the second particle in space
    std::vector<bool> fixed_particles = {false, true};
    afm_rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> custom_handler_instance(fixed_particles, step_handler_instance);

    // Create an instance of contact force model
    // Using field type Eigen::Vector3d with real type double
    contact_force_functor<Eigen::Vector3d, double> contact_force_model(x0.size(),
                                                                       k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                                                       r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    // Create an instance of Hamaker model
    hamaker_functor<Eigen::Vector3d, double> hamaker_model(A, h0,
                                                           r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    binary_force_functor_container<Eigen::Vector3d, double, contact_force_functor<Eigen::Vector3d, double>, hamaker_functor<Eigen::Vector3d, double>>
            binary_force_functors{contact_force_model, hamaker_model};

    no_unary_force_functor<Eigen::Vector3d, double> no_unary_force(Eigen::Vector3d::Zero());

    unary_force_functor_container<Eigen::Vector3d, double, no_unary_force_functor<Eigen::Vector3d, double>>
            unary_force_functors(no_unary_force);

    // Create an instance of granular_system using the contact force model
    // Using velocity Verlet integrator for rotational systems and a default
    // step handler for rotational systems
    granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
            afm_rotational_step_handler, binary_force_functor_container<Eigen::Vector3d, double, contact_force_functor<Eigen::Vector3d, double>, hamaker_functor<Eigen::Vector3d, double>>,
            unary_force_functor_container<Eigen::Vector3d, double, no_unary_force_functor<Eigen::Vector3d, double>>> system(x0,
                                                                                                v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0, custom_handler_instance, binary_force_functors, unary_force_functors);

    std::vector<double> separation_span, t_span, ke_trs_span, ke_rot_span, ke_tot_span, lm_span, am_span;

    bool has_collided = false;

    double separation0 = (system.get_x()[1] - system.get_x()[0]).norm() - 2.0 * r_part;
    double hamaker_U0 = -A / 6.0 * (2.0*pow(r_part, 2.0)/(4.0*r_part+separation0)/separation0 +
            2.0*pow(r_part, 2.0)/pow(2.0*r_part+separation0, 2.0) +
            log((4.0*r_part+separation0)*separation0/pow(2.0*r_part+separation0, 2.0)));

    std::ofstream ofs("../plots/verification/05_hamaker.dat");

    if (!ofs.good()) {
        std::cerr << "Unable to create a data file" << std::endl;
        return EXIT_FAILURE;
    }

    ofs << "D\tKE\tU\n";

    for (size_t n = 0; n < n_steps; n ++) {
        double separation = (system.get_x()[1] - system.get_x()[0]).norm() - 2.0 * r_part;
        if (separation <= h0) {
            has_collided = true;
            return 0;
        }
        if (n % dump_period == 0) {
            double hamaker_U = -A / 6.0 * (2.0*pow(r_part, 2.0)/(4.0*r_part+separation)/separation +
                                            2.0*pow(r_part, 2.0)/pow(2.0*r_part+separation, 2.0) +
                                            log((4.0*r_part+separation)*separation/pow(2.0*r_part+separation, 2.0)));

            double dU = hamaker_U - hamaker_U0;

            ofs << separation * 1e9 << "\t" << compute_translational_kinetic_energy(system.get_v(), mass) << "\t" << -dU << "\n";

        }
        system.do_step(dt);
    }

    return 0;
}
