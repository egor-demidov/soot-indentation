//
// Created by egor on 2/3/24.
//

#include <vector>
#include <iostream>
#include <random>

#include <Eigen/Eigen>

#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/sinter_bridge/alt_sinter_bridge.h>
#include <libgran/granular_system/granular_system.h>

#include "../reader/mackowski_reader.h"
#include "../writer/writer.h"
#include "../energy/compute_energy.h"
#include "../afm/afm_step_handler.h"
#include "../afm/anchor_plane.h"

// Van der Waals attraction model
using vdw_functor_t = hamaker_functor<Eigen::Vector3d, double>;
// Non-bonded contact model
using non_bonded_contact_functor = contact_force_functor<Eigen::Vector3d, double>;
// Combined bonded and non-bonded contact models
using bonded_contact_functor_t = alt_sinter_functor<Eigen::Vector3d, double>;
// We will wrap the generic step handler into the afm step handler
template <typename field_container_t, typename field_value_t>
using afm_step_handler_t = afm_step_handler<field_container_t, field_value_t, rotational_step_handler, double>;
// Container for binary force functors
using binary_functors_t = binary_force_functor_container<Eigen::Vector3d, double, vdw_functor_t, bonded_contact_functor_t>;
// Container for unary force functors (empty)
using unary_functors_t = unary_force_functor_container<Eigen::Vector3d, double>;
// Assemble the granular system representation
using granular_system_t = granular_system<Eigen::Vector3d, double, rotational_velocity_verlet_half,
    afm_step_handler_t, binary_functors_t, unary_functors_t>;

int main() {
    // General simulation parameters
    const double dt = 1e-14;
    const double t_tot = 3.0e-7;
    const auto n_steps = size_t(t_tot / dt);
    const size_t n_dumps = 300;
    const size_t dump_period = n_steps / n_dumps;
    const size_t n_thermo_dumps = 10000;
    const size_t thermo_dump_period = n_steps / n_thermo_dumps;

    // General parameters
    const double rho = 1700.0;
    const double r_part = 14.0e-9;
    const double mass = 4.0 / 3.0 * M_PI * pow(r_part, 3.0) * rho;
    const double inertia = 2.0 / 5.0 * mass * pow(r_part, 2.0);

    // Parameters for the non-bonded contact model
    const double k = 10000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the bonded contact model
    const double k_bond = 1000.0;
    const double gamma_n_bond = 1.25e-8;
    const double gamma_t_bond = 0.2 * gamma_n;
    const double gamma_r_bond = 0.05 * gamma_n;
    const double gamma_o_bond = 0.05 * gamma_n;
    const double d_crit = 1.0e-9;

    // Parameters for the Van der Waals model
    const double A = 1.0e-20;
    const double h0 = 1.0e-9;

    // Load the aggregate from file
    auto x0 = load_mackowski_aggregate("../aggregates/size_50/aggregate_1.txt");
    std::cout << "Loaded " << x0.size() << " particles from the aggregate file\n";
    // Pre-process the loaded aggregate by removing overlap
    remove_mackowski_overlap(x0, 100, 0.01, 1.0);
    // Re-size the monomers
    for (auto & point : x0) {
        point *= r_part;
    }

    // Pick points that will be anchored to the substrate
    auto anchor_points = find_anchor_plane_points(x0, 0);
    // Pick the tip point
    size_t tip_point = 45;
    // Solve for the contact plane normal
    Eigen::Matrix3d mat;
    mat.row(0) = x0[std::get<0>(anchor_points)];
    mat.row(1) = x0[std::get<1>(anchor_points)];
    mat.row(2) = x0[std::get<2>(anchor_points)];

    Eigen::Vector3d c = {1, 1, 1};
    Eigen::Vector3d normal = (mat.inverse() * c).normalized();

    // Make sure normal is pointing inward
    Eigen::Vector3d distance = x0[tip_point] - x0[std::get<0>(anchor_points)];
    std::cout << "Anchor points are: " <<
        std::get<0>(anchor_points) << " " <<
        std::get<1>(anchor_points) << " " <<
        std::get<2>(anchor_points) << "\n";

    if (distance.dot(normal) < 0.0)
        normal = -normal;

    // Initialize the remaining initial conditions
    std::vector<Eigen::Vector3d> v0, theta0, omega0;
    v0.resize(x0.size());
    theta0.resize(x0.size());
    omega0.resize(x0.size());
    std::fill(v0.begin(), v0.end(), Eigen::Vector3d::Zero());
    std::fill(theta0.begin(), theta0.end(), Eigen::Vector3d::Zero());
    std::fill(omega0.begin(), omega0.end(), Eigen::Vector3d::Zero());

    // Prescribe the retraction velocity to the AFM tip
    v0[tip_point] = normal * 0.1;

    // Initialize the force buffer
    std::vector<double> force_buffer(n_thermo_dumps);
    std::fill(force_buffer.begin(), force_buffer.end(), 0.0);

    rotational_step_handler<std::vector<Eigen::Vector3d>, Eigen::Vector3d> base_step_handler;
    afm_step_handler_t<std::vector<Eigen::Vector3d>, Eigen::Vector3d> step_handler(anchor_points,
               tip_point, x0.size(), thermo_dump_period, base_step_handler, mass, force_buffer.begin());

    non_bonded_contact_functor non_bonded_model(x0.size(),
                                                k, gamma_n, k, gamma_t, mu, phi, k, gamma_r, mu_o, phi, k, gamma_o, mu_o, phi,
                                                r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0);

    bonded_contact_functor_t sinter_model(x0.size(), x0,
                                          k_bond, gamma_n_bond, k_bond, gamma_t_bond, k_bond, gamma_r_bond, k_bond,
                                          gamma_o_bond, r_part, mass, inertia, dt, Eigen::Vector3d::Zero(), 0.0, d_crit, non_bonded_model);

    vdw_functor_t hamaker_model(A, h0,
                                      r_part, mass, Eigen::Vector3d::Zero(), 0.0);

    binary_functors_t binary_forces{hamaker_model, sinter_model};
    unary_functors_t unary_forces{};

    granular_system_t system(x0,
         v0, theta0, omega0, 0.0, Eigen::Vector3d::Zero(), 0.0,
            step_handler, binary_forces, unary_forces);

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump #" << n / dump_period << std::endl;
            write_particles("run", system.get_x(), system.get_theta(), r_part);
            write_boundary_particles("run", system.get_x(), r_part, anchor_points, tip_point);
            write_necks("run", system.get_x(), r_part, sinter_model.bonded_contacts);

            if (n == 2 * dump_period) {
                // Break a bond
//                size_t bond_number = 35; // This is a looped branch
                size_t bond_number = 27;
                size_t bond_counter = 0;
                for (size_t i = 0; i < x0.size() - 1; i ++) {
                    for (size_t j = i + 1; j < x0.size(); j ++) {
                        if (sinter_model.bonded_contacts[x0.size() * i + j]) {
                            // This is a bonded contact
                            if (bond_counter == bond_number) {
                                sinter_model.bonded_contacts[x0.size() * i + j] = false;
                                sinter_model.bonded_contacts[x0.size() * j + i] = false;
                            }
                            bond_counter ++;
                        }
                    }
                }
            }
        }
        system.do_step(dt);
    }

    return 0;
}
