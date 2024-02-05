//
// Created by egor on 2/3/24.
//

#include <vector>
#include <iostream>
#include <random>

#include <cfenv>

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

int main(int argc, char ** argv) {
    feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);

    if (argc < 7) {
        std::cerr << "Not enough arguments provided\n" <<
            "Required args (6): face number to use as anchor plane, tip point, bond to break, aggregate file path, dump directory, data file path\n";
        return EXIT_FAILURE;
    }

    // Parse the CL args
    size_t face_num = std::stoul(argv[1]);
    size_t tip_point = std::stoul(argv[2]);
    size_t bond_to_break = std::stoul(argv[3]);
    const std::string aggregate_file_path = argv[4];
    const std::string dump_dir = argv[5];
    const std::string data_file = argv[6];

    // General simulation parameters
    const double dt = 1e-14;
    const double t_tot = 1.5e-6;
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
    const double k = 10'000.0;
    const double gamma_n = 5.0e-9;
    const double mu = 1.0;
    const double phi = 1.0;
    const double mu_o = 0.1;
    const double gamma_t = 0.2 * gamma_n;
    const double gamma_r = 0.05 * gamma_n;
    const double gamma_o = 0.05 * gamma_n;

    // Parameters for the bonded contact model
    const double k_bond = 10'000.0;
    const double gamma_n_bond = 2.0*sqrt(2.0*mass*k_bond);
    const double gamma_t_bond = 0.2 * gamma_n;
    const double gamma_r_bond = 0.05 * gamma_n;
    const double gamma_o_bond = 0.05 * gamma_n;
    const double d_crit = 1.0e-9;

    // Parameters for the Van der Waals model
    const double A = 1.0e-20;
    const double h0 = 1.0e-9;

    // Load the aggregate from file
    auto x0 = load_mackowski_aggregate(aggregate_file_path);
    std::cout << "Loaded " << x0.size() << " particles from the aggregate file\n";
    // Pre-process the loaded aggregate by removing overlap
    remove_mackowski_overlap(x0, 100, 0.01, 1.0);
    // Re-size the monomers
    for (auto & point : x0) {
        point *= r_part;
    }

    // Pick points that will be anchored to the substrate
    auto anchor_points = find_anchor_plane_points(x0, face_num);
    // Pick the tip point
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
    std::vector<double> force_buffer(n_thermo_dumps+1);
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

    // Buffers for thermo data
    std::vector<double> t_span, ke_trs_span, ke_rot_span, ke_tot_span, lm_span, am_span, lm_span_norm, am_span_norm;

    // Break a bond
    size_t bond_counter = 0;
    for (size_t i = 0; i < x0.size() - 1; i ++) {
        for (size_t j = i + 1; j < x0.size(); j ++) {
            if (sinter_model.bonded_contacts[x0.size() * i + j]) {
                // This is a bonded contact
                if (bond_counter == bond_to_break) {
                    sinter_model.bonded_contacts[x0.size() * i + j] = false;
                    sinter_model.bonded_contacts[x0.size() * j + i] = false;
                }
                bond_counter ++;
            }
        }
    }

    for (size_t n = 0; n < n_steps; n ++) {
        if (n % dump_period == 0) {
            std::cout << "Dump #" << n / dump_period << std::endl;
            write_particles(dump_dir, system.get_x(), system.get_theta(), r_part);
            write_boundary_particles(dump_dir, system.get_x(), r_part, anchor_points, tip_point);
            write_necks(dump_dir, system.get_x(), r_part, sinter_model.bonded_contacts);
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
    std::ofstream ofs(data_file);

    if (!ofs.good()) {
        std::cerr << "Unable to create a data file" << std::endl;
        return EXIT_FAILURE;
    }

    double ke_trs_max = *std::max_element(ke_trs_span.begin(), ke_trs_span.end());
    double ke_rot_max = *std::max_element(ke_rot_span.begin(), ke_rot_span.end());
    double ke_tot_max = *std::max_element(ke_tot_span.begin(), ke_tot_span.end());
    double lm_max = *std::max_element(lm_span.begin(), lm_span.end());
    double am_max = *std::max_element(am_span.begin(), am_span.end());

    lm_span_norm.resize(lm_span.size());
    am_span_norm.resize(am_span.size());

    // Normalize the buffers
    std::transform(ke_trs_span.begin(), ke_trs_span.end(), ke_trs_span.begin(), [ke_tot_max] (auto ke) {
        return ke / ke_tot_max;
    });
    std::transform(ke_rot_span.begin(), ke_rot_span.end(), ke_rot_span.begin(), [ke_tot_max] (auto ke) {
        return ke / ke_tot_max;
    });
    std::transform(lm_span.begin(), lm_span.end(), lm_span_norm.begin(), [lm_max] (auto lm) {
        return lm / lm_max;
    });
    std::transform(am_span.begin(), am_span.end(), am_span_norm.begin(), [am_max] (auto am) {
        return am / am_max;
    });

    // Write the data
    ofs << "t\tE_trs\tE_rot\tE_tot\tP\tL\tPnorm\tLnorm\tf\n";
    for (size_t i = 0; i < ke_trs_span.size(); i ++) {
        ofs << t_span[i] << "\t"
            << ke_trs_span[i] << "\t"
            << ke_rot_span[i] << "\t"
            << ke_tot_span[i] << "\t"
            << lm_span[i] << "\t"
            << am_span[i] << "\t"
            << lm_span_norm[i] << "\t"
            << am_span_norm[i] << "\t"
            << force_buffer[i] << "\n";
    }

    return 0;
}
