//
// Created by egor on 1/28/24.
//

#ifndef SOOT_INDENTATION_AFM_STEP_HANDLER_H
#define SOOT_INDENTATION_AFM_STEP_HANDLER_H

#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

#include <libtimestep/step_handler/step_handler.h>
#include <libtimestep/integrator/integrator.h>
#include <libtimestep/system/system.h>

template <typename real_t>
struct afm_filter_ode : public unary_system<real_t, real_t, forward_euler, step_handler, afm_filter_ode<real_t>> {
    afm_filter_ode(std::vector<real_t> f0, std::vector<real_t> df0, real_t omega0) :
            unary_system<real_t, real_t, forward_euler, step_handler, afm_filter_ode<real_t>>(std::move(f0),
                    std::move(df0), 0.0, 0.0, 0.0, *this, step_handler_instance), g(0.0), omega0(omega0) {}

    real_t compute_acceleration(size_t i,
                                std::vector<real_t> const & f,
                                std::vector<real_t> const & df,
                                real_t t [[maybe_unused]]) {

        return omega0*omega0*(g - f[i] - 2.0/omega0*df[i]);
    }

    void set_g(real_t new_g) {
        g = new_g;
    }

private:
    real_t g;
    const real_t omega0;
    step_handler<std::vector<real_t>, real_t> step_handler_instance;
};

// The user must ensure that the force buffer is sufficiently large to store all dumps
template <typename field_container_t,
        typename field_value_t,
        template <
            typename _field_container_t,
            typename _field_value_t>
        typename base_step_handler_t,
        typename real_t>
struct afm_step_handler {
    afm_step_handler(std::tuple<size_t, size_t, size_t> const & anchor_points,
                     size_t tip_point,
                     size_t n_part, size_t dump_period,
                     base_step_handler_t<field_container_t, field_value_t> base_step_handler,
                     real_t mass, typename std::vector<real_t>::iterator tip_force_buffer,
                     afm_filter_ode<real_t> & afm_filter_ode_instance, real_t dt) :

        dt(dt),
        filter(afm_filter_ode_instance),
        tip_point(tip_point), dump_period(dump_period), mass(mass),
        base_step_handler(std::move(base_step_handler)),
        tip_force_buffer(tip_force_buffer) {

        fixed_particles.resize(n_part);
        std::fill(fixed_particles.begin(), fixed_particles.end(), false);

        // Fix the anchor points
        fixed_particles[std::get<0>(anchor_points)] = true;
        fixed_particles[std::get<1>(anchor_points)] = true;
        fixed_particles[std::get<2>(anchor_points)] = true;

        // Fix the tip point
        fixed_particles[tip_point] = true;
    }

    // This method increments the specified value in the x buffer
    void increment_x(size_t n,                                                                              // index of the value to increment
                     field_value_t const & dx,                                                              // value of the position increment
                     typename field_container_t::iterator x_begin_itr,                                      // iterator pointing to the start of the x buffer
                     typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) {   // iterator pointing to the start of the alpha buffer

        base_step_handler.increment_x(n, dx, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    // This method increments the specified value in the v buffer
    void increment_v(size_t n,                                                                              // index of the value to increment
                     field_value_t const & dv,                                                              // value of the velocity increment
                     typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                     typename field_container_t::iterator v_begin_itr,                                      // iterator pointing to the start of the v buffer
                     typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                     typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                     typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                     typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) {   // iterator pointing to the start of the alpha buffer

        // If this is the tip, store the force acting on it
        if (n == tip_point) {
            if (counter % dump_period == 0) {
                *tip_force_buffer = filter.get_x()[0];
                tip_force_buffer ++;  // Increment the iterator
            }
            filter.set_g((*(a_begin_itr + n)).norm() * mass);
            filter.do_step(dt);
            counter ++;
        }

        if (!fixed_particles[n])
            base_step_handler.increment_v(n, dv, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    // This method increments the specified value in the theta buffer
    void increment_theta(size_t n,                                                                              // index of the value to increment
                         field_value_t const & dtheta,                                                          // value of the angle increment
                         typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                         typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                         typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                         typename field_container_t::iterator theta_begin_itr,                                  // iterator pointing to the start of the theta buffer
                         typename field_container_t::const_iterator omega_begin_itr [[maybe_unused]],           // iterator pointing to the start of the omega buffer
                         typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) {   // iterator pointing to the start of the alpha buffer

        base_step_handler.increment_theta(n, dtheta, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

    // This method increments the specified value in the omega buffer
    void increment_omega(size_t n,                                                                              // index of the value to increment
                         field_value_t const & domega,                                                          // value of the angular velocity increment
                         typename field_container_t::const_iterator x_begin_itr [[maybe_unused]],               // iterator pointing to the start of the x buffer
                         typename field_container_t::const_iterator v_begin_itr [[maybe_unused]],               // iterator pointing to the start of the v buffer
                         typename field_container_t::const_iterator a_begin_itr [[maybe_unused]],               // iterator pointing to the start of the a buffer
                         typename field_container_t::const_iterator theta_begin_itr [[maybe_unused]],           // iterator pointing to the start of the theta buffer
                         typename field_container_t::iterator omega_begin_itr,                                  // iterator pointing to the start of the omega buffer
                         typename field_container_t::const_iterator alpha_begin_itr [[maybe_unused]]) {   // iterator pointing to the start of the alpha buffer

        base_step_handler.increment_omega(n, domega, x_begin_itr, v_begin_itr, a_begin_itr, theta_begin_itr, omega_begin_itr, alpha_begin_itr);
    }

private:
    const real_t dt;
    afm_filter_ode<real_t> & filter;
    const size_t tip_point, dump_period;
    size_t counter = 0;
    const real_t mass;
    std::vector<bool> fixed_particles;
    base_step_handler_t<field_container_t, field_value_t> base_step_handler;
    typename std::vector<real_t>::iterator tip_force_buffer;
};

#endif //SOOT_INDENTATION_AFM_STEP_HANDLER_H
