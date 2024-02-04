//
// Created by egor on 2/3/24.
//

#include <vector>
#include <iostream>

#include <Eigen/Eigen>

#include <libgran/contact_force/contact_force.h>
#include <libgran/hamaker_force/hamaker_force.h>
#include <libgran/sinter_bridge/alt_sinter_bridge.h>
#include <libgran/granular_system/granular_system.h>

#include "../reader/mackowski_reader.h"
#include "../writer/writer.h"
#include "../energy/compute_energy.h"

int main() {

    /* SIMULATION PARAMETERS WILL BE HERE */

    const double r_part = 14.0e-9;

    // Load the aggregate from file
    auto x0 = load_mackowski_aggregate("../aggregates/size_50/aggregate_1.txt");
    // Pre-process the loaded aggregate by removing overlap
    remove_mackowski_overlap(x0, 100, 0.01, 1.0);
    // Re-size the monomers
    for (auto & point : x0) {
        point *= r_part;
    }

    return 0;
}
