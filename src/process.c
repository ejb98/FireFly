#include "wing.h"
#include "solve.h"
#include "get_size.h"
#include "shed_wake.h"
#include "rollup_wake.h"
#include "compute_velocities.h"

void process(Wing *wing, double delta_time) {
    wing->iteration++;

    if (wing->iteration) {
        compute_velocities(wing, delta_time);
    }

    shed_wake(wing);
    solve(wing);

    if (wing->iteration) {
        rollup_wake(wing, delta_time);
    }

    wing->position_prev = wing->position;
    wing->rotation_prev = wing->rotation;
}