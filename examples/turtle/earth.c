/*
 * This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * For more information, please refer to <http://unlicense.org>
 */

/* This a basic example illustrating the usage of TURTLE for representing an
 * Earth geometry. Note that you'll need TURTLE along with PUMAS for this
 * example to compile.
 *
 * Note that for this example to work a PUMAS physics dump must have been
 * generated first e.g. by running the `example-dump` program (under
 * examples/pumas/dump.c).
 */

/* Standard library includes */
#include <math.h>
#include <stdlib.h>
/* The PUMAS API */
#include "pumas.h"
#include "turtle.h"

/* Handles for PUMAS Physics & simulation context */
static struct pumas_physics * physics = NULL;
static struct pumas_context * context = NULL;

/* Handle for TURTLE Monte Carlo stepper (ray tracer) */
static struct turtle_stepper * stepper = NULL;

/* The media container. The locals callback are set to `NULL` resulting in
 * the default material densities being used with a null magnetic field
 */
#define NUMBER_OF_LAYERS 3
static struct pumas_medium earth_media[NUMBER_OF_LAYERS] = {
    {-1, NULL}, {-1, NULL}, {-1, NULL}};

/* A medium callback using a `turtle_stepper` in order to represent an
 * Earth geometry
 */
static enum pumas_step earth_medium(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr,
    double * step_ptr)
{
        /* Do a `turtle_stepper` step. The stepper returns a tentative step
         * length as well as the index of the topography layer at the end step
         * position. Note that the particle is not moved however.
         *
         * The `turtle_stepper` can provide extra informations that is not used
         * in this example, e.g. the particle latitude, longitude and altitude,
         * or the topography layer bottom and top elevations.
         */
        double step;
        int index[2]; /* Note that while index[1] is not used in this example
                       * the index array must at least be of size 2.
                       */
        turtle_stepper_step(stepper, state->position, NULL, NULL, NULL,
            NULL, NULL, &step, index);

        if (step_ptr != NULL) *step_ptr = step;

        if (medium_ptr != NULL) {
                if ((index[0] >= 0) && (index[0] < NUMBER_OF_LAYERS)) {
                        *medium_ptr = earth_media + index[0];
                } else {
                        *medium_ptr = NULL;
                }
        }

        return PUMAS_STEP_CHECK;
}

/* Print the given Monte Carlo state to stdout */
static void print_state(struct pumas_state * state)
{
        static int step = 0;

        /* Get the current medium and material */
        struct pumas_medium * medium;
        earth_medium(context, state, &medium, NULL);
        const char * material;
        if (medium == NULL) {
                material = "(void)";
        } else {
                pumas_physics_material_name(
                    physics, medium->material, &material);
        }

        /* Get the altitude w.r.t. the ellipsoid (GPS altitude) */
        double altitude;
        turtle_ecef_to_geodetic(state->position, NULL, NULL, &altitude);
        printf("%2d. energy = %.3E, altitude = %8.2f, material = %s\n",
            step++, state->energy, altitude, material);
}

/* The executable main entry point */
int main(int narg, char * argv[])
{
        /* Check the number of arguments */
        if (narg < 4) {
                fprintf(stderr,
                    "Usage: %s AZIMUTH ELEVATION KINETIC_ENERGY\n",
                    argv[0]);
                exit(EXIT_FAILURE);
        }

        /* Parse the arguments */
        const double azimuth = strtod(argv[1], NULL);
        const double elevation = strtod(argv[2], NULL);
        const double energy = strtod(argv[3], NULL);

        /* Initialise PUMAS physics from a binary dump, e.g. generated by the
         * `pumas/loader` example
         */
        const char * dump_file = "examples/data/materials.pumas";
        FILE * fid = fopen(dump_file, "rb");
        if (fid == NULL) {
                perror(dump_file);
                fputs(
                    "Please run pumas-loader (example/pumas/loader.c) to "
                    "generate the physics dump.\n", stderr);
                exit(EXIT_FAILURE);
        }
        pumas_physics_load(&physics, fid);
        fclose(fid);

        /* Map the PUMAS materials indices */
        pumas_physics_material_index(
            physics, "StandardRock", &earth_media[0].material);
        pumas_physics_material_index(
            physics, "Water", &earth_media[1].material);
        pumas_physics_material_index(physics, "Air", &earth_media[2].material);

        /* Create a new PUMAS simulation context */
        pumas_context_create(&context, physics, 0);

        /* Set the medium callback */
        context->medium = &earth_medium;

        /* Configure the transport for stopping at each change of medium */
        context->event |= PUMAS_EVENT_MEDIUM;

        /* Create and configure a TURTlE stepper
         *
         * The Earth geometry is made of three layers with a flat topography.
         * These layers are mapped to the media array in the `earth_medium`
         * callback. The bottom layer is mapped to `StandardRock`, the middle
         * one to `Water` and the the top one to air. The layers are spaced by
         * 1 km each. Hence, this represents a fictious Earth covered with a
         * 1 km deep ocean and a 1 km high atmosphere all of uniform density.
         *
         * Note that more complex topographies can be used at this step rather
         * than flat ones. This is done by adding `turtle_map` and/or
         * `turtle_stack` objects to the topography layer. Please refer to the
         * TURTLE examples and doc for more details.
         */
        turtle_stepper_create(&stepper);
        turtle_stepper_add_flat(stepper, -1E+03);
        turtle_stepper_add_layer(stepper);
        turtle_stepper_add_flat(stepper, 0);
        turtle_stepper_add_layer(stepper);
        turtle_stepper_add_flat(stepper, 1E+03);

        /* Initialise the Monte-Carlo state. Geodetic coordinates (latitude,
         * longitude) are transformed to Earth-Centered Earth-Fixed (ECEF) ones
         * using TURTLE. The particle is located 0.5 m below the top of the
         * lowest layer, i.e. rocks.
         *
         * Note that the state is initialised using C99 designated initializers.
         * Therefore, unspecified fields are set to zero, e.g. the travelled
         * *distance* and the *decayed* flag.
         *
         * Please take care to properly initialize all fields when creating a
         * new state. In particular the *decay* flag must be zero (false) and
         * the Monte Carlo weight strictly positive.
         */
        struct pumas_state state = { .charge = -1.,
            .energy = energy,
            .weight = 1
        };

        const double latitude = 45, longitude = 3;
        turtle_stepper_position(stepper, latitude, longitude, -0.5,
            0, state.position, NULL);
        turtle_ecef_from_horizontal(
            latitude, longitude, azimuth, elevation, state.direction);

        /* Transport the particle */
        print_state(&state);
        for (;;) {
                enum pumas_event event;
                struct pumas_medium * media[2];
                 pumas_context_transport(context, &state, &event, media);

                print_state(&state);

                if ((state.energy == 0) || (media[1] == NULL)) break;
        }

        /* Clean and exit to the OS */
        turtle_stepper_destroy(&stepper);
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);

        exit(EXIT_SUCCESS);
}
