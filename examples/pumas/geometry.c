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

/* This example illustrates the backward computation of a transmitted muon flux
 * through a simple geometry composed of two layers: Standard Rock and Air. The
 * Air medium has an exponential density profile. If a maximum kinetic energy is
 * provided the flux is integrated between energy_min and energy_max.  Otherwise
 * a point estimate of the flux is done, at the provided kinetic energy.
 *
 * Note that for this example to work a PUMAS physics dump must have been
 * generated first e.g. by running the `example-dump` program (under
 * examples/pumas/dump.c).
 */

/* Standard library includes */
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
/* The PUMAS API */
#include "pumas.h"
/* The atmospheric muon fluxes library */
#include "flux.h"

#ifndef M_PI
/* Define pi, if unknown */
#define M_PI 3.14159265358979323846
#endif

/* Altitude at which the primary flux is sampled */
#define PRIMARY_ALTITUDE 1E+03

/* Handles for PUMAS Physics and simulation context */
static struct pumas_physics * physics = NULL;
static struct pumas_context * context = NULL;

/* The uniform rock medium. Note that the geomagnetic field can be neglected
 * in rocks
 */
static double locals_rock(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals)
{
        /* Set the medium density, in kg / m3. Setting zero or less results in
         * the default material density being used */
        locals->density = 2.65E+03;

        /* Propose a maximum stepping distance. Returning zero or less indicates
         * a uniform medium
         */
        return 0.;
}

/* The non uniform atmosphere, using an exponential model */
static double locals_air(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals)
{
        /* Set the geomagnetic field, assumed uniform */
        locals->magnet[0] = 0.;
        locals->magnet[1] = 2E-05;
        locals->magnet[2] = -4E-05;

        /* Set the atmosphere density, depending on the altitude a.s.l. */
        const double rho0 = 1.205;
        const double h = 12E+03;
        locals->density = rho0 * exp(-state->position[2] / h);

        /* Propose a maximum stepping distance as 1 percent of the projected
         * attenuation length, for the density
         */
        const double eps = 5E-02;
        const double uz = fabs(state->direction[2]);
        return 1E-02 * h / ((uz <= eps) ? eps : uz);
}

/* The media container */
static struct pumas_medium media[2] = { { 0, &locals_rock },
        { 1, &locals_air } };

/* A simple medium with a flat rock layer and a flat atmosphere */
static double rock_thickness = 0.;

static enum pumas_step medium2(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr,
    double * step_ptr)
{
        if ((medium_ptr == NULL) && (step_ptr == NULL))
                return PUMAS_STEP_RAW;

        /* Check the muon position and direction */
        const double z = state->position[2];
        const double sgn = (context->mode.direction == PUMAS_MODE_FORWARD) ?
            1 : -1;
        const double uz = state->direction[2] * sgn; /* Note that in backward
                                                      * mode the muon propagates
                                                      * reverse to the state
                                                      * direction
                                                      */
        double step = 1E+03;
        if (z < 0.) {
                /* The muon is outside of the simulation area */
                if (medium_ptr != NULL) *medium_ptr = NULL;
                step = -1.;
        } else if (z < rock_thickness) {
                if (medium_ptr != NULL) *medium_ptr = media;
                if (uz > FLT_EPSILON)
                        /* The muon is (backward) upgoing. The next boundary is
                         * the rock-air interface. */
                        step = (rock_thickness - z) / uz;
                else if (uz < -FLT_EPSILON)
                        /* The muon is (backward) downgoing. The next boundary
                         * is the rock bottom. */
                        step = -z / uz;
        } else if (z < PRIMARY_ALTITUDE) {
                if (medium_ptr != NULL) *medium_ptr = media + 1;
                if (uz > FLT_EPSILON)
                        /* The muon is (backward) upgoing. The next boundary is
                         * the air top. */
                        step = (PRIMARY_ALTITUDE - z) / uz;
                else if (uz < -FLT_EPSILON)
                        /* The muon is (backward) downgoing. The next boundary
                         * is the rock-air interface. */
                        step = (rock_thickness - z) / uz;
        } else {
                /* The muon is outside of the simulation area */
                if (medium_ptr != NULL) *medium_ptr = NULL;
                step = -1.;
        }

        if (step_ptr != NULL) {
#define STEP_EPSILON 1E-07
                /* Offset slightly the step length in order to end in the new
                 * medium
                 */
                if (step > 0) step += STEP_EPSILON;
                *step_ptr = step;
        }

        return PUMAS_STEP_RAW;
}

/* The executable main entry point */
int main(int narg, char * argv[])
{
        /* Check the number of arguments */
        if (narg < 4) {
                fprintf(stderr,
                    "Usage: %s ROCK_THICKNESS ELEVATION "
                    "KINETIC_ENERGY[_MIN] [KINETIC_ENERGY_MAX]\n",
                    argv[0]);
                exit(EXIT_FAILURE);
        }

        /* Parse the arguments */
        rock_thickness = strtod(argv[1], NULL);
        if ((rock_thickness < 0.) || (rock_thickness > PRIMARY_ALTITUDE)) {
                errno = EINVAL;
                perror("rock thickness");
                exit(EXIT_FAILURE);
        }
        const double elevation = strtod(argv[2], NULL);
        const double energy_min = strtod(argv[3], NULL);
        const double energy_max =
            (narg >= 5) ? strtod(argv[4], NULL) : energy_min;

        /* Initialise PUMAS physics from a binary dump, e.g. generated by the
         * `dump` example
         */
        const char * dump_file = "examples/data/materials.pumas";
        FILE * fid = fopen(dump_file, "rb");
        if (fid == NULL) {
                perror(dump_file);
                exit(EXIT_FAILURE);
        }
        pumas_physics_load(&physics, fid);
        fclose(fid);

        /* Map the PUMAS material indices */
        pumas_physics_material_index(
            physics, "StandardRock", &media[0].material);
        pumas_physics_material_index(
            physics, "Air", &media[1].material);

        /* Create a new PUMAS simulation context */
        pumas_context_create(&context, physics, 0);

        /* Configure the context for a backward transport */
        context->mode.direction = PUMAS_MODE_BACKWARD;

        /* Set the medium callback */
        context->medium = &medium2;

        /* Enable external limit on the kinetic energy */
        context->event |= PUMAS_EVENT_LIMIT_ENERGY;

        /* Run the Monte-Carlo */
        const double cos_theta = cos((90. - elevation) / 180. * M_PI);
        const double sin_theta = sqrt(1. - cos_theta * cos_theta);
        const double rk = log(energy_max / energy_min);
        double w = 0., w2 = 0.;
        const int n = 10000;
        int i;
        for (i = 0; i < n; i++) {
                /* Set the muon final state */
                double kf, wf;
                if (rk) {
                        /* The final state kinetic energy is randomised over
                         * a log-uniform distribution. The Monte Carlo weight is
                         * initialised according to this generating bias PDF,
                         * i.e. wf = 1 / PDF(kf).
                         */
                        kf = energy_min * exp(rk * context->random(context));
                        wf = kf * rk;
                } else {
                        /* A point estimate is computed, for a fixed final
                         * state energy.
                         */
                        kf = energy_min;
                        wf = 1;
                }
                const double cf = (context->random(context) > 0.5) ? 1 : -1;
                wf *= 2; /* Update the Monte Carlo weight according to the bias
                          * PDF used for the charge randomisation, i.e. 1 / 0.5
                          */

                /* Create a new Monte Carlo state.
                 *
                 * Note that the state is initialised using C99 designated
                 * initializers. Therefore, unspecified fields are set to zero,
                 * e.g. the travelled *distance* and the *decayed* flag.
                 *
                 * Please take care to properly initialize all fields when
                 * creating a new state. In particular the *decay* flag must be
                 * zero (false) and the Monte Carlo weight strictly positive.
                 */
                struct pumas_state state = { .charge = cf,
                        .energy = kf,
                        .weight = wf,
                        .direction = { -sin_theta, 0., -cos_theta } };

                /* Transport the muon backwards */
                const double energy_threshold = energy_max * 1E+03;
                while (state.energy < energy_threshold - FLT_EPSILON) {
                        if (state.energy < 1E+02 - FLT_EPSILON) {
                                /* Below 100 GeV do a detailed simulation
                                 * à la Geant4, including transverse transport
                                 */
                                context->mode.energy_loss =
                                    PUMAS_MODE_DETAILED;
                                context->mode.scattering =
                                    PUMAS_MODE_FULL_SPACE;
                                context->limit.energy = 1E+02;
                        } else {
                                /* Do a fast simulation à la MUM */
                                context->mode.energy_loss =
                                    PUMAS_MODE_HYBRID;
                                context->mode.scattering =
                                    PUMAS_MODE_LONGITUDINAL;
                                context->limit.energy = energy_threshold;
                        }
                        enum pumas_event event;
                        struct pumas_medium * medium[2];
                        pumas_context_transport(
                            context, &state, &event, medium);

                        /* Check if the muon has exit the simulation area */
                        if (event == PUMAS_EVENT_MEDIUM) {
                                if (medium[1] == NULL) {
                                        /* Update the integrated flux */
                                        const double wi = state.weight *
                                        flux_gccly(-state.direction[2],
                                            state.energy, state.charge);
                                        w += wi;
                                        w2 += wi * wi;
                                        break;
                                }
                        } else if (event != PUMAS_EVENT_LIMIT_ENERGY) {
                                /* This should not happen */
                                fprintf(stderr,
                                    "error: unexpected PUMAS event `%d`\n",
                                    event);
                                exit(EXIT_FAILURE);
                        }
                }
        }

        /* Print the (integrated) flux */
        w /= n;
        const double sigma =
            (rock_thickness <= 0.) ? 0. : sqrt(((w2 / n) - w * w) / n);
        const char * unit = rk ? "" : "GeV^{-1} ";
        printf("Flux : %.5lE \\pm %.5lE %sm^{-2} s^{-2} sr^{-1}\n", w, sigma,
            unit);

        /* Clean and exit to the OS */
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);

        exit(EXIT_SUCCESS);
}
