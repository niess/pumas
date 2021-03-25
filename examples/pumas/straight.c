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

/* This a basic example illustrating the backward computation of a transmitted
 * muon flux through a constant thickness of a uniform material, e.g.
 * Standard Rock. If a maximum kinetic energy is provided the flux is
 * integrated between energy_min and energy_max. Otherwise a point estimate
 * of the flux is done, at the provided kinetic energy.
 *
 * Note that for this example to work you need the corresponding MDF and energy
 * loss tables. Those can be downloaded with git, as following:
 *
 * git clone https://gitub.com/niess/pumas-materials materials
 */

/* Standard library includes */
#include <math.h>
#include <stdlib.h>
/* The PUMAS API */
#include "pumas.h"
/* The atmospheric muon fluxes library */
#include "flux.h"

/* The name of the medium's material */
#define MATERIAL_NAME "StandardRock"

#ifndef M_PI
/* Define pi, if unknown */
#define M_PI 3.14159265358979323846
#endif

/* Handles for PUMAS Physics & simulation context */
static struct pumas_physics * physics = NULL;
static struct pumas_context * context = NULL;

/* The medium container. The locals callback is set to `NULL` resulting in
 * the default material density being used with a null magnetic field
 */
static struct pumas_medium medium = { 0, NULL };

/* A basic medium callback providing an infinite single medium */
static enum pumas_step medium1(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr,
    double * step_ptr)
{
        /* Set the medium */
        if (medium_ptr != NULL) *medium_ptr = &medium;

        /* Propose a maximum stepping distance. Returning zero or less indicates
         * an infinite medium
         */
        if (step_ptr != NULL) *step_ptr = 0.;

        return PUMAS_STEP_CHECK;
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
        const double rock_thickness = strtod(argv[1], NULL);
        const double elevation = strtod(argv[2], NULL);
        const double energy_min = strtod(argv[3], NULL);
        const double energy_max =
            (narg >= 5) ? strtod(argv[4], NULL) : energy_min;

        /* Initialise PUMAS physics from a Material Description File (MDF). This
         * can a few seconds, depending on the number of materials in the MDF.
         */
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/mdf/examples/standard.xml", "materials/dedx/muon", NULL);

        /* Map the PUMAS material index */
        pumas_physics_material_index(physics, MATERIAL_NAME, &medium.material);

        /* Create a new PUMAS simulation context */
        pumas_context_create(&context, physics, 0);

        /* Configure the context for simulating the detailed energy loss, à
         * la Geant4
         */
        context->mode.energy_loss = PUMAS_MODE_DETAILED;

        /* Do a backward transport */
        context->mode.direction = PUMAS_MODE_BACKWARD;

        /* Disable any transverse transport */
        context->mode.scattering = PUMAS_MODE_LONGITUDINAL;

        /* Set the medium callback */
        context->medium = &medium1;

        /* Set a distance limit for the transport as the total rock depth */
        context->limit.distance =
            (rock_thickness <= 0.) ? 1E-06 : rock_thickness;
        context->event |= PUMAS_EVENT_LIMIT_DISTANCE;

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
                         * a log-uniform distribution. The Monte-Carlo weight is
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
                struct pumas_state state = { .charge = cf,
                        .energy = kf,
                        .weight = wf,
                        .direction = { -sin_theta, 0., -cos_theta } };

                /* Transport the muon backwards */
                pumas_context_transport(context, &state, NULL, NULL);

                /* Update the integrated flux */
                const double wi = state.weight *
                    flux_gccly(-state.direction[2], state.energy, state.charge);
                w += wi;
                w2 += wi * wi;
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
