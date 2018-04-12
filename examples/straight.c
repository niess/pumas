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
 * integrated between kinetic_min and kinetic_max. Otherwise a point estimate
 * of the flux is done, at the provided kinetic energy.
 */

/* Standard library includes */
#include <math.h>
#include <stdlib.h>
/* The PUMAS API */
#include "pumas.h"

/* The name of the medium's material */
#define MATERIAL_NAME "StandardRock"

/* The medium uniform density, in kg / m^3 */
#define MEDIUM_DENSITY 2.65E+03

#ifndef M_PI
/* Define pi, if unknown */
#define M_PI 3.14159265358979323846
#endif

/* A handle for the PUMAS simulation context */
static struct pumas_context * context = NULL;

/* Gracefully exit to the OS */
static int exit_gracefully(int rc)
{
        pumas_context_destroy(&context);
        pumas_finalise();
        exit(rc);
}

/* Error handler for PUMAS with a graceful exit */
static void error_handler(
    enum pumas_return rc, pumas_function_t * caller, struct pumas_error * error)
{
        /* Dump the error summary */
        fprintf(stderr, "error : ");
        pumas_error_print(stderr, rc, caller, error);
        fprintf(stderr, "\n");

        /* Exit to the OS */
        exit_gracefully(EXIT_FAILURE);
}

/* A uniform rock medium */
static double locals_rock(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals)
{
        /* Set the medium density */
        locals->density = MEDIUM_DENSITY;

        /* Propose a maximum stepping distance. Returning zero or less indicates
         * a uniform medium
         */
        return 0.;
}

/* The medium container */
static struct pumas_medium medium = { 0, &locals_rock };

/* A basic medium callback providing an infinite single medium */
static double medium1(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr)
{
        /* Set the medium */
        *medium_ptr = &medium;

        /* Propose a maximum stepping distance. Returning zero or less indicates
         * an infinite medium
         */
        return 0.;
}

/* A basic Pseudo Random Number Generator (PRNG) providing a uniform
 * distribution over [0, 1]
 */
static double uniform01(struct pumas_context * context)
{
        return rand() / (double)RAND_MAX;
}

/* Gaisser's spectrum, from the PDG */
double spectrum_gaisser(double cos_theta, double kinetic)
{
        const double r_pi = 1.;
        const double E_pi = 115.;
        const double r_K = 0.054;
        const double E_K = 850.;
        const double E = kinetic + 0.10566;
        const double E_star = E * cos_theta;
        return 1.4 * pow(E, -2.7) * (r_pi / (1. + 1.1 * E_star / E_pi) +
                                        r_K / (1. + 1.1 * E_star / E_K));
}

/* The executable main entry point */
int main(int narg, char * argv[])
{
        /* Check the number of arguments */
        if (narg < 4) {
                fprintf(stderr, "Usage: %s ROCK_THICKNESS ELEVATION "
                                "KINETIC_ENERGY[_MIN] [KINETIC_ENERGY_MAX]\n",
                    argv[0]);
                exit_gracefully(EXIT_FAILURE);
        }

        /* Parse the arguments */
        const double rock_thickness = strtod(argv[1], NULL);
        const double elevation = strtod(argv[2], NULL);
        const double kinetic_min = strtod(argv[3], NULL);
        const double kinetic_max =
            (narg >= 5) ? strtod(argv[4], NULL) : kinetic_min;

        /* Set the error handler callback. Whenever an error occurs during a
         * PUMAS function call, the supplied error handler will be evaluated,
         * resulting in an exit to the OS
         */
        pumas_error_handler_set(&error_handler);

        /* Initialise PUMAS from a Material Description File (MDF). This can
         * a few seconds, depending on the number of materials in the MDF.
         */
        pumas_initialise(PUMAS_PARTICLE_MUON, "materials/mdf/standard.xml",
            "../dedx/muon", NULL);

        /* Map the PUMAS material index */
        pumas_material_index(MATERIAL_NAME, &medium.material);

        /* Create a new PUMAS simulation context */
        pumas_context_create(0, &context);

        /* Configure the context for simulating the detailed energy loss, à
         * la Geant4
         */
        context->scheme = PUMAS_SCHEME_DETAILED;

        /* Do a backward transport */
        context->forward = 0;

        /* Disable any transverse transport */
        context->longitudinal = 1;

        /* Set the medium callback */
        context->medium = &medium1;

        /* Provide a PRNG for the Monte-Carlo simulation */
        context->random = &uniform01;

        /* Set a distance limit for the transport as the total rock depth */
        context->distance_max = (rock_thickness <= 0.) ? 1E-06 : rock_thickness;

        /* Run the Monte-Carlo */
        const double cos_theta = cos((90. - elevation) / 180. * M_PI);
        const double sin_theta = sqrt(1. - cos_theta * cos_theta);
        const double rk = log(kinetic_max / kinetic_min);
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
                        kf = kinetic_min * exp(rk * uniform01(context));
                        wf = kf * rk;
                } else {
                        /* A point estimate is computed, for a fixed final
                         * state energy.
                         */
                        kf = kinetic_min;
                        wf = 1;
                }
                struct pumas_state state = {.charge = -1.,
                        .kinetic = kf,
                        .weight = wf,
                        .direction = { -sin_theta, 0., -cos_theta } };

                /* Transport the muon backwards */
                pumas_transport(context, &state);

                /* Update the integrated flux */
                const double wi = state.weight *
                    spectrum_gaisser(-state.direction[2], state.kinetic);
                w += wi;
                w2 += wi * wi;
        }

        /* Print the (integrated) flux */
        w /= n;
        const double sigma =
            (rock_thickness <= 0.) ? 0. : sqrt(((w2 / n) - w * w) / n);

        const char * unit = rk ? "" : "GeV^{-1} ";
        printf("Integrated flux : %.5lE \\pm %.5lE %sm^{-2} s^{-2} sr^{-1}\n",
            w, sigma, unit);

        /* Exit to the OS */
        exit_gracefully(EXIT_SUCCESS);
}
