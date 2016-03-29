/**
 * This is a basic example illustrating the propagation of muons in the
 * Continuously Slowing Down Approximation (CSDA). The CSDA muons are propagated
 * in an infinite and uniform medium. A stopping condition is given by a kinetic
 * energy threshold and/or a total path length limitation.
 *
 * The present file implements and exports a set of routines for test and
 * validation tools. The corresponding interface is declared in the file
 * [csda.h](https://github.com/niess/pumas/blob/master/examples/csda.h).
 * Examples of steering are given in files
 * [csda-propagation.c](csda-propagation.html) and
 * [csda-tracking.c](csda-tracking.html).
 */
/* C89 standard library. */
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
/* The PuMAS API. */
#include "pumas.h"
/* The CSDA interface. */
#include "csda.h"

/**
 * The settings for the propagation of CSDA muons are stored localy in a static
 * structure and exported by a getter/setter function.
 */
/* The propagation settings initialised with defaults. */
static struct csda_settings settings = {1E+01, 1E-01, 0., 0., 5E-05, 90.,
        2.00E+03, "Carbon"};

/* Get/setter for the propagation settings. */
struct csda_settings * csda_settings(void)
{
        return &settings;
}

/**
 * First we need to set a basic uniform medium. In PUMAS the local, i.e.
 * mutable, medium properties are provided by the user to the API with a
 * `pumas_locals_setter_t` callback. In this case, since the local medium
 * properties are constant, we can define them with a static structure that is
 * copied by the medium's setter callback.
 */
/* The uniform medium local properties. */
static struct pumas_locals medium_locals = {0., {0., 0., 0.}};

/* Setter for the uniform medium. */
static double locals_setter(const struct pumas_state * state,
        struct pumas_locals * locals)
{
        memcpy(locals, &medium_locals, sizeof(medium_locals));
        return settings.step_max;
}

/**
 * Secondly, we need to configure a context for the simulation, i.e. a thread
 * safe running environment. We proceed as following:
 *
 * + The context is configured with a single medium and a uniform magnetic
 * field. The latter is computed at the context configuration. Note that it if
 * multiple simulation contexts are used they de facto share the same simulation
 * medium.
 *
 * + The context is configured for pure CSDA muons. That for stochastic energy
 * losses and transverse transportation are disabled.
 *
 * + A kinetic limit and maximum propagation distance are registered.
 */
/* Configure a simulation context for CSDA propagation. */
enum pumas_return csda_configure(struct pumas_context * context)
{
        /* Protect against rounding errors in sine and cosine. */
        #define PROTECT(x) (fabs(x) < FLT_EPSILON ? 0. : x)

        #ifndef M_PI
        /* Define pi, if unknown. */
        #define M_PI 3.14159265358979323846
        #endif

        /* Define a single medium. */
        static struct pumas_medium medium;
        context->media = &medium;

        const int material = pumas_material_index(settings.material);
        if (material < 0) return PUMAS_ERROR;
        medium.material = material;
        medium.setter = &locals_setter;

        /* Precompute the uniform magnetic field. */
        const double deg = M_PI/180.;
        const double angle = settings.magnet_angle*deg;
        medium_locals.density = settings.density;
        medium_locals.magnet[0] = 0.;
        medium_locals.magnet[1] = PROTECT(settings.magnet_amplitude*sin(angle));
        medium_locals.magnet[2] = PROTECT(settings.magnet_amplitude*cos(angle));

        /* Configure the propagation mode. */
        context->forward = 1;
        context->longitudinal = 1;
        context->scheme = PUMAS_SCHEME_CSDA;
        context->kinetic_limit = settings.kinetic_limit;
        context->distance_max = settings.distance_max;

        return PUMAS_SUCCESS;
}

/**
 * Finaly we can run the simulation. That for we call the `pumas_propagate` API
 * function given a simulation context and an initial state. At return the
 * provided `pumas_state` is updated with the final state properties.
 */
/* Propagate a PuMAS state with CSDA. */
enum pumas_return csda_propagate(struct pumas_context * context,
        struct pumas_state * state)
{
        /* Initialise the state. */
        struct pumas_state initial = {-1, settings.kinetic_energy,
                0., 0., 0., 1., {0., 0., 0.}, {0., 0., 1.}};
        memcpy(state, &initial, sizeof(*state));

        /* Propagate. */
        return pumas_propagate(context, state);
}

/**
 * The following routine is used by the test tools. It dumps a summary of
 * the muon state to a stream.
 */
/* Dump a PuMAS state to a stream. */
void csda_dump(FILE * stream, struct pumas_state * state)
{
        if (stream == NULL) return;

        /* Dump the state. */
        const double * p = state->position;
        const double * u = state->direction;
        fprintf(stream, "%10.3E %10.3E %10.3E  %10.3E %10.3E  "
                "%10.3E %10.3E\n", state->kinetic, state->distance,
                state->time, p[0], p[2], u[0], u[2]);
}
