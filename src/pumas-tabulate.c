/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C library whose purpose is to transport high energy
 * muons or taus in various media.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "pumas-tabulate.h"
#define OPTPARSE_IMPLEMENTATION
#define OPTPARSE_API static
#include "optparse.h"

/* Clean and exit. */
void exit_gracefully(int code)
{
        tabulation_finalise();
        pumas_finalise();
        exit(code);
}

/* Show help and exit. */
static void exit_with_help(int code)
{
        /* This is the help string formated as it should output on the screen.
         * Therefore it follows a specific formating rule.
         */
        // clang-format off
	fprintf(stderr,
"Usage: pumas-tabulate [PARTICLE] [MDF] [OPTION]... \n"
"Tabulate the energy loss of the given PARTICLE (`muon` or `tau`) in one or\n"
"more material defined in MDF.\n"
"\n"
"Material options:\n"
"  -d, --density=DEN       set the material density to DEN\n"
"  -I, --mee=MEE           set the mean excitation energy to MEE\n"
"  -m, --material=NAME     switch to the material with the provided NAME\n"
"      --state=STATE       set the material state to STATE\n"
"      --sternheimer=STR   set the sternheimer coefficients for the density\n"
"                            effect according to the values read from STR"
"\n"
"Global options:\n"
"  -f, --force             overwrite any existing output file(s)\n"
"  -h, --help              show this help and exit\n"
"      --kinetic-max=K     set the maximum kinetic energy to K [1E-03]\n"
"      --kinetic-min=K     set the minimum kinetic energy to K [1E+06]\n"
"  -k, --kinetic-n=N       switch to a log sampling of the kinetic energy and\n"
"                            set the number of kinetic energy values to N\n"
"  -o, --outdir            set the output directory for the tabulation(s)\n"
"\n"
"Material's properties must be specified as `-m [NAME] [OPTION]` where the\n"
"subsequent material [OPTION]s refer to the prepending material. Several\n"
"materials can be specified. Note that at least one material must be provi-\n"
"ded in order to proceed.\n"
"\n"
"Densities (DEN) must be given in g/cm^3, kinetic energies (K) in GeV and\n"
"Mean Excitation Energies (MEE) in eV. The material [STATE] must be one of\n"
"`gaz`, `liquid` or `solid`. If not specified, materials with densities below\n"
"0.1 g/cm^3 are considered as `gaz` and above as `liquid`. The Sternheimer\n"
"coefficients must be given as a list of values as \"[A] [K] [X_0] [X_1]\n"
"[CBAR] [DELTA0]\".\n"
"\n"
"The default behaviour is to use the PDG sampling for the kinetic energy\n"
"values. Use the `-k` option in order to switch to a log sampling. Note that\n"
"for the log sampling the minimum kinetic energy must be strictly positive.\n"
	);
	exit_gracefully(code);
        // clang-format on
}

/* PDG sampling of kinetic energy values. */
static double pdg_kinetic[145] = { 1.E-03, 1.2E-03, 1.400E-03, 1.700E-03,
        2.000E-03, 2.500E-03, 3.000E-03, 3.500E-03, 4.000E-03, 4.500E-03,
        5.000E-03, 5.500E-03, 6.000E-03, 7.000E-03, 8.000E-03, 9.000E-03,
        1.E-02, 1.2E-02, 1.400E-02, 1.700E-02, 2.000E-02, 2.500E-02, 3.000E-02,
        3.500E-02, 4.000E-02, 4.500E-02, 5.000E-02, 5.500E-02, 6.000E-02,
        7.000E-02, 8.000E-02, 9.000E-02, 1.E-01, 1.2E-01, 1.400E-01, 1.700E-01,
        2.000E-01, 2.500E-01, 3.000E-01, 3.500E-01, 4.000E-01, 4.500E-01,
        5.000E-01, 5.500E-01, 6.000E-01, 7.000E-01, 8.000E-01, 9.000E-01,
        1.E+00, 1.2E+00, 1.400E+00, 1.700E+00, 2.000E+00, 2.500E+00, 3.000E+00,
        3.500E+00, 4.000E+00, 4.500E+00, 5.000E+00, 5.500E+00, 6.000E+00,
        7.000E+00, 8.000E+00, 9.000E+00, 1.E+01, 1.2E+01, 1.400E+01, 1.700E+01,
        2.000E+01, 2.500E+01, 3.000E+01, 3.500E+01, 4.000E+01, 4.500E+01,
        5.000E+01, 5.500E+01, 6.000E+01, 7.000E+01, 8.000E+01, 9.000E+01,
        1.E+02, 1.2E+02, 1.400E+02, 1.700E+02, 2.000E+02, 2.500E+02, 3.000E+02,
        3.500E+02, 4.000E+02, 4.500E+02, 5.000E+02, 5.500E+02, 6.000E+02,
        7.000E+02, 8.000E+02, 9.000E+02, 1.E+03, 1.2E+03, 1.400E+03, 1.700E+03,
        2.000E+03, 2.500E+03, 3.000E+03, 3.500E+03, 4.000E+03, 4.500E+03,
        5.000E+03, 5.500E+03, 6.000E+03, 7.000E+03, 8.000E+03, 9.000E+03,
        1.E+04, 1.2E+04, 1.400E+04, 1.700E+04, 2.000E+04, 2.500E+04, 3.000E+04,
        3.500E+04, 4.000E+04, 4.500E+04, 5.000E+04, 5.500E+04, 6.000E+04,
        7.000E+04, 8.000E+04, 9.000E+04, 1.E+05, 1.2E+05, 1.400E+05, 1.700E+05,
        2.000E+05, 2.500E+05, 3.000E+05, 3.500E+05, 4.000E+05, 4.500E+05,
        5.000E+05, 5.500E+05, 6.000E+05, 7.000E+05, 8.000E+05, 9.000E+05,
        1.E+06 };

/* Handle for temporary tabulations data. */
static struct tabulation_data data = { sizeof(pdg_kinetic) /
            sizeof(*pdg_kinetic),
        pdg_kinetic, NULL, NULL, 0, NULL, NULL };

/* Handler for PUMAS errors. */
static void handle_error(
    enum pumas_return rc, pumas_function_t * caller, struct pumas_error * error)
{
        /* Dump the error summary. */
        fprintf(stderr, "pumas: library error. See details below.\n");
        fprintf(stderr, "error : ");
        pumas_error_print(stderr, rc, caller, error);
        fprintf(stderr, "\n");

        /* Exit to the OS. */
        exit_gracefully(EXIT_FAILURE);
}

/* Clear the temporary data for the tabulations. */
void tabulation_finalise()
{
        free(data.path);
        if (data.kinetic != pdg_kinetic) free(data.kinetic);

        while (data.e_stack != NULL) {
                struct tabulation_element * prev = data.e_stack->prev;
                free(data.e_stack);
                data.e_stack = prev;
        }

        while (data.m_stack != NULL) {
                struct tabulation_material * prev = data.m_stack->prev;
                free(data.m_stack);
                data.m_stack = prev;
        }
}

/*
 * Create a new material entry and add it on top of the stack.
 */
struct tabulation_material * tabulation_material_create(
    struct tabulation_data * data, int material)
{
        /* Allocate memory for the new material. */
        struct tabulation_material * m = malloc(sizeof(*m));
        if (m == NULL) return NULL;
        memset(m, 0x0, sizeof(*m));
        m->index = material;
        m->state = TABULATION_STATE_UNKNOWN;

        /* Add the material on top of the stack. */
        if (data->m_stack != NULL) data->m_stack->next = m;
        m->prev = data->m_stack;
        m->next = NULL;
        data->m_stack = m;

        return m;
}

/* Drop the material on top of the stack. */
void tabulation_material_drop(struct tabulation_data * data)
{
        if (data->m_stack == NULL) return;
        struct tabulation_material * prev = data->m_stack->prev;
        free(data->m_stack);
        data->m_stack = prev;
        if (prev != NULL) prev->next = NULL;
}

/*
 * Get the requested material from the temporary data and put it on top of
 * the stack.
 */
struct tabulation_material * tabulation_material_get(
    struct tabulation_data * data, int material)
{
        struct tabulation_material * m;
        for (m = data->m_stack; m != NULL; m = m->next) {
                if (m->index == material) {
                        struct tabulation_material * next = m->next;
                        if (next != NULL) {
                                /* Put the material on top of the stack. */
                                struct tabulation_material * prev = m->prev;
                                if (prev != NULL) prev->next = next;
                                next->prev = prev;
                                data->m_stack->next = m;
                                m->prev = data->m_stack;
                                data->m_stack = m;
                        }
                        return m;
                }
        }

        /* The material wasn't found, return `NULL`. */
        return NULL;
}

int main(int argc, char * argv[])
{
        /*
         * Initialise the library in dry mode, i.e. without loading and
         * processing the energy loss tables.
         */
        if (argc <= 2) exit_with_help(EXIT_FAILURE);
        enum pumas_particle particle = -1;
        if (strcmp(argv[1], "muon") == 0)
                particle = PUMAS_PARTICLE_MUON;
        else if (strcmp(argv[1], "tau") == 0)
                particle = PUMAS_PARTICLE_TAU;
        else {
                fprintf(stderr, "pumas-tabulate: error. "
                                "Unknown particle `%s`. "
                                "Call with -h, --help for usage.\n",
                    argv[1]);
                exit_gracefully(EXIT_FAILURE);
        }
        pumas_error_handler_set(handle_error);
        _pumas_initialise_dry(particle, argv[2], NULL);

        /* Parse the input parameters. */
        double kinetic_min = 1E-03, kinetic_max = 1E+06;
        int kinetic_n = 0;
        struct tabulation_material * material = NULL;
        const int n_base = pumas_material_length() - pumas_composite_length();

        struct optparse options;
        optparse_init(&options, argv + 2);
        for (;;) {
                struct optparse_long long_options[] = { /* Material options. */
                        { "density", 'd', OPTPARSE_REQUIRED },
                        { "mee", 'I', OPTPARSE_REQUIRED },
                        { "material", 'm', OPTPARSE_REQUIRED },
                        { "state", 's', OPTPARSE_REQUIRED },
                        { "sternheimer", 'S', OPTPARSE_REQUIRED },

                        /* Global options. */
                        { "force", 'f', OPTPARSE_NONE },
                        { "kinetic-max", 'x', OPTPARSE_REQUIRED },
                        { "kinetic-min", 'i', OPTPARSE_REQUIRED },
                        { "kinetic-n", 'k', OPTPARSE_REQUIRED },
                        { "help", 'h', OPTPARSE_REQUIRED },
                        { "outdir", 'o', OPTPARSE_REQUIRED },

                        { 0 }
                };

                /* Process the next argument. */
                char * endptr = NULL; /* For parsing with str* functions. */
                int option_index;
                if ((option_index = optparse_long(&options, long_options, NULL))
                    == -1) break; /* No more options to parse. */

                /* Material options. */
                else if (option_index == 'd') {
                        if (material == NULL) goto error_no_material;
                        material->density = strtod(options.optarg, &endptr);
                        material->density *= 1E+03;
                } else if (option_index == 'I') {
                        if (material == NULL) goto error_no_material;
                        material->I = strtod(options.optarg, &endptr);
                        material->I *= 1E-09;
                } else if (option_index == 'm') {
                        /* Fetch the material data. */
                        int index;
                        pumas_material_index(options.optarg, &index);
                        if (index >= n_base) {
                                fprintf(stderr,
                                    "pumas-tabulate: error. "
                                    "Only base materials can be tabulated. "
                                    "Call with -h, --help for usage.\n");
                                exit_gracefully(EXIT_FAILURE);
                        }

                        material = tabulation_material_get(&data, index);
                        if (material == NULL) {
                                /* Create the new data entry. */
                                material =
                                    tabulation_material_create(&data, index);
                                if (material == NULL) {
                                        perror("pumas-tabulate");
                                        exit_gracefully(EXIT_FAILURE);
                                }
                        }
                } else if (option_index == 's') {
                        if (material == NULL) goto error_no_material;
                        if (strcmp(options.optarg, "gaz") == 0)
                                material->state = TABULATION_STATE_GAZ;
                        else if (strcmp(options.optarg, "liquid") == 0)
                                material->state = TABULATION_STATE_LIQUID;
                        else if (strcmp(options.optarg, "solid") == 0)
                                material->state = TABULATION_STATE_SOLID;
                        else {
                                fprintf(stderr,
                                    "pumas-tabulate: error. "
                                    "Unknown state `%s`. "
                                    "Call with -h, --help for usage.\n",
                                    options.optarg);
                                exit_gracefully(EXIT_FAILURE);
                        }
                } else if (option_index == 'S') {
                        if (material == NULL) goto error_no_material;
                        if (sscanf(options.optarg, "%lf%lf%lf%lf%lf%lf",
                                &material->a, &material->k, &material->x_0,
                                &material->x_1, &material->Cbar,
                                &material->delta0) != 6) {
                                fprintf(stderr,
                                    "pumas-tabulate: error. "
                                    "Invalid Sternheimer coefficients. "
                                    "Call with -h, --help for usage.\n");
                                exit_gracefully(EXIT_FAILURE);
                        }
                        continue;
                }

                /* Global options. */
                else if (option_index == 'f') {
                        data.overwrite = 1;
                        continue;
                } else if ((option_index == 'h') || (option_index == '?')) {
                        exit_with_help(EXIT_SUCCESS);
                } else if (option_index == 'x') {
                        kinetic_max = strtod(options.optarg, &endptr);
                } else if (option_index == 'i') {
                        kinetic_min = strtod(options.optarg, &endptr);
                } else if (option_index == 'k') {
                        kinetic_n = (int)strtol(options.optarg, &endptr, 0);
                } else if (option_index == 'o') {
                        data.outdir = options.optarg;
                        continue;
                } else {
                        /*
                         * getopt should already have reported
                         * an error.
                         */
                        exit_gracefully(EXIT_FAILURE);
                }

                /* Check the parsing. */
                if (endptr == options.optarg) errno = EINVAL;
                if (errno != 0) {
                        perror("pumas-tabulate");
                        exit_gracefully(EXIT_FAILURE);
                }
                continue;

        error_no_material:
                /* No material was defined. */
                fprintf(stderr, "pumas-tabulate: error. No material specified "
                                "Call with -h, --help for usage.\n");
                exit_gracefully(EXIT_FAILURE);
        }

        /* Check the consistency of the inputs. */
        if (material == NULL) goto error_no_material;

        for (material = data.m_stack; material != NULL;
             material = material->prev) {
                if (material->a <= 0.) {
                        /* No explicit Sternheimer coefficients. We need at
                         * at least the density to proceed.
                         */
                        if (material->density <= 0.) {
                                fprintf(stderr,
                                    "pumas-tabulate: error."
                                    " No density and Sternheimer "
                                    "coefficient specified "
                                    "Call with -h, --help for usage.\n");
                                exit_gracefully(EXIT_FAILURE);
                        }

                        if (material->state == TABULATION_STATE_UNKNOWN) {
                                /* Guess the state from the density. */
                                if (material->density < 0.1E+03)
                                        material->state = TABULATION_STATE_GAZ;
                                else
                                        material->state =
                                            TABULATION_STATE_LIQUID;
                        }
                }
        }

        /* Set the kinetic energies. */
        if (kinetic_n > 0) {
                if ((kinetic_min <= 0.) || (kinetic_max < kinetic_min)) {
                        fprintf(stderr, "pumas-tabulate: error."
                                        " Invalid kinetic energy limits. "
                                        "Call with -h, --help for usage.\n");
                        exit_gracefully(EXIT_FAILURE);
                }

                data.kinetic = malloc(kinetic_n * sizeof(double));
                if (data.kinetic == NULL) {
                        perror("pumas-tabulate");
                        exit_gracefully(EXIT_FAILURE);
                }
                int i;
                const double dlnk = (kinetic_n > 1) ?
                    log(kinetic_max / kinetic_min) / (kinetic_n - 1) :
                    0.;
                for (i = 0; i < kinetic_n; i++)
                        data.kinetic[i] = kinetic_min * exp(dlnk * i);
                data.n_kinetics = kinetic_n;
        }

        /* Tabulate the energy loss. */
        while (data.m_stack != NULL) {
                enum pumas_return rc = _pumas_tabulate(&data);
                if (rc == PUMAS_RETURN_MEMORY_ERROR) {
                        perror("pumas-tabulate");
                        exit_gracefully(EXIT_FAILURE);
                } else if (rc == PUMAS_RETURN_PATH_ERROR) {
                        fprintf(stderr, "pumas-tabulate: error."
                                        " Couldn't write to output file.\n"
                                        "  -> %s\n",
                            data.path);
                        exit_gracefully(EXIT_FAILURE);
                } else if (rc == PUMAS_RETURN_IO_ERROR) {
                        fprintf(stderr,
                            "pumas-tabulate: error."
                            " Output file already exists.\n  -> %s\n"
                            "Call with -f, --force to overwrite it.\n",
                            data.path);
                        exit_gracefully(EXIT_FAILURE);
                }

                /* Drop the processed material. */
                tabulation_material_drop(&data);
        }

        /* Exit to the OS. */
        exit_gracefully(EXIT_SUCCESS);
}
