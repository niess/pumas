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

/* The PUMAS API. */
#include "pumas.h"

#ifdef _WIN32
/* For rand_s on Windows */
#define _CRT_RAND_S
#endif

/* The C standard library. */
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*  For debugging with gdb, on linux. */
#define GDB_MODE 0
#if (GDB_MODE)
#ifndef __USE_GNU
#define __USE_GNU
#endif
#include <fenv.h>
#endif

/* For the versioning. */
#ifndef PUMAS_VERSION
#define PUMAS_VERSION 1.1
#endif

/* Some tuning factors as macros. */
/**
 * Number of schemes to tabulate for the computation of the energy loss.
 *
 * The NO_LOSS scheme has no tabulation and the detailed one uses the same
 * tabulation than the hybrid scheme.
 */
#define N_SCHEMES 2
/**
 * Order of expansion for the computation of the magnetic deflection,
 * when using CSDA in a uniform medium of infinite extension.
 */
#define N_LARMOR_ORDERS 8
/**
 * Number of inelastic processes for which Discrete Energy Losses (DELs) are
 * simulated.
 */
#define N_DEL_PROCESSES 4
/**
 * Default cutoff between Continuous Energy Loss (CEL) and DELs.
 */
#define DEFAULT_CUTOFF 5E-02
/**
 * Exponents of the differential cross section approximation in Backward
 * Monte-Carlo (BMC).
 */
#define RMC_ALPHA_LOW 2.5
#define RMC_ALPHA_HIGH 1.4
/**
 * Ratio of the Coulomb interaction length, restricted to catastrophic events,
 * to the multiple scattering 1st transport path length.
 */
#define EHS_OVER_MSC 1E-04
/**
 * Maximum path length for Elastic Hard Scattering (EHS) events, in kg/m^-2.
 */
#define EHS_PATH_MAX 1E+09
/**
 * Default accuracy (relative step length limit) of the Monte Carlo integration.
 */
#define DEFAULT_ACCURACY 1E-02
/**
 * Maximum deflection angle for a soft scattering event, in degrees.
 */
#define MAX_SOFT_ANGLE 1E+00
/**
 * Minimum step size.
 */
#define STEP_MIN 1E-07
/**
 * Order of the polynomials for the Differential Cross-Section (DCS) model.
 *
 * P is a polynomial of log(nu) while Q is in log(1-nu) with a null constant
 * term.
 */
#define DCS_MODEL_ORDER_P 6
#define DCS_MODEL_ORDER_Q 2
/**
 * Minimum kinetic energy for using the DCS model.
 */
#define DCS_MODEL_MIN_KINETIC 10.
/**
 * Maximum allowed energy tranfer for using the DCS model.
 */
#define DCS_MODEL_MAX_FRACTION 0.95
/**
 * Number of samples for the DCS tabulation.
 */
#define DCS_SAMPLING_N 11

/* Some constants, as macros. */
/**
 * The muon decay length in m.
 */
#define MUON_C_TAU 658.654
/**
 * The tau decay length in m.
 */
#define TAU_C_TAU 87.03E-06
/**
 * Larmor magnetic factor in m^-1 GeV/c T^-1.
 */
#define LARMOR_FACTOR 0.299792458
/**
 * The electron mass in GeV/c^2.
 */
#define ELECTRON_MASS 0.510998910E-03
/**
 * The muon mass in GeV/c^2.
 */
#define MUON_MASS 0.10565839
/**
 * The tau mass in GeV/c^2.
 */
#define TAU_MASS 1.77682
/**
 * The proton mass in GeV/c^2.
 */
#define PROTON_MASS 0.938272
/**
 * The neutron mass in GeV/c^2.
 */
#define NEUTRON_MASS 0.939565
#ifndef M_PI
/**
 * Define pi, if unknown.
 */
#define M_PI 3.14159265358979323846
#endif
/**
 * Avogadro's number
 */
#define AVOGADRO_NUMBER 6.02214076E+23
/**
 * Default Bremsstrahlung model
 */
#define DEFAULT_BREMSSTRAHLUNG "KKP"
/**
 * Default pair production model
 */
#define DEFAULT_PAIR_PRODUCTION "KKP"
/**
 * Default photonuclear model
 */
#define DEFAULT_PHOTONUCLEAR "DRSS"

/* Helper macros for managing errors. */
#define ERROR_INITIALISE(caller)                                               \
        struct error_context error_data = {.code = PUMAS_RETURN_SUCCESS,       \
                .function = (pumas_function_t *)caller };                      \
        struct error_context * error_ = &error_data;

#define ERROR_MESSAGE(rc, message)                                             \
        error_format(error_, rc, __FILE__, __LINE__, message),                 \
            error_raise(error_)

#define ERROR_FORMAT(rc, format, ...)                                          \
        error_format(error_, rc, __FILE__, __LINE__, format, __VA_ARGS__),     \
            error_raise(error_)

#define ERROR_REGISTER(rc, message)                                            \
        error_format(error_, rc, __FILE__, __LINE__, message)

#define ERROR_VREGISTER(rc, format, ...)                                       \
        error_format(error_, rc, __FILE__, __LINE__, format, __VA_ARGS__)

#define ERROR_RAISE() error_raise(error_)

#define ERROR_NULL_PHYSICS()                                                   \
        ERROR_MESSAGE(PUMAS_RETURN_PHYSICS_ERROR,                              \
            "a NULL physics pointer was provided")

#define ERROR_NOT_INITIALISED()                                                \
        ERROR_MESSAGE(PUMAS_RETURN_PHYSICS_ERROR,                              \
            "the Physics has not been initialised")

#define ERROR_INVALID_SCHEME(scheme)                                           \
        ERROR_FORMAT(PUMAS_RETURN_INDEX_ERROR,                                 \
            "invalid energy loss scheme [%d]", scheme)

#define ERROR_INVALID_ELEMENT(element)                                         \
        ERROR_FORMAT(                                                          \
            PUMAS_RETURN_INDEX_ERROR, "invalid element index [%d]", element)

#define ERROR_INVALID_MATERIAL(material)                                       \
        ERROR_FORMAT(                                                          \
            PUMAS_RETURN_INDEX_ERROR, "invalid material index [%d]", material)

#define ERROR_REGISTER_MEMORY()                                                \
        ERROR_REGISTER(PUMAS_RETURN_MEMORY_ERROR, "could not allocate memory")

#define ERROR_REGISTER_EOF(path)                                               \
        ERROR_VREGISTER(PUMAS_RETURN_END_OF_FILE,                              \
            "abnormal end of file when parsing `%s'", path)

#define ERROR_REGISTER_UNEXPECTED_TAG(tag, path, line)                         \
        ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,                             \
            "unexpected XML tag `%s' [@%s:%d]", tag, path, line)

#define ERROR_REGISTER_TOO_LONG(path, line)                                    \
        ERROR_VREGISTER(PUMAS_RETURN_TOO_LONG,                                 \
            "XML node is too long [@%s:%d]", path, line)

#define ERROR_REGISTER_INVALID_XML_VALUE(value, path, line)                    \
        ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,                             \
            "invalid XML value `%s' [@%s:%d]", value, path, line)

#define ERROR_REGISTER_INVALID_XML_ATTRIBUTE(attribute, node, path, line)      \
        ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,                             \
            "invalid XML attribute `%s' for element %s [@%s:%d]", attribute,   \
            node, path, line)

#define ERROR_REGISTER_NEGATIVE_DENSITY(material)                              \
        ERROR_VREGISTER(                                                       \
            PUMAS_RETURN_DENSITY_ERROR, "negative density for `%s'", material)

/* Function prototypes for the DCS implementations. */
/**
 * Handle for a DCS computation function.
 */
struct atomic_element;
typedef double(dcs_function_t)(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q);
/**
 * Handle for a polar angle sampling function.
 */
typedef double(polar_function_t)(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf);

/* A collection of low level flags. */
/**
 * Node keys for the MDF files.
 */
enum mdf_key {
        /** An atomic component. */
        MDF_KEY_ATOMIC_COMPONENT,
        /** A composite material. */
        MDF_KEY_COMPOSITE,
        /** A component of a composite material. */
        MDF_KEY_COMPOSITE_COMPONENT,
        /** The root PUMAS node. */
        MDF_KEY_PUMAS,
        /** An atomic element. */
        MDF_KEY_ELEMENT,
        /** A base material. */
        MDF_KEY_MATERIAL,
        /** Any other node. */
        MDF_KEY_OTHER
};
/**
 * Nodes hierarchy in MDFs.
 */
enum mdf_depth {
        /** Outside of the root PUMAS node. */
        MDF_DEPTH_EXTERN = 0,
        /** Inside the root pumas node. */
        MDF_DEPTH_ROOT,
        /** Inside an element node. */
        MDF_DEPTH_ELEMENT,
        /** Inside a material node. */
        MDF_DEPTH_MATERIAL,
        /** Inside a composite material node. */
        MDF_DEPTH_COMPOSITE
};
/**
 * Tags for operations relative to the parsing of materials in MDFs.
 */
enum mdf_settings_operation {
        /** Free the names table. */
        MDF_INDEX_FREE = -1,
        /** Initialise the names table. */
        MDF_INDEX_INITIALISE,
        /** Add a new base material to the table. */
        MDF_INDEX_WRITE_MATERIAL,
        /** Finalise a base material info. */
        MDF_INDEX_FINALISE_MATERIAL,
        /** Add a new composite material to the table. */
        MDF_INDEX_INITIALISE_COMPOSITE,
        /** Update a composite material with a new component. */
        MDF_INDEX_UPDATE_COMPOSITE,
        /** Finalise the composite material. */
        MDF_INDEX_FINALISE_COMPOSITE
};

/* Low level data structures. */
/**
 * Temporary data for the computation of the Coulomb scattering.
 *
 * The DCS is given by a product of 3 Wentzel distributions with 1 atomic
 * screening parameter and 2 nuclear screening ones. A pole reduction allows
 * the analytical computation of the various momenta: total cross-section, first
 * transport path length, ...
 */
struct coulomb_data {
        /** The partial cross section for hard events. */
        double cs_hard;
        /** The inverse of the mean free grammage. */
        double invlambda;
        /** The atomic and nuclear screening parameters. */
        double screening[3];
        /** The 1st order coefficients of the pole reduction of the DCS. */
        double a[3];
        /** The 2nd order coefficients of the pole reduction of the DCS. */
        double b[3];
        /** The spin correction factor. */
        double fspin;
        /**
         * The parameters of the relativistic transform from Center of Mass (CM)
         * frame to the Laboratory one.
         */
        double fCM[2];
};
/**
 * Temporary data for the rendering of a Coulomb DEL. The event is randomised
 * using the inverse transform method. A root bracketing algorithm is used,
 * seeded with an approximate solution.
 */
struct coulomb_workspace {
        /** The index of the propagation material. */
        int material;
        /** The index of the hard scatterer element. */
        int ihard;
        /** The targeted cross section. */
        double cs_h;
        /** Placeholder for the table of Coulomb scattering data. */
        struct coulomb_data data[];
};
/**
 * Low level container for the local properties of a propagation medium.
 */
struct medium_locals {
        /** The public API properties exposed to the end user. */
        struct pumas_locals api;
        /** A flag telling if the material has a magnetic field or not. */
        int magnetized;
        /** The local's physics */
        const struct pumas_physics * physics;
};
/**
 *  Data for the default per context PRNG
 */
struct pumas_random_data {
/*
 * Version tag for the random data format. Increment whenever the
 * structure changes.
 */
#define RANDOM_BINARY_DUMP_TAG 0
        /** The initial seed */
        unsigned long seed;
        /** Index in the PRNG buffer */
        int index;
#define MT_PERIOD 624
        /** PRNG buffer (Mersenne Twister) */
        unsigned long buffer[MT_PERIOD];
};
/**
 * The local data managed by a simulation context.
 */
struct simulation_context {
        /** The public API settings exposed to the end user. */
        struct pumas_context api;
        /** Handle for Physics tables. */
        const struct pumas_physics * physics;
        /** Lifetime limit for the decay process. */
        double lifetime;
        /**
         * Last kinetic energy indices used in the tables.
         *
         * We keep a memory of the last indices accessed in tables. The memory
         * has a depth of 2 which should be a good compromise since most use
         * cases consider differences between an initial and final state.
         */
        int index_K_last[2];
        /** Last grammage indices used in the tables. */
        int index_X_last[2];
        /** Flag for the first step, for integration of various quantities. */
        int step_first;
        /** Tracking of stepping events. */
        enum pumas_event step_event;
        /** The expected next event during the stepping. */
        enum pumas_event step_foreseen;
        /** The kinetic limit converted to grammage. */
        double step_X_limit;
        /** The scaterring 1st transport path length of the previous step. */
        double step_invlb1;
        /** The larmor radius of the previous step. */
        double step_rLarmor;
        /** The magnetic transverse direction of the previous step. */
        double step_uT[3];
        /** Data for the default PRNG. */
        struct pumas_random_data * random_data;
        /** Flag for the parity check of the Gaussian random generator.
         *
         * Gaussian variates are generated in pair using the Box-Muller
         * transform.
         */
        int randn_done;
        /** The next Gaussian variate. */
        double randn_next;
        /**
         * Pointer to the worspace for the temporary storage of intermediary
         * computations.
         */
        struct coulomb_workspace * workspace;
        /** Size of the user extended memory. */
        int extra_memory;
        /**
         * Placeholder for variable data storage with -double- memory alignment.
         *
         * Extra bytes are allocated for the workspace, the error stack and
         * any extended memory for end user usage.
         */
        double data[];
};
/**
 * Handle for a stack of recorded frames.
 *
 * The frames are allocated in bunches. The bunches are managed as a chained
 * list.
 */
struct frame_stack {
        /** The memory size left, in bytes. */
        int size;
        /** Pointer to the next memory segment. */
        struct frame_stack * next;
        /** Pointer to the first frame in the stack. */
        struct pumas_frame * frame;
        /** Placeholder for frames */
        struct pumas_frame frames[];
};
/**
 * Low level container for a frame recorder.
 */
struct frame_recorder {
        /** The public API data exposed to the end user. */
        struct pumas_recorder api;
        /** Link to the last record. */
        struct pumas_frame * last;
        /** Link to the 1st entry of the chained list of stacks. */
        struct frame_stack * stack;
        /** Placeholder for extra data. */
        double data[];
};
/**
 * Data relative to an atomic element.
 */
struct atomic_element {
        /** The element atomic number. */
        double Z;
        /** The element atomic mass, in g/mol. */
        double A;
        /** The element Mean Excitation Energy (MEE), in eV. */
        double I;
        /** The element name. */
        char * name;
        /** Tabulation for the DCS model: polynomials and precomputed values. */
        float * dcs_data;
        /** Placeholder for user data with -double- memory alignment. */
        double data[];
};
/**
 * An atomic component of a material.
 */
struct material_component {
        /** The atomic element index. */
        int element;
        /** The element mass fraction in the material. */
        double fraction;
};
/**
 * A component of a composite material.
 */
struct composite_component {
        /** The constituent base material index. */
        int material;
        /** The component mass fraction in the composite. */
        double fraction;
};
/**
 * Handle for a composite material.
 */
struct composite_material {
        /** The number of sub components. */
        int n_components;
        /** Placeholder for the sub components' data. */
        struct composite_component component[];
};
/**
 * Temporary data for the parsing of a MDF.
 */
struct mdf_buffer {
        /** Flag for the dry initialisation mode. */
        int dry_mode;
        /** Handle to the MDF. */
        FILE * fid;
        /** Current position in the read buffer. */
        char * pos;
        /** The number of bytes left to be read. */
        int left;
        /** The total size of the read buffer. */
        int size;
        /** The current line number in the MDF file. */
        int line;
        /** Current node depth during the MDF parsing. */
        enum mdf_depth depth;
        /** Counter for the number of base materials in a composite. */
        int materials_in;
        /** Counter for the number of elements in a material. */
        int elements_in;
        /** Pointer to the current MDF. */
        const char * mdf_path;
        /** The number of kinetic energy rows in a dE/dX file. */
        int n_energies;
        /** The total number of materials, base and composites. */
        int n_materials;
        /** The number of composite materials. */
        int n_composites;
        /** The number of atomic elements. */
        int n_elements;
        /** The header length of dE/dX files. */
        int n_energy_loss_header;
        /** The total number of atomic element components. */
        int n_components;
        /** The maximum number of atomic elements in a single material. */
        int max_components;
        /** The total byte size for the storage of composite materials. */
        int size_composite;
        /** The total size of the path to the dE/dX files. */
        int size_dedx_path;
        /** The total size of elements names. */
        int size_elements_names;
        /** The total size of materials names. */
        int size_materials_names;
        /**
         * The offset for the tabulations of DCS and of the coefficients of the
         * polynomial approximations.
         */
        int dcs_model_offset;
        /** Placeholder for the read buffer. */
        char data[];
};
/*!
 * Pointers to the data fields of a node in a MDF file.
 *
 * The pointers refer to an mdf_buffer object. If the buffer is refilled the
 * links are no mode valid.
 */
struct mdf_node {
        /** The node key. */
        enum mdf_key key;
        /** Flag telling if this is a head node. */
        int head;
        /** Flag telling if this is a tail node. */
        int tail;
        /** The first attribute names. */
        union attribute1 {
                /** The node name. */
                char * name;
        } at1;
        /** The second attribute names. */
        union attribute2 {
                /** The energy loss file name. */
                char * file;
                /** The atomic number. */
                char * Z;
                /** The mass fraction. */
                char * fraction;
        } at2;
        /** The third attribute names. */
        union attribute3 {
                /** The atomic mass. */
                char * A;
                /** The composite component's density. */
                char * density;
        } at3;
        /** The fourth attribute names. */
        union attribute4 {
                /** The mean excitation energy (MEE). */
                char * I;
        } at4;
};
/**
 * Temporary data for a DEL event.
 */
struct del_info {
        union {
                /** The energy transfer. */
                double Q;
                /** The Monte-Carlo weight factor. */
                double weight;
        }
        /** Data specific to reverse Monte-Carlo. */
        reverse;
        /** The index of the inelastic sub-process. */
        int process;
        /** The index of the target element. */
        int element;
};
/**
 * Global data shared by all simulation contexts.
 */
struct pumas_physics {
/*
 * Version tag for the physics data format. Increment whenever the
 * structure changes.
 */
#define PHYSICS_BINARY_DUMP_TAG 5

        /** The total byte size of the shared data. */
        int size;
        /** The number of kinetic energy values in the dE/dX tables. */
        int n_energies;
        /** The total number of materials, basic and composites. */
        int n_materials;
        /** The total number of composite materials. */
        int n_composites;
        /** The total number of declared atomic_elements. */
        int n_elements;
        /** The total number of atomic components in materials. */
        int n_components;
        /** The maximum number of atomic components in a single material. */
        int max_components;
        /** The number of header lines in a dE/dX table. */
        int n_energy_loss_header;
        /**
         * Offset for the tabulation of DCS and their polynomial
         * approximations.
         */
        int dcs_model_offset;
        /** The transported particle type. */
        enum pumas_particle particle;
        /** The transported particle decay length, in m. */
        double ctau;
        /** The transported particle rest mass, in GeV. */
        double mass;
        /** The relative cutoff between CEL and DELs. */
        double cutoff;
        /** Path to the current MDF. */
        char * mdf_path;
        /** Path where the dE/dX files are stored. */
        char * dedx_path;
        /** Names of the dE/dX files. */
        char ** dedx_filename;
        /** The tabulated values of the kinetic energy. */
        double * table_K;
        /** The tabulated values of the total grammage (CSDA range). */
        double * table_X;
        /** The tabulated values of the total proper time. */
        double * table_T;
        /** The tabulated values of the average energy loss. */
        double * table_dE;
        /** The tabulated values of the EHS number of interaction lengths. */
        double * table_NI_el;
        /**
         * The tabulated values of the number of interaction lengths for
         * inelastic DELs.
         */
        double * table_NI_in;
        /** The tabulated cross section values. */
        double * table_CS;
        /** The tabulated cross section fractions. */
        double * table_CSf;
        /** The tabulated cross section normalisation. */
        double * table_CSn;
        /** The element wise fractional threshold for DELs. */
        double * table_Xt;
        /** The total kinetic threshold for DELs. */
        double * table_Kt;
        /**
         * The tabulated values of the magnetic deflection momenta, within
         * the CSDA.
         */
        double * table_Li;
        /** The last tabulated value of the ionisation energy loss. */
        double * table_a_max;
        /** The last tabulated value of the radiative energy loss. */
        double * table_b_max;
        /**
         * The tabulated angular cutoff values for the splitting of Coulomb
         * scattering.
         */
        double * table_Mu0;
        /** The tabulated interaction lengths for DEL Coulomb events. */
        double * table_Lb;
        /** The tabulated multiple scattering 1st moment. */
        double * table_Ms1;
        /** The number of elements in a material. */
        int * elements_in;
        /** The reference density of a material. */
        double * material_density;
        /** The relative electronic density of a material. */
        double * material_ZoA;
        /** The mean excitation energy of a base material. */
        double * material_I;
        /** The density effect parameters of a base material. */
        struct pumas_physics_density_effect * material_density_effect;
        /** The properties of an atomic element . */
        struct atomic_element ** element;
        /** The composition of a base material. */
        struct material_component ** composition;
        /** The composition of a composite material. */
        struct composite_material ** composite;
        /** The material names. */
        char ** material_name;
        /** The Bremsstrahlung model. */
        char * model_bremsstrahlung;
        /** The pair_production model. */
        char * model_pair_production;
        /** The photonuclear model. */
        char * model_photonuclear;
        /** The Bremsstrahlung DCS. */
        pumas_dcs_t * dcs_bremsstrahlung;
        /** The pair_production DCS. */
        pumas_dcs_t * dcs_pair_production;
        /** The photonuclear DCS. */
        pumas_dcs_t * dcs_photonuclear;
        /**
         * Placeholder for shared data storage with -double- memory alignment.
         */
        double data[];
};

struct error_context {
        enum pumas_return code;
        pumas_function_t * function;
#define ERROR_MSG_LENGTH 1024
        char message[ERROR_MSG_LENGTH];
};

/**
 * Default error handler.
 */
static void default_error_handler(
    enum pumas_return rc, pumas_function_t * caller, const char * message)
{
        /* Dump the error summary */
        fputs("pumas: library error. See details below\n", stderr);
        fprintf(stderr, "error: %s\n", message);

        /* Exit to the OS */
        exit(EXIT_FAILURE);
}

/**
 * Shared data for the error handling.
 */
static struct {
        pumas_handler_cb * handler;
        int catch;
        struct error_context catch_error;
} s_error = { &default_error_handler, 0 };

/* Prototypes of low level static functions. */
/**
 * Encapsulations of the tabulated CEL and DEL properties.
 */
static double cel_grammage(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double kinetic);
static double cel_grammage_as_time(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double time);
static double cel_proper_time(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double kinetic);
static double cel_kinetic_energy(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double grammage);
static double cel_energy_loss(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double kinetic);
static double cel_magnetic_rotation(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic);
static double del_cross_section(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic);
static double del_interaction_length(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic);
static double del_kinetic_from_interaction_length(
    const struct pumas_physics * physics, struct pumas_context * context,
    int material, double nI);
static double ehs_interaction_length(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double kinetic);
static double ehs_kinetic_from_interaction_length(
    const struct pumas_physics * physics, struct pumas_context * context,
    enum pumas_mode scheme, int material, double nI);
/**
 * Routines related to DCS: implementation and handling.
 */
static inline dcs_function_t * dcs_get(int process);
static inline int dcs_get_index(dcs_function_t * dcs_func);
static double dcs_bremsstrahlung(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q);
static double dcs_pair_production(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q);
static double dcs_photonuclear(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q);
static inline int dcs_photonuclear_check(double K, double q);
static double dcs_ionisation(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q);
static double dcs_ionisation_integrate(const struct pumas_physics * physics,
    int mode, const struct atomic_element * element, double K, double xlow);
static double dcs_ionisation_randomise(const struct pumas_physics * physics,
    struct pumas_context * context, const struct atomic_element * element,
    double K, double xlow);
static double dcs_evaluate(const struct pumas_physics * physics,
    struct pumas_context * context, dcs_function_t * dcs_func,
    const struct atomic_element * element, double K, double q);
static void dcs_model_fit(int m, int n, const double * x, const double * y,
    const double * w, double * c);

/**
 * Implementations of polar angle distributions and accessor.
 */
static inline polar_function_t * polar_get(int process);
static double polar_bremsstrahlung(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf);
static double polar_pair_production(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf);
static double polar_photonuclear(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf);
static double polar_ionisation(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf);
/**
 * Low level routines for the propagation in matter.
 */
static enum pumas_event transport_with_csda(
    const struct pumas_physics * physics, struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    struct medium_locals * locals, struct error_context * error_);
static enum pumas_return transport_csda_deflect(
    const struct pumas_physics * physics, struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    struct medium_locals * locals, double ki, double distance,
    struct error_context * error_);
static enum pumas_return csda_magnetic_transport(
    const struct pumas_physics * physics, struct pumas_context * context,
    int material, double density, double magnet, double charge, double kinetic,
    double phase, double * x, double * y, double * z,
    struct error_context * error_);
static enum pumas_event transport_with_stepping(
    const struct pumas_physics * physics, struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr,
    struct medium_locals * locals, double step_max_medium,
    enum pumas_step step_max_type, double step_max_locals,
    struct error_context * error_);
static double transport_set_locals(const struct pumas_context * context,
    struct pumas_medium * medium, struct pumas_state * state,
    struct medium_locals * locals);
static void transport_limit(const struct pumas_physics * physics,
    struct pumas_context * context, const struct pumas_state * state,
    int material, double di, double Xi, double * distance_max);
static void transport_do_del(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material);
static void transport_do_ehs(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material);
/**
 * Low level routines for randomising DELs.
 */
static polar_function_t * del_randomise_forward(
    const struct pumas_physics * physics, struct pumas_context * context,
    struct pumas_state * state, int material, int * process);
static polar_function_t * del_randomise_reverse(
    const struct pumas_physics * physics, struct pumas_context * context,
    struct pumas_state * state, int material, int * process);
static void del_randomise_power_law(struct pumas_context * context,
    double alpha, double xmin, double xmax, double * p_r, double * p_w);
static void del_randomise_ziggurat(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state,
    dcs_function_t * dcs_func, const struct atomic_element * element,
    double xmin, double xmax, float * cdf_sampling);
static void del_randomise_target(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material,
    struct del_info * info);
/**
 * Helper routine for recording a state.
 */
static void record_state(struct pumas_context * context,
    struct pumas_medium * medium, enum pumas_event event,
    struct pumas_state * state);
/**
 * For memory padding.
 */
static int memory_padded_size(int size, int pad_size);
/**
 * For error handling.
 */
static enum pumas_return error_raise(struct error_context * context);
static enum pumas_return error_format(struct error_context * context,
    enum pumas_return rc, const char * file, int line, const char * format,
    ...);
/**
 * Routines for the Coulomb scattering and Transverse Transport (TT).
 */
static void coulomb_screening_parameters(const struct pumas_physics * physics,
    struct pumas_context * context, double kinetic, int element,
    double * screening);
static double coulomb_wentzel_path(const struct pumas_physics * physics,
    double kinetic, double Z, double A, double screening);
static double coulomb_ehs_length(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic);
static double coulomb_spin_factor(
    const struct pumas_physics * physics, double kinetic);
static void coulomb_frame_parameters(const struct pumas_physics * physics,
    double kinetic, double Ma, double * kinetic0, double * parameters);
static void coulomb_pole_decomposition(
    double * screening, double * a, double * b);
static double coulomb_restricted_cs(
    double mu0, double fspin, double * screening, double * a, double * b);
static void coulomb_transport_coefficients(double mu, double fspin,
    double * screening, double * a, double * b, double * coefficient);
static double transverse_transport_ionisation(
    const struct pumas_physics * physics, const struct atomic_element * element,
    double kinetic);
static double transverse_transport_photonuclear(
    const struct pumas_physics * physics, const struct atomic_element * element,
    double kinetic);
/**
 * Routines for handling tables: interpolation and utility accessors.
 */
static void table_bracket(
    const double * table, double value, int * p1, int * p2);
static int table_index(const struct pumas_physics * physics,
    struct pumas_context * context, const double * table, double value);
static double table_interpolate(const struct pumas_physics * physics,
    struct pumas_context * context, const double * table_X,
    const double * table_Y, double x);
static void table_get_msc(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic, double * mu0,
    double * invlb1);
static inline double * table_get_K(
    const struct pumas_physics * physics, int row);
static inline double * table_get_X(
    const struct pumas_physics * physics, int scheme, int material, int row);
static inline double * table_get_T(
    const struct pumas_physics * physics, int scheme, int material, int row);
static inline double * table_get_dE(
    const struct pumas_physics * physics, int scheme, int material, int row);
static inline double * table_get_NI_el(
    const struct pumas_physics * physics, int scheme, int material, int row);
static inline double * table_get_NI_in(
    const struct pumas_physics * physics, int material, int row);
static inline double * table_get_CS(
    const struct pumas_physics * physics, int material, int row);
static inline double * table_get_CSf(
    const struct pumas_physics * physics, int process, int component, int row);
static inline double * table_get_CSn(
    const struct pumas_physics * physics, int process, int element, int row);
static inline double * table_get_Xt(
    const struct pumas_physics * physics, int process, int element, int row);
static inline double * table_get_Kt(
    const struct pumas_physics * physics, int material);
static inline double * table_get_cel(const struct pumas_physics * physics,
    int process, int element, int row, double * table);
static inline double * table_get_Li(
    const struct pumas_physics * physics, int material, int order, int row);
static inline double * table_get_a_max(
    const struct pumas_physics * physics, int material);
static inline double * table_get_b_max(
    const struct pumas_physics * physics, int scheme, int material);
static inline double * table_get_Mu0(
    const struct pumas_physics * physics, int material, int row);
static inline double * table_get_Lb(
    const struct pumas_physics * physics, int material, int row);
static inline double * table_get_Ms1(
    const struct pumas_physics * physics, int material, int row);
static inline double * table_get_ms1(
    const struct pumas_physics * physics, int element, int row, double * table);
static inline float * table_get_dcs_coeff(const struct pumas_physics * physics,
    const struct atomic_element * element, int process, int kinetic);
static inline float * table_get_dcs_value(const struct pumas_physics * physics,
    const struct atomic_element * element, int process, int kinetic);
/**
 * Low level routines for the stepping.
 */
static enum pumas_return step_transport(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int straight,
    struct pumas_medium * medium, struct medium_locals * locals,
    double grammage_max, double step_max_medium, enum pumas_step step_max_type,
    double * step_max_locals, struct pumas_medium ** out_medium,
    struct error_context * error_);
static void step_fluctuate(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material,
    double Xtot, double dX, double * kf, double * dE);
static double step_fluctuations2(
    const struct pumas_physics * physics, int material, double kinetic);
static double step_randn(struct pumas_context * context);
static double transport_hard_coulomb_objective(
    const struct pumas_physics * physics, double mu, void * parameters);
static void step_rotate_direction(struct pumas_context * context,
    struct pumas_state * state, double cos_theta);
/**
 * I/O utility routines.
 */
static enum pumas_return io_parse_dedx_file(struct pumas_physics * physics,
    FILE * fid, int material, const char * filename,
    struct error_context * error_);
static enum pumas_return io_parse_dedx_row(struct pumas_physics * physics,
    char * buffer, int material, int * row, const char * filename, int line,
    struct error_context * error_);
static enum pumas_return io_read_line(FILE * fid, char ** buffer,
    const char * filename, int line, struct error_context * error_);
/**
 * Routines for the parsing of MDFs.
 */
static enum pumas_return mdf_parse_settings(
    const struct pumas_physics * physics, struct mdf_buffer * mdf,
    const char * dedx_path, struct error_context * error_);
static int mdf_settings_index(
    int operation, int value, struct error_context * error_);
static int mdf_settings_name(
    int size, char prefix, const char * name, struct error_context * error_);
static enum pumas_return mdf_parse_kinetic(
    struct mdf_buffer * mdf, const char * path, struct error_context * error_);
static enum pumas_return mdf_parse_elements(
    const struct pumas_physics * physics, struct mdf_buffer * mdf,
    struct error_context * error_);
static enum pumas_return mdf_parse_materials(struct pumas_physics * physics,
    struct mdf_buffer * mdf, struct error_context * error_);
static enum pumas_return mdf_parse_composites(struct pumas_physics * physics,
    struct mdf_buffer * mdf, struct error_context * error_);
static enum pumas_return mdf_get_node(struct mdf_buffer * mdf,
    struct mdf_node * node, struct error_context * error_);
static enum pumas_return mdf_skip_pattern(struct mdf_buffer * mdf,
    const char * pattern, struct error_context * error_);
static enum pumas_return mdf_format_path(const char * directory,
    const char * mdf_path, char ** filename, int * offset_dir, int * size_name,
    struct error_context * error_);
/**
 * Routines for the pre-computation of various properties: CEL, DCS, ...
 */
static enum pumas_return compute_composite(struct pumas_physics * physics,
    int material, struct error_context * error_);
static enum pumas_return compute_composite_density(
    struct pumas_physics * physics, int material,
    struct error_context * error_);
static void compute_composite_weights(
    struct pumas_physics * physics, int material);
static void compute_composite_tables(
    struct pumas_physics * physics, int material);
static void compute_cel_integrals(struct pumas_physics * physics, int imed);
static void compute_kinetic_integral(
    struct pumas_physics * physics, double * table);
static void compute_time_integrals(
    struct pumas_physics * physics, int material);
static void compute_cel_grammage_integral(
    struct pumas_physics * physics, int scheme, int material);
static void compute_csda_magnetic_transport(
    struct pumas_physics * physics, int imed);
static enum pumas_return compute_coulomb_parameters(
    struct pumas_physics * physics, int medium_index, int row,
    struct error_context * error_);
static enum pumas_return compute_coulomb_soft(struct pumas_physics * physics,
    int row, double ** data, struct error_context * error_);
static double compute_cutoff_objective(
    const struct pumas_physics * physics, double mu, void * workspace);
static double * compute_cel_and_del(struct pumas_physics * physics, int row);
static void compute_regularise_del(
    struct pumas_physics * physics, int material);
static double compute_dcs_integral(struct pumas_physics * physics, int mode,
    const struct atomic_element * element, double kinetic, dcs_function_t * dcs,
    double xlow, int nint);
static void compute_ZoA(struct pumas_physics * physics, int material);
static enum pumas_return compute_dcs_model(struct pumas_physics * physics,
    dcs_function_t * dcs_func, struct atomic_element * element,
    struct error_context * error_);
/**
 * Helper function for mapping an atomic element from its name.
 */
static int element_index(
    const struct pumas_physics * physics, const char * name);
/**
 * Helper function for mapping a material from its name.
 */
static enum pumas_return material_index(const struct pumas_physics * physics,
    const char * material, int * index, struct error_context * error_);
/**
 * Various math utilities, for integration, root finding and SVD.
 */
static int math_find_root(
    double (*f)(const struct pumas_physics * physics, double x, void * params),
    const struct pumas_physics * physics, double xa, double xb,
    const double * fa_p, const double * fb_p, double xtol, double rtol,
    int iter, void * params, double * x0);
static int math_gauss_quad(int n, double * p1, double * p2);
static int math_svd(
    int m, int n, double * a, double * w, double * v, double * work);
static void math_svdsol(int m, int n, double * b, double * u, double * w,
    double * v, double * work, double * x);
static double math_rms(double a, double b);

/* Below is the implementation of the public API functions. See pumas.h for a
 * detailed description of each function.
 */
#if (GDB_MODE)
/**
 * A flag for floating point exceptions.
 */
static int fe_status;
#endif

/* Getter for library constants */
PUMAS_API enum pumas_return pumas_constant(
    enum pumas_constant index, double * value)
{
        ERROR_INITIALISE(pumas_constant);

        const double values[] = {AVOGADRO_NUMBER, ELECTRON_MASS, MUON_C_TAU,
            MUON_MASS, NEUTRON_MASS, PROTON_MASS, TAU_C_TAU, TAU_MASS};

        if (value == NULL) {
                return ERROR_MESSAGE(PUMAS_RETURN_VALUE_ERROR,
                    "NULL value pointer");
        }

        if ((index < 0) || (index >= PUMAS_N_CONSTANTS)) {
                *value = 0;
                return ERROR_FORMAT(PUMAS_RETURN_INDEX_ERROR,
                    "invalid `constant' index [%d]", index);
        }

        *value = values[index];

        return PUMAS_RETURN_SUCCESS;
}

/*
 * Public library functions: user supplied memory allocation.
 */
static pumas_allocate_cb * allocate = malloc;
static pumas_reallocate_cb * reallocate = realloc;
static pumas_deallocate_cb * deallocate = free;

/**
 * Set the memory allocation function for the PUMAS library.
 */
void pumas_memory_allocator(pumas_allocate_cb * allocator)
{
        allocate = (allocator == NULL) ? malloc : allocator;
}

/**
 * Set the memory reallocation function for the PUMAS library.
 */
void pumas_memory_reallocator(pumas_reallocate_cb * reallocator)
{
        reallocate = (reallocator == NULL) ? realloc : reallocator;
}

/**
 * Set the memory deallocation function for the PUMAS library.
 */
void pumas_memory_deallocator(pumas_deallocate_cb * deallocator)
{
        deallocate = (deallocator == NULL) ? free : deallocator;
}

/*
 * Public library functions: initialisation and termination.
 */

/* Routine for checking the validity of a DCS model name (forward decl.) */
static enum pumas_return dcs_check_model(enum pumas_process process,
     const char * model, struct error_context * error_);

/*
 * Low level initialisation. If *dry_mode* is not null the energy loss tables
 * are not loaded and processed.
 */
static enum pumas_return _initialise(struct pumas_physics ** physics_ptr,
    enum pumas_particle particle, const char * mdf_path, const char * dedx_path,
    int dry_mode, const struct pumas_physics_settings * settings_)
{
        ERROR_INITIALISE(pumas_physics_create);
        if (dry_mode) {
                error_data.function =
                    (pumas_function_t *)&pumas_physics_create_tabulation;
        }

        /* Check if the Physics pointer is NULL. */
        if (physics_ptr == NULL) {
                return ERROR_NULL_PHYSICS();
        }
        *physics_ptr = NULL;
#if (GDB_MODE)
        /* Save the floating points exceptions status and enable them. */
        fe_status = fegetexcept();
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
        FILE * fid_mdf = NULL;
        struct mdf_buffer * mdf = NULL;
        const int pad_size = sizeof(*((*physics_ptr)->data));
#define N_DATA_POINTERS 32
        int size_data[N_DATA_POINTERS];

        /* Check the particle type. */
        if ((particle != PUMAS_PARTICLE_MUON) &&
            (particle != PUMAS_PARTICLE_TAU)) {
                return ERROR_FORMAT(PUMAS_RETURN_UNKNOWN_PARTICLE,
                    "invalid particle index `%d'", (int)particle);
        }

        /* Check and unpack any extra settings */
        struct pumas_physics_settings opts = {
                DEFAULT_CUTOFF,
                DEFAULT_BREMSSTRAHLUNG,
                DEFAULT_PAIR_PRODUCTION,
                DEFAULT_PHOTONUCLEAR
        };
        if (settings_ != NULL) {
                if (settings_->cutoff >= 1) {
                        return ERROR_FORMAT(PUMAS_RETURN_CUTOFF_ERROR,
                            "bad cutoff value (expected a value in ]0, 1[, "
                            " got %g)", settings_->cutoff);
                } else if (settings_->cutoff > 0) {
                        opts.cutoff = settings_->cutoff;
                }

                if (settings_->bremsstrahlung != NULL) {
                        if (dcs_check_model(PUMAS_PROCESS_BREMSSTRAHLUNG,
                            settings_->bremsstrahlung, error_) ==
                            PUMAS_RETURN_SUCCESS) {
                                opts.bremsstrahlung =
                                    settings_->bremsstrahlung;
                        } else {
                                return ERROR_RAISE();
                        }
                }

                if (settings_->pair_production != NULL) {
                        if (dcs_check_model(PUMAS_PROCESS_PAIR_PRODUCTION,
                            settings_->pair_production, error_) ==
                            PUMAS_RETURN_SUCCESS) {
                                opts.pair_production =
                                    settings_->pair_production;
                        } else {
                                return ERROR_RAISE();
                        }
                }

                if (settings_->photonuclear != NULL) {
                        if (dcs_check_model(PUMAS_PROCESS_PHOTONUCLEAR,
                            settings_->photonuclear, error_) ==
                            PUMAS_RETURN_SUCCESS) {
                                opts.photonuclear =
                                    settings_->photonuclear;
                        } else {
                                return ERROR_RAISE();
                        }
                }
        }

        /* Check the path to energy loss tables. */
        struct pumas_physics * physics = NULL;
        if (!dry_mode) {
                if (dedx_path == NULL) dedx_path = getenv("PUMAS_DEDX");
                if (dedx_path == NULL) {
                        ERROR_REGISTER(PUMAS_RETURN_UNDEFINED_DEDX,
                            "missing path to energy loss tables");
                        goto clean_and_exit;
                }
        }

        /* Parse the MDF. */
        const int size_mdf = 2048;
        const char * file_mdf =
            (mdf_path != NULL) ? mdf_path : getenv("PUMAS_MDF");
        if (file_mdf == NULL) {
                ERROR_REGISTER(PUMAS_RETURN_UNDEFINED_MDF,
                    "missing materials description file");
                goto clean_and_exit;
        }
        int size_path = strlen(file_mdf) + 1;
        fid_mdf = fopen(file_mdf, "r");
        if (fid_mdf == NULL) {
                ERROR_VREGISTER(PUMAS_RETURN_PATH_ERROR,
                    "could not open MDF file `%s'", mdf_path);
                goto clean_and_exit;
        }
        mdf = allocate(size_mdf);
        if (mdf == NULL) {
                ERROR_REGISTER_MEMORY();
                goto clean_and_exit;
        }
        mdf->dry_mode = dry_mode;
        mdf->mdf_path = file_mdf;
        mdf->fid = fid_mdf;
        mdf->size = size_mdf - sizeof(*mdf);
        if ((mdf_parse_settings(*physics_ptr, mdf, dedx_path, error_)) !=
            PUMAS_RETURN_SUCCESS)
                goto clean_and_exit;

        /* Backup the parsed settings. */
        struct mdf_buffer settings;
        memcpy(&settings, mdf, sizeof(settings));

        /* Compute the memory mapping. */
        int imem = 0;
        /* mdf_path. */
        size_data[imem++] =
            memory_padded_size(sizeof(char) * size_path, pad_size);
        /* dedx_path. */
        size_data[imem++] = memory_padded_size(
            sizeof(char) * settings.size_dedx_path, pad_size);
        /* dedx_filename. */
        size_data[imem++] = memory_padded_size(
            sizeof(char *) * (settings.n_materials - settings.n_composites),
            pad_size);
        /* table_K. */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_energies, pad_size);
        /* table_X. */
        size_data[imem++] = memory_padded_size(sizeof(double) * N_SCHEMES *
                settings.n_materials * settings.n_energies,
            pad_size);
        /* table_T. */
        size_data[imem++] = memory_padded_size(sizeof(double) * N_SCHEMES *
                settings.n_materials * settings.n_energies,
            pad_size);
        /* table_dE. */
        size_data[imem++] = memory_padded_size(sizeof(double) * N_SCHEMES *
                settings.n_materials * settings.n_energies,
            pad_size);
        /* table_NI_el. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * 2 * settings.n_materials * settings.n_energies,
            pad_size);
        /* table_NI_in. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_energies,
            pad_size);
        /* table_CS. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_energies,
            pad_size);
        /* table_CSf. */
        size_data[imem++] = memory_padded_size(sizeof(double) *
                N_DEL_PROCESSES * settings.n_components * settings.n_energies,
            pad_size);
        /* table_CSn. */
        size_data[imem++] = memory_padded_size(sizeof(double) *
                N_DEL_PROCESSES * settings.n_elements * settings.n_energies,
            pad_size);
        /* table_Xt */
        size_data[imem++] = memory_padded_size(sizeof(double) *
                N_DEL_PROCESSES * settings.n_elements * settings.n_energies,
            pad_size);
        /* table_Kt. */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials, pad_size);
        /* table_Li. */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials *
                    (N_LARMOR_ORDERS + 1) * settings.n_energies,
                pad_size);
        /* table_a_max. */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials, pad_size);
        /* table_b_max. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * N_SCHEMES * settings.n_materials, pad_size);
        /* table_Mu0. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_energies,
            pad_size);
        /* table_Lb. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_energies,
            pad_size);
        /* table_Ms1. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_energies,
            pad_size);
        /* elements_in. */
        size_data[imem++] =
            memory_padded_size(sizeof(int) * settings.n_materials, pad_size);
        /* material_density */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials, pad_size);
        /* material_ZoA */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials, pad_size);
        /* material_I */
        size_data[imem++] = memory_padded_size(sizeof(double) *
            (settings.n_materials - settings.n_composites), pad_size);
        /* material_density_effect */
        size_data[imem++] = memory_padded_size(
            sizeof(struct pumas_physics_density_effect) *
            (settings.n_materials - settings.n_composites), pad_size);
        /* element. */
        const int n_dcs = (N_DEL_PROCESSES - 1) *
            (DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1 + DCS_SAMPLING_N) *
            (settings.n_energies - settings.dcs_model_offset);
        size_data[imem++] = memory_padded_size(sizeof(struct atomic_element *) *
                                    settings.n_elements,
                                pad_size) +
            (sizeof(struct atomic_element) +
                memory_padded_size(n_dcs * sizeof(float), pad_size)) *
                settings.n_elements +
            settings.size_elements_names;
        /* Atomic composition. */
        size_data[imem++] =
            memory_padded_size(
                sizeof(struct material_component *) * settings.n_materials,
                pad_size) +
            memory_padded_size(
                settings.n_components * sizeof(struct material_component),
                pad_size);
        /* Composite material. */
        size_data[imem++] =
            memory_padded_size(
                sizeof(struct composite_material *) * settings.n_composites,
                pad_size) +
            settings.size_composite;
        /* material_name. */
        size_data[imem++] =
            memory_padded_size(sizeof(char *) * settings.n_materials +
                    settings.size_materials_names,
                pad_size);
        /* Bremsstrahlung model name. */
        size_data[imem++] = memory_padded_size(
            sizeof(char) * (strlen(opts.bremsstrahlung) + 1), pad_size);
        /* Pair production model name. */
        size_data[imem++] = memory_padded_size(
            sizeof(char) * (strlen(opts.pair_production) + 1), pad_size);
        /* Photonuclear model name. */
        size_data[imem++] = memory_padded_size(
            sizeof(char) * (strlen(opts.photonuclear) + 1), pad_size);

        /* Allocate the shared memory. */
        int size_total = 0;
        for (imem = 0; imem < N_DATA_POINTERS; imem++)
                size_total += size_data[imem];
        const int size_shared = sizeof(**physics_ptr) + size_total;
        void * tmp_ptr = reallocate(mdf, size_shared);
        if (tmp_ptr == NULL) {
                ERROR_REGISTER_MEMORY();
                goto clean_and_exit;
        }
        physics = tmp_ptr;
        *physics_ptr = physics;
        memset(physics, 0x0, size_shared);
        physics->size = size_shared;
        mdf = NULL;

        /* Map the data pointers. */
        double * p = physics->data;
        void ** ptr = (void **)(&(physics->mdf_path));
        for (imem = 0; imem < N_DATA_POINTERS; imem++) {
                *ptr = p;
                ptr++;
                p += size_data[imem] / pad_size;
        }

        /* Set the DCS's. */
        pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG, opts.bremsstrahlung,
            &physics->dcs_bremsstrahlung);
        strcpy(physics->model_bremsstrahlung, opts.bremsstrahlung);
        pumas_dcs_get(PUMAS_PROCESS_PAIR_PRODUCTION, opts.pair_production,
            &physics->dcs_pair_production);
        strcpy(physics->model_pair_production, opts.pair_production);
        pumas_dcs_get(PUMAS_PROCESS_PHOTONUCLEAR, opts.photonuclear,
            &physics->dcs_photonuclear);
        strcpy(physics->model_photonuclear, opts.photonuclear);

        /* Copy the global settings. */
        physics->particle = particle;
        if (particle == PUMAS_PARTICLE_MUON) {
                physics->ctau = MUON_C_TAU;
                physics->mass = MUON_MASS;
        } else {
                physics->ctau = TAU_C_TAU;
                physics->mass = TAU_MASS;
        }
        physics->n_energies = settings.n_energies;
        physics->n_materials = settings.n_materials;
        physics->n_composites = settings.n_composites;
        physics->n_elements = settings.n_elements;
        physics->n_components = settings.n_components;
        physics->max_components = settings.max_components;
        physics->n_energy_loss_header = settings.n_energy_loss_header;
        physics->dcs_model_offset = settings.dcs_model_offset;
        strcpy(physics->mdf_path, file_mdf);

        /* Set the cutoff */
        physics->cutoff = opts.cutoff;

        /* Allocate a new MDF buffer. */
        if ((mdf = allocate(sizeof(struct mdf_buffer) + size_mdf)) == NULL) {
                ERROR_REGISTER_MEMORY();
                goto clean_and_exit;
        }
        memcpy(mdf, &settings, sizeof(settings));

        /* Set the path to the dE/dX files. */
        if (!dry_mode)
                strcpy(physics->dedx_path, dedx_path);
        else
                physics->dedx_path = NULL;

        /* Parse the elements. */
        if ((mdf_parse_elements(physics, mdf, error_)) != PUMAS_RETURN_SUCCESS)
                goto clean_and_exit;

        /* Parse the base materials. */
        if ((mdf_parse_materials(physics, mdf, error_)) != PUMAS_RETURN_SUCCESS)
                goto clean_and_exit;

        /* Parse the composite materials. */
        if ((mdf_parse_composites(physics, mdf, error_)) !=
            PUMAS_RETURN_SUCCESS)
                goto clean_and_exit;

        /* All done if in dry mode. */
        if (dry_mode) goto clean_and_exit;

        /* Precompute the CEL integrals and the TT parameters. */
        int imat;
        for (imat = 0; imat < physics->n_materials - physics->n_composites;
             imat++) {
                int ikin;
                for (ikin = 0; ikin < physics->n_energies; ikin++) {
                        if (compute_coulomb_parameters(physics, imat, ikin,
                                error_) != PUMAS_RETURN_SUCCESS)
                                goto clean_and_exit;
                }
                compute_cel_integrals(physics, imat);
                compute_csda_magnetic_transport(physics, imat);
        }

        /* Precompute the same properties for composite materials. */
        for (imat = physics->n_materials - physics->n_composites;
             imat < physics->n_materials; imat++) {
                if (((compute_composite(physics, imat, error_)) !=
                        PUMAS_RETURN_SUCCESS) ||
                    ((compute_composite_density(physics, imat, error_)) !=
                        PUMAS_RETURN_SUCCESS))
                        goto clean_and_exit;
        }

        /* Tabulate and fit the DCS for atomic elements. */
        int iel;
        for (iel = 0; iel < physics->n_elements; iel++) {
                int ip;
                for (ip = 0; ip < N_DEL_PROCESSES - 1; ip++) {
                        if (compute_dcs_model(physics, dcs_get(ip),
                                physics->element[iel],
                                error_) != PUMAS_RETURN_SUCCESS)
                                goto clean_and_exit;
                }
        }

clean_and_exit:
        if (fid_mdf != NULL) fclose(fid_mdf);
        deallocate(mdf);
        io_read_line(NULL, NULL, NULL, 0, error_);
        compute_coulomb_parameters(physics, -1, -1, error_);
        compute_coulomb_soft(physics, -1, NULL, error_);
        compute_cel_and_del(physics, -1);
        compute_dcs_model(physics, NULL, NULL, error_);
        if ((error_->code != PUMAS_RETURN_SUCCESS) && (physics != NULL)) {
                deallocate(physics);
                *physics_ptr = NULL;
        }

        return ERROR_RAISE();
}

/* The standard API initialisation. */
enum pumas_return pumas_physics_create(struct pumas_physics ** physics,
    enum pumas_particle particle, const char * mdf_path, const char * dedx_path,
    const struct pumas_physics_settings * settings)
{
        return _initialise(physics, particle, mdf_path, dedx_path, 0, settings);
}

enum pumas_return pumas_physics_load(
    struct pumas_physics ** physics_ptr, FILE * stream)
{
        ERROR_INITIALISE(pumas_physics_load);

        /* Check the physics pointer. */
        if (physics_ptr == NULL) {
                return ERROR_NULL_PHYSICS();
        }

        /* Check the input stream */
        if (stream == NULL)
                return ERROR_MESSAGE(
                    PUMAS_RETURN_PATH_ERROR, "invalid input stream (null)");
#if (GDB_MODE)
        /* Save the floating points exceptions status and enable them. */
        fe_status = fegetexcept();
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
        /* Check the binary dump tag. */
        struct pumas_physics * physics = NULL;
        int tag;
        if (fread(&tag, sizeof(tag), 1, stream) != 1) goto error;
        if (tag != PHYSICS_BINARY_DUMP_TAG) {
                ERROR_REGISTER(PUMAS_RETURN_FORMAT_ERROR,
                    "incompatible version of binary dump");
                goto error;
        }

        /* Allocate the container. */
        int size;
        if (fread(&size, sizeof(size), 1, stream) != 1) goto error;
        physics = allocate(size);
        *physics_ptr = physics;
        if (physics == NULL) {
                ERROR_REGISTER_MEMORY();
                goto error;
        }

        /* Load the data and remap the addresses. */
        if (fread(physics, size, 1, stream) != 1) goto error;

        void ** ptr = (void **)(&(physics->mdf_path));
        ptrdiff_t delta = (char *)(physics->data) - (char *)(*ptr);
        int i;
        for (i = 0; i < N_DATA_POINTERS; i++, ptr++)
                *ptr = ((char *)(*ptr)) + delta;

        struct atomic_element ** element = physics->element;
        for (i = 0; i < physics->n_elements; i++) {
                element[i] =
                    (struct atomic_element *)(((char *)element[i]) + delta);
                element[i]->name += delta;
                element[i]->dcs_data =
                    (float *)(((char *)(element[i]->dcs_data)) + delta);
        }

        struct material_component ** composition = physics->composition;
        for (i = 0; i < physics->n_materials; i++)
                composition[i] =
                    (struct material_component *)(((char *)composition[i]) +
                        delta);

        struct composite_material ** composite = physics->composite;
        for (i = 0; i < physics->n_composites; i++)
                composite[i] =
                    (struct composite_material *)(((char *)composite[i]) +
                        delta);

        char ** material_name = physics->material_name;
        for (i = 0; i < physics->n_materials; i++) material_name[i] += delta;

        /* Set the DCS models */
        if (dcs_check_model(PUMAS_PROCESS_BREMSSTRAHLUNG,
            physics->model_bremsstrahlung, error_) == PUMAS_RETURN_SUCCESS) {
                pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG,
                    physics->model_bremsstrahlung,
                    &physics->dcs_bremsstrahlung);
        } else {
                goto error;
        }

        if (dcs_check_model(PUMAS_PROCESS_PAIR_PRODUCTION,
            physics->model_pair_production, error_) == PUMAS_RETURN_SUCCESS) {
                pumas_dcs_get(PUMAS_PROCESS_PAIR_PRODUCTION,
                    physics->model_pair_production,
                    &physics->dcs_pair_production);
        } else {
                goto error;
        }

        if (dcs_check_model(PUMAS_PROCESS_PHOTONUCLEAR,
            physics->model_photonuclear, error_) == PUMAS_RETURN_SUCCESS) {
                pumas_dcs_get(PUMAS_PROCESS_PHOTONUCLEAR,
                    physics->model_photonuclear,
                    &physics->dcs_photonuclear);
        } else {
                goto error;
        }

        return PUMAS_RETURN_SUCCESS;

error:
        deallocate(physics);
        *physics_ptr = NULL;
        return ERROR_RAISE();

#undef N_DATA_POINTERS
}

enum pumas_return pumas_physics_dump(
    const struct pumas_physics * physics, FILE * stream)
{
        ERROR_INITIALISE(pumas_physics_dump);

        /* Check if the Physics is initialised. */
        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }

        /* Check the output stream */
        if (stream == NULL)
                return ERROR_MESSAGE(
                    PUMAS_RETURN_PATH_ERROR, "invalid output stream (null)");

        /* Dump the configuration. */
        int tag = PHYSICS_BINARY_DUMP_TAG;
        if (fwrite(&tag, sizeof(tag), 1, stream) != 1) goto error;
        if (fwrite(&physics->size, sizeof(physics->size), 1, stream) != 1)
                goto error;
        if (fwrite(physics, physics->size, 1, stream) != 1) goto error;

        return PUMAS_RETURN_SUCCESS;

error:
        return ERROR_MESSAGE(
            PUMAS_RETURN_IO_ERROR, "could not write to dump file");

#undef PHYSICS_BINARY_DUMP_TAG
}

void pumas_physics_destroy(struct pumas_physics ** physics_ptr)
{
        if ((physics_ptr == NULL) || (*physics_ptr == NULL)) return;

        /* Free the shared data. */
        struct pumas_physics * physics = *physics_ptr;
        int i;
        for (i = 0; i < physics->n_materials - physics->n_composites; i++) {
                deallocate(physics->dedx_filename[i]);
                physics->dedx_filename[i] = NULL;
        }
        deallocate(physics);
        *physics_ptr = NULL;
}

double pumas_physics_cutoff(const struct pumas_physics * physics)
{
        return (physics == NULL) ? -1 : physics->cutoff;
}

const char * pumas_error_function(pumas_function_t * caller)
{
#define TOSTRING(function)                                                     \
        if (caller == (pumas_function_t *)function) return #function;

        /* Library functions with an error code. */
        TOSTRING(pumas_physics_create)
        TOSTRING(pumas_physics_create_tabulation)
        TOSTRING(pumas_physics_dump)
        TOSTRING(pumas_physics_load)
        TOSTRING(pumas_physics_tabulate)
        TOSTRING(pumas_context_transport)
        TOSTRING(pumas_physics_particle)
        TOSTRING(pumas_context_create)
        TOSTRING(pumas_context_random_dump)
        TOSTRING(pumas_context_random_load)
        TOSTRING(pumas_context_random_seed_get)
        TOSTRING(pumas_context_random_seed_set)
        TOSTRING(pumas_recorder_create)
        TOSTRING(pumas_physics_dcs)
        TOSTRING(pumas_physics_element_name)
        TOSTRING(pumas_physics_element_index)
        TOSTRING(pumas_physics_element_properties)
        TOSTRING(pumas_physics_material_name)
        TOSTRING(pumas_physics_material_index)
        TOSTRING(pumas_physics_material_properties)
        TOSTRING(pumas_physics_composite_update)
        TOSTRING(pumas_physics_composite_properties)
        TOSTRING(pumas_physics_print)
        TOSTRING(pumas_error_raise)
        TOSTRING(pumas_physics_property_grammage)
        TOSTRING(pumas_physics_property_proper_time)
        TOSTRING(pumas_physics_property_magnetic_rotation)
        TOSTRING(pumas_physics_property_kinetic_energy)
        TOSTRING(pumas_physics_property_energy_loss)
        TOSTRING(pumas_physics_property_scattering_length)
        TOSTRING(pumas_physics_property_cross_section)
        TOSTRING(pumas_physics_table_value)
        TOSTRING(pumas_physics_table_index)
        TOSTRING(pumas_dcs_get)
        TOSTRING(pumas_dcs_register)

        /* Other library functions. */
        TOSTRING(pumas_constant)
        TOSTRING(pumas_dcs_default)
        TOSTRING(pumas_physics_cutoff)
        TOSTRING(pumas_physics_destroy)
        TOSTRING(pumas_physics_tabulation_clear)
        TOSTRING(pumas_context_destroy)
        TOSTRING(pumas_context_physics_get)
        TOSTRING(pumas_recorder_clear)
        TOSTRING(pumas_recorder_destroy)
        TOSTRING(pumas_version)
        TOSTRING(pumas_error_function)
        TOSTRING(pumas_error_handler_set)
        TOSTRING(pumas_error_handler_get)
        TOSTRING(pumas_error_catch)
        TOSTRING(pumas_physics_element_length)
        TOSTRING(pumas_physics_material_length)
        TOSTRING(pumas_physics_composite_length)
        TOSTRING(pumas_physics_table_length)
        TOSTRING(pumas_memory_allocator)
        TOSTRING(pumas_memory_reallocator)
        TOSTRING(pumas_memory_deallocator)

        return NULL;
#undef TOSTRING
}

void pumas_error_handler_set(pumas_handler_cb * handler)
{
        s_error.handler = handler;
}

pumas_handler_cb * pumas_error_handler_get(void) { return s_error.handler; }

void pumas_error_catch(int catch)
{
        if (catch) {
                s_error.catch = 1;
                s_error.catch_error.code = PUMAS_RETURN_SUCCESS;
                s_error.catch_error.function = NULL;
        } else {
                s_error.catch = 0;
        }
}

enum pumas_return pumas_error_raise(void)
{
        ERROR_INITIALISE(pumas_error_raise);

        if (s_error.catch == 0)
                ERROR_MESSAGE(
                    PUMAS_RETURN_RAISE_ERROR, "`raise' called without `catch'");
        s_error.catch = 0;
        memcpy(error_, &s_error.catch_error, sizeof(*error_));
        return ERROR_RAISE();
}

/* Set the MT initial state */
static enum pumas_return random_initialise(struct pumas_context * context,
    const unsigned long * seed_ptr, struct error_context * error_)
{
        unsigned long seed;
        if (seed_ptr == NULL) {
                /* Sample the seed from the OS */
#ifdef _WIN32
                unsigned int tmp;
                if (rand_s(&tmp) != 0) {
                        return ERROR_REGISTER(PUMAS_RETURN_PATH_ERROR,
                                "could not read from `rand_s'");
                } else {
                        seed = tmp;
                }
#else
                size_t count = 0;
                FILE * fid = fopen("/dev/urandom", "r");
                if (fid != NULL) {
                        count = fread(&seed, sizeof seed, 1, fid);
                        fclose(fid);
                }

                if (count == 0) {
                        return ERROR_REGISTER(PUMAS_RETURN_PATH_ERROR,
                                "could not read from `/dev/urandom'");
                }
#endif
        } else {
                seed = *seed_ptr;
        }

        struct simulation_context * context_ = (void *)context;
        if (context_->random_data == NULL) {
                context_->random_data =
                    allocate(sizeof(*context_->random_data));
                if (context_->random_data == NULL) {
                        return ERROR_REGISTER_MEMORY();
                }
        }
        struct pumas_random_data * data = context_->random_data;

        memset(data, 0x0, sizeof (*data));
        data->buffer[0] = seed & 0xffffffffUL;
        int j;
        for (j = 1; j < MT_PERIOD; j++) {
                data->buffer[j] = (1812433253UL *
                        (data->buffer[j - 1] ^
                        (data->buffer[j - 1] >> 30)) +
                    j);
                data->buffer[j] &= 0xffffffffUL;
        }
        data->seed = seed;
        data->index = MT_PERIOD;

        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_context_random_seed_set(
    struct pumas_context * context, const unsigned long * seed_ptr)
{
        ERROR_INITIALISE(pumas_context_random_seed_set);

        random_initialise(context, seed_ptr, error_);

        return ERROR_RAISE();
}

enum pumas_return pumas_context_random_seed_get(
    struct pumas_context * context, unsigned long * seed_ptr)
{
        ERROR_INITIALISE(pumas_context_random_seed_get);

        struct simulation_context * context_ = (void *)context;
        if (context_->random_data == NULL) {
                enum pumas_return rc = random_initialise(context, NULL, error_);
                if (rc != PUMAS_RETURN_SUCCESS) {
                        *seed_ptr = 0;
                        return ERROR_RAISE();
                }
        }
        *seed_ptr = context_->random_data->seed;

        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_context_random_dump(
    struct pumas_context * context, FILE * stream)
{
        ERROR_INITIALISE(pumas_context_random_dump);

        /* Check the input stream */
        if (stream == NULL)
                return ERROR_MESSAGE(
                    PUMAS_RETURN_PATH_ERROR, "invalid input stream (null)");

        /* Write the version tag */
        const int tag = RANDOM_BINARY_DUMP_TAG;
        if (fwrite(&tag, sizeof(tag), 1, stream) != 1) goto error;

        /* Initialise the random engine if needed */
        struct simulation_context * context_ = (void *)context;
        if (context_->random_data == NULL) {
                enum pumas_return rc = random_initialise(context, NULL, error_);
                if (rc != PUMAS_RETURN_SUCCESS)
                        return ERROR_RAISE();
        }

        /* Dump the data */
        if (fwrite(context_->random_data, sizeof (*context_->random_data), 1,
            stream) != 1) goto error;

        return PUMAS_RETURN_SUCCESS;

error:
        return ERROR_MESSAGE(
            PUMAS_RETURN_IO_ERROR, "could not not write to stream");
}

enum pumas_return pumas_context_random_load(
    struct pumas_context * context, FILE * stream)
{
        ERROR_INITIALISE(pumas_context_random_load);

        /* Check the input stream */
        if (stream == NULL)
                return ERROR_MESSAGE(
                    PUMAS_RETURN_PATH_ERROR, "invalid input stream (null)");

        /* Check the binary dump tag */
        int tag;
        if (fread(&tag, sizeof(tag), 1, stream) != 1) goto error;
        if (tag != RANDOM_BINARY_DUMP_TAG) {
                ERROR_REGISTER(PUMAS_RETURN_FORMAT_ERROR,
                    "incompatible version of binary dump");
        }

        /* Allocate memory if needed */
        struct simulation_context * context_ = (void *)context;
        if (context_->random_data == NULL) {
                context_->random_data =
                    allocate(sizeof(*context_->random_data));
                if (context_->random_data == NULL) {
                        return ERROR_MESSAGE(PUMAS_RETURN_MEMORY_ERROR,
                            "could not allocate memory");
                }
        }

        /* Load the data */
        if (fread(context_->random_data, sizeof (*context_->random_data), 1,
            stream) != 1) goto error;

        return PUMAS_RETURN_SUCCESS;

error:
        return ERROR_MESSAGE(
            PUMAS_RETURN_IO_ERROR, "could not not read from stream");

#undef RANDOM_BINARY_DUMP_TAG
}

/* Uniform pseudo random distribution from a Mersenne Twister */
static double random_uniform01(struct pumas_context * context)
{
        /* Lazy initialisation of the MT if not already done */
        struct simulation_context * context_ = (void *)context;
        if (context_->random_data == NULL) {
                if (random_initialise(context, NULL, NULL) !=
                    PUMAS_RETURN_SUCCESS) return -1.;

        }
        struct pumas_random_data * data = context_->random_data;

        /* Check the buffer */
        if (data->index < MT_PERIOD - 1) {
                data->index++;
        } else {
                /* Update the MT state */
                const int M = 397;
                const unsigned long UPPER_MASK = 0x80000000UL;
                const unsigned long LOWER_MASK = 0x7fffffffUL;
                static unsigned long mag01[2] = { 0x0UL, 0x9908b0dfUL };
                unsigned long y;
                int kk;
                for (kk = 0; kk < MT_PERIOD - M; kk++) {
                        y = (data->buffer[kk] & UPPER_MASK) |
                            (data->buffer[kk + 1] & LOWER_MASK);
                        data->buffer[kk] = data->buffer[kk + M] ^
                            (y >> 1) ^ mag01[y & 0x1UL];
                }
                for (; kk < MT_PERIOD - 1; kk++) {
                        y = (data->buffer[kk] & UPPER_MASK) |
                            (data->buffer[kk + 1] & LOWER_MASK);
                        data->buffer[kk] =
                            data->buffer[kk + (M - MT_PERIOD)] ^
                            (y >> 1) ^ mag01[y & 0x1UL];
                }
                y = (data->buffer[MT_PERIOD - 1] & UPPER_MASK) |
                    (data->buffer[0] & LOWER_MASK);
                data->buffer[MT_PERIOD - 1] =
                    data->buffer[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
                data->index = 0;
        }

        /* Tempering */
        unsigned long y = data->buffer[data->index];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);

        /* Convert to a floating point and return */
        return y * (1.0 / 4294967295.0);
}

/* Public library functions: simulation context management. */
enum pumas_return pumas_context_create(struct pumas_context ** context_,
    const struct pumas_physics * physics, int extra_memory)
{
        ERROR_INITIALISE(pumas_context_create);
        *context_ = NULL;

        /* Check the Physics initialisation. */
        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }

        /* Allocate the new context. */
        struct simulation_context * context;
        const int pad_size = sizeof(*(context->data));
        const int work_size =
            memory_padded_size(sizeof(struct coulomb_workspace) +
                    physics->max_components * sizeof(struct coulomb_data),
                pad_size);
        if (extra_memory < 0)
                extra_memory = 0;
        else
                extra_memory = memory_padded_size(extra_memory, pad_size);
        context = allocate(sizeof(*context) + work_size + extra_memory);
        if (context == NULL) {
                ERROR_REGISTER_MEMORY();
                return ERROR_RAISE();
        }

        /* Set the default configuration. */
        *context_ = (struct pumas_context *)context;
        context->physics = physics;
        context->extra_memory = extra_memory;
        if (extra_memory > 0)
                (*context_)->user_data = context->data + work_size / pad_size;
        else
                (*context_)->user_data = NULL;

        int imax = physics->n_energies - 2;
        context->index_K_last[0] = context->index_K_last[1] = imax;
        context->index_X_last[0] = context->index_X_last[1] = imax;

        context->random_data = NULL;
        (*context_)->random = &random_uniform01;

        (*context_)->medium = NULL;
        (*context_)->recorder = NULL;

        (*context_)->mode.decay = (physics->particle == PUMAS_PARTICLE_MUON) ?
            PUMAS_MODE_WEIGHT :
            PUMAS_MODE_DECAY;
        (*context_)->mode.direction = PUMAS_MODE_FORWARD;
        (*context_)->mode.energy_loss = PUMAS_MODE_DETAILED;
        (*context_)->mode.scattering = PUMAS_MODE_FULL_SPACE;
        (*context_)->event = PUMAS_EVENT_NONE;

        (*context_)->limit.energy = 0.;    /* GeV */
        (*context_)->limit.distance = 0.;  /* m */
        (*context_)->limit.grammage = 0.;  /* kg/m^2 */
        (*context_)->limit.time = 0.;      /* m/c */

        (*context_)->accuracy = DEFAULT_ACCURACY;

        /* Initialise the Gaussian transform of the random stream. */
        context->randn_done = 0;
        context->randn_next = 0.;

        /* Initialise the work space. */
        context->workspace = (struct coulomb_workspace *)context->data;

        return PUMAS_RETURN_SUCCESS;
}

void pumas_context_destroy(struct pumas_context ** context)
{
        /* Check that the context hasn't already been destroyed */
        if ((context == NULL) || (*context == NULL)) return;

        /* Release the memory */
        struct simulation_context * context_ = (void *)(*context);
        deallocate(context_->random_data);
        deallocate(*context);
        *context = NULL;
}

const struct pumas_physics * pumas_context_physics_get(
    const struct pumas_context * context)
{
        if (context == NULL) return NULL;

        const struct simulation_context * context_ =
            (const struct simulation_context *)context;
        return context_->physics;
}

/* Public library functions: global print routines */
enum pumas_return pumas_physics_print(const struct pumas_physics * physics,
    FILE * stream, const char * tabulation, const char * newline)
{
        ERROR_INITIALISE(pumas_physics_print);

        const char * tab = (tabulation == NULL) ? "" : tabulation;
        const char * cr = (newline == NULL) ? "" : newline;

        /* Check the Physics initialisation */
        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }

        /* Check the output stream */
        if (stream == NULL) goto error;

        /* Print the particle info */
        if (fprintf(stream, "{%s%s\"particle\" : {%s%s%s\"mass (GeV/c^2)\""
                            " : %.6lf",
                cr, tab, cr, tab, tab, physics->mass) < 0)
                goto error;
        if (fprintf(stream, ",%s%s%s\"lifetime (m/c)\" : %.3lf%s%s}", cr, tab,
                tab, physics->ctau, cr, tab) < 0)
                goto error;

        /* Print the atomic elements */
        if (fprintf(stream, ",%s%s\"elements\" : {", cr, tab) < 0) goto error;
        int iel = 0;
        for (; iel < physics->n_elements; iel++) {
                const char * head = (iel == 0) ? "" : ",";
                const struct atomic_element * element = physics->element[iel];
                if (fprintf(stream, "%s%s%s%s\"%s\" : {", head, cr, tab, tab,
                        element->name) < 0)
                        goto error;
                if (fprintf(stream, "%s%s%s%s\"Z\" : %.0lf", cr, tab, tab, tab,
                        element->Z) < 0)
                        goto error;
                if (fprintf(stream, ",%s%s%s%s\"A (g/mol)\" : %.5lg", cr, tab,
                        tab, tab, element->A) < 0)
                        goto error;
                if (fprintf(stream, ",%s%s%s%s\"I (eV)\" : %.1lf%s%s%s}", cr,
                        tab, tab, tab, element->I * 1E+09, cr, tab, tab) < 0)
                        goto error;
        }

        /* Print the materials */
        if (fprintf(stream, ",%s%s\"materials\" : {", cr, tab) < 0) goto error;
        int material = 0;
        for (; material < physics->n_materials - physics->n_composites;
             material++) {
                const char * head = (material == 0) ? "" : ",";
                if (fprintf(stream, "%s%s%s%s\"%s\" : {", head, cr, tab, tab,
                        physics->material_name[material]) < 0)
                        goto error;
                if (fprintf(stream, "%s%s%s%s\"density\" : %.5lg", cr, tab, tab,
                        tab, physics->material_density[material] * 1E-03) < 0)
                        goto error;
                if (fprintf(stream, ",%s%s%s%s\"elements\" : {", cr, tab, tab,
                        tab) < 0)
                        goto error;
                int iel = 0;
                for (; iel < physics->elements_in[material]; iel++) {
                        const char * head2 = (iel == 0) ? "" : ",";
                        int element =
                            physics->composition[material][iel].element;
                        if (fprintf(stream, "%s%s%s%s%s%s\"%s (%%)\" : %.5lg",
                                head2, cr, tab, tab, tab, tab,
                                physics->element[element]->name, 100. *
                                    physics->composition[material][iel]
                                        .fraction) < 0)
                                goto error;
                }
                if (fprintf(stream, "%s%s%s%s}%s%s%s}",
                    cr, tab, tab, tab, cr, tab, tab) < 0)
                        goto error;
        }
        if (fprintf(stream, "%s%s}", cr, tab) < 0) goto error;
        if (physics->n_composites <= 0) goto closure;

        /* Print the composites */
        if (fprintf(stream, ",%s%s\"composites\" : {", cr, tab) < 0) goto error;
        const int material0 = physics->n_materials - physics->n_composites;
        material = material0;
        for (; material < physics->n_materials; material++) {
                const char * head = (material == material0) ? "" : ",";
                struct composite_material * composite =
                    physics->composite[material - material0];
                if (fprintf(stream, "%s%s%s%s\"%s\" : {", head, cr, tab, tab,
                        physics->material_name[material]) < 0)
                        goto error;
                if (fprintf(stream, "%s%s%s%s\"density\" : %.5lg", cr, tab, tab,
                        tab, physics->material_density[material] * 1E-03) < 0)
                        goto error;
                if (fprintf(stream, ",%s%s%s%s\"materials\" : {", cr, tab, tab,
                        tab) < 0)
                        goto error;

                int imat = 0;
                for (; imat < composite->n_components; imat++) {
                        const char * head2 = (imat == 0) ? "" : ",";
                        struct composite_component * c =
                            composite->component + imat;
                        if (fprintf(stream, "%s%s%s%s%s%s\"%s\" : {", head2, cr,
                                tab, tab, tab, tab,
                                physics->material_name[c->material]) < 0)
                                goto error;
                        if (fprintf(stream, ",%s%s%s%s%s%s\"fraction (%%)\" "
                                            ": %.5lg%s%s%s%s%s}",
                                cr, tab, tab, tab, tab, tab, 100. * c->fraction,
                                cr, tab, tab, tab, tab) < 0)
                                goto error;
                }
                if (fprintf(stream, "%s%s%s%s}%s%s%s}", cr, tab, tab, tab, cr,
                        tab, tab) < 0)
                        goto error;
        }
        if (fprintf(stream, "%s%s}", cr, tab) < 0) goto error;

closure:
        if (fprintf(stream, "%s}", cr) < 0) goto error;

        return PUMAS_RETURN_SUCCESS;
error:
        return ERROR_MESSAGE(
            PUMAS_RETURN_IO_ERROR, "could not write to stream");
}

int pumas_version() { return 100 * PUMAS_VERSION; }

/* Public library functions: recorder handling. */
enum pumas_return pumas_recorder_create(
    struct pumas_recorder ** recorder_, int extra_memory)
{
        ERROR_INITIALISE(pumas_recorder_create);
        *recorder_ = NULL;

        /*  Allocate memory for the new recorder. */
        struct frame_recorder * recorder = NULL;
        if (extra_memory < 0) extra_memory = 0;
        recorder = allocate(sizeof(*recorder) + extra_memory);
        if (recorder == NULL) {
                ERROR_REGISTER_MEMORY();
                return ERROR_RAISE();
        }

        /* Configure the context in order to use the new recorder. */
        *recorder_ = (struct pumas_recorder *)recorder;

        /* configure the new recorder and return it. */
        (*recorder_)->period = 1;
        (*recorder_)->record = NULL;
        (*recorder_)->length = 0;
        (*recorder_)->first = NULL;
        (*recorder_)->user_data = (extra_memory > 0) ? recorder->data : NULL;
        recorder->last = NULL;
        recorder->stack = NULL;

        return PUMAS_RETURN_SUCCESS;
}

void pumas_recorder_destroy(struct pumas_recorder ** recorder)
{
        if ((recorder == NULL) || (*recorder == NULL)) return;
        pumas_recorder_clear(*recorder);
        deallocate(*recorder);
        *recorder = NULL;
}

void pumas_recorder_clear(struct pumas_recorder * recorder)
{
        if (recorder == NULL) return;

        struct frame_recorder * const rec =
            (struct frame_recorder * const)recorder;
        struct frame_stack * current = rec->stack;
        while (current != NULL) {
                struct frame_stack * next = current->next;
                deallocate(current);
                current = next;
        }
        rec->stack = NULL;
        rec->last = NULL;
        recorder->first = NULL;
        recorder->length = 0;
}

/* Public library functions: properties accessors. */
enum pumas_return pumas_physics_property_grammage(
    const struct pumas_physics * physics, enum pumas_mode scheme,
    int material, double kinetic, double * grammage)
{
        ERROR_INITIALISE(pumas_physics_property_grammage);
        *grammage = 0.;

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        } else if ((scheme <= PUMAS_MODE_VIRTUAL) ||
            (scheme >= PUMAS_MODE_DETAILED)) {
                return ERROR_INVALID_SCHEME(scheme);
        } else if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        *grammage = cel_grammage(physics, NULL, scheme, material, kinetic);
        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_physics_property_proper_time(
    const struct pumas_physics * physics, enum pumas_mode scheme,
    int material, double kinetic, double * time)
{
        ERROR_INITIALISE(pumas_physics_property_proper_time);
        *time = 0.;

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        } else if ((scheme <= PUMAS_MODE_VIRTUAL) ||
            (scheme >= PUMAS_MODE_DETAILED)) {
                return ERROR_INVALID_SCHEME(scheme);
        } else if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        *time = cel_proper_time(physics, NULL, scheme, material, kinetic);
        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_physics_property_magnetic_rotation(
    const struct pumas_physics * physics, int material, double kinetic,
    double * angle)
{
        ERROR_INITIALISE(pumas_physics_property_magnetic_rotation);
        *angle = 0.;

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        } else if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        *angle = cel_magnetic_rotation(physics, NULL, material, kinetic);
        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_physics_property_kinetic_energy(
    const struct pumas_physics * physics, enum pumas_mode scheme,
    int material, double grammage, double * kinetic)
{
        ERROR_INITIALISE(pumas_physics_property_kinetic_energy);
        *kinetic = 0.;

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        } else if ((scheme <= PUMAS_MODE_VIRTUAL) ||
            (scheme >= PUMAS_MODE_DETAILED)) {
                return ERROR_INVALID_SCHEME(scheme);
        } else if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        *kinetic =
            cel_kinetic_energy(physics, NULL, scheme, material, grammage);
        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_physics_property_energy_loss(
    const struct pumas_physics * physics, enum pumas_mode scheme,
    int material, double kinetic, double * dedx)
{
        ERROR_INITIALISE(pumas_physics_property_energy_loss);
        *dedx = 0.;

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        } else if ((scheme <= PUMAS_MODE_VIRTUAL) ||
            (scheme >= PUMAS_MODE_DETAILED)) {
                return ERROR_INVALID_SCHEME(scheme);
        } else if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        *dedx = cel_energy_loss(physics, NULL, scheme, material, kinetic);
        return PUMAS_RETURN_SUCCESS;
}

/* Public library function: elastic scattering 1st transport path length. */
enum pumas_return pumas_physics_property_scattering_length(
    const struct pumas_physics * physics, int material, double kinetic,
    double * length)
{
        ERROR_INITIALISE(pumas_physics_property_scattering_length);
        *length = 0.;

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        } else if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        double path = 0.;
        int i = 0;
        for (; i < physics->elements_in[material]; i++) {
                const struct material_component * component =
                    &physics->composition[material][i];
                const struct atomic_element * element =
                    physics->element[component->element];
                double kinetic0, screening[3], coefficient[2], fCM[2];
                coulomb_frame_parameters(
                    physics, kinetic, element->A, &kinetic0, fCM);
                coulomb_screening_parameters(
                    physics, NULL, kinetic0, component->element, screening);
                const double fspin = coulomb_spin_factor(physics, kinetic0);
                double a[3], b[3];
                coulomb_pole_decomposition(screening, a, b);
                coulomb_transport_coefficients(
                    1., fspin, screening, a, b, coefficient);
                double d = 1. / (fCM[0] * (1. + fCM[1]));
                d *= d;
                coefficient[1] *= d;
                path += component->fraction /
                    coulomb_wentzel_path(physics, kinetic, element->Z,
                            element->A, screening[0]) *
                    coefficient[1];
        }

        if (path == 0.) {
                *length = DBL_MAX;
                return ERROR_MESSAGE(
                    PUMAS_RETURN_VALUE_ERROR, "infinite scattering length");
        }

        *length = 1. / path;
        return PUMAS_RETURN_SUCCESS;
}

/* Public library function: total inelastic cross section. */
enum pumas_return pumas_physics_property_cross_section(
    const struct pumas_physics * physics, int material, double kinetic,
    double * cross_section)
{
        ERROR_INITIALISE(pumas_physics_property_cross_section);
        *cross_section = 0.;

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        } else if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        *cross_section = (kinetic < *table_get_Kt(physics, material)) ?
            0. :
            del_cross_section(physics, NULL, material, kinetic);
        return PUMAS_RETURN_SUCCESS;
}

/* Public library function: the main transport routine. */
enum pumas_return pumas_context_transport(struct pumas_context * context,
    struct pumas_state * state, enum pumas_event * event,
    struct pumas_medium * media[2])
{
        ERROR_INITIALISE(pumas_context_transport);

        /* Check the context and state */
        if (context == NULL)
                return ERROR_MESSAGE(
                    PUMAS_RETURN_VALUE_ERROR, "no context (null)");
        if (state == NULL)
                return ERROR_MESSAGE(
                    PUMAS_RETURN_VALUE_ERROR, "no state (null)");

        /* Check the Physics initialisation */
        struct simulation_context * context_ =
            (struct simulation_context *)context;
        const struct pumas_physics * physics = context_->physics;
        if (physics == NULL) return ERROR_NOT_INITIALISED();

        /* Check the initial state. */
        if (state->decayed) {
                if (event != NULL) *event = PUMAS_EVENT_VERTEX_DECAY;
                return PUMAS_RETURN_SUCCESS;
        }

        /* Check the direction norm */
        {
                const double * const u = state->direction;
                const double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
                if (fabs(norm2 - 1) > FLT_EPSILON) {
                        return ERROR_FORMAT(
                            PUMAS_RETURN_DIRECTION_ERROR,
                            "bad norm for state direction (norm^2 - 1 = %g)",
                            norm2 - 1);
                }
        }

        /* Check the configuration. */
        if (event != NULL) *event = PUMAS_EVENT_NONE;
        if (context->medium == NULL) {
                return ERROR_MESSAGE(
                    PUMAS_RETURN_MEDIUM_ERROR, "no medium specified");
        } else if ((physics->particle == PUMAS_PARTICLE_TAU) &&
            (context->mode.direction == PUMAS_MODE_FORWARD) &&
            (context->mode.decay == PUMAS_MODE_WEIGHT)) {
                return ERROR_MESSAGE(PUMAS_RETURN_DECAY_ERROR,
                    "`PUMAS_MODE_WEIGHT' mode is not valid for forward taus");
        }

        if ((context->accuracy <= 0) || (context->accuracy > 1)) {
                return ERROR_FORMAT(PUMAS_RETURN_ACCURACY_ERROR,
                    "bad accuracy value (expected a value in ]0,1], got %g)",
                    context->accuracy);
        }

        /* Get the start medium. */
        struct pumas_medium * medium;
        double step_max_medium;
        enum pumas_step step_max_type = context->medium(
            context, state, &medium, &step_max_medium);
        if (media != NULL) {
                media[0] = medium;
                media[1] = NULL;
        }
        if (medium == NULL) {
                if (event != NULL) *event = PUMAS_EVENT_MEDIUM;
                /* Register the start of the the track, if recording. */
                if (context->recorder != NULL)
                        record_state(context, medium, PUMAS_EVENT_MEDIUM |
                                PUMAS_EVENT_START | PUMAS_EVENT_STOP,
                            state);
                return PUMAS_RETURN_SUCCESS;
        } else if ((step_max_medium > 0.) &&
            (step_max_type == PUMAS_STEP_CHECK))
                step_max_medium += 0.5 * STEP_MIN;
        struct medium_locals locals = { { 0., { 0., 0., 0. }}, 0, physics };
        const double step_max_locals =
            transport_set_locals(context, medium, state, &locals);
        if ((step_max_locals > 0.) && (step_max_locals < step_max_medium))
                step_max_medium = step_max_locals;
        if (locals.api.density <= 0.) {
                ERROR_REGISTER_NEGATIVE_DENSITY(
                    physics->material_name[medium->material]);
                return ERROR_RAISE();
        }

        /* Randomise the lifetime, if required. */
        if (context->mode.decay == PUMAS_MODE_DECAY) {
                if (context->random == NULL) {
                        return ERROR_MESSAGE(PUMAS_RETURN_MISSING_RANDOM,
                            "no random engine specified");
                }
                const double u = context->random(context);
                context_->lifetime = state->time - physics->ctau * log(u);
        }

        /* Call the relevant transport engine. */
        int do_stepping = 1;
        enum pumas_event e = PUMAS_EVENT_NONE;
        if ((step_max_medium <= 0.) && (step_max_locals <= 0.) &&
            (context->mode.energy_loss <= PUMAS_MODE_CSDA)) {
                /* This is an infinite and uniform medium. */
                if ((context->mode.energy_loss == PUMAS_MODE_VIRTUAL) &&
                    ((context->event & PUMAS_EVENT_LIMIT) == 0)) {
                        return ERROR_MESSAGE(PUMAS_RETURN_MISSING_LIMIT,
                            "infinite medium without external limit(s)");
                } else if (
                    (context->mode.scattering == PUMAS_MODE_LONGITUDINAL) &&
                    (context->mode.energy_loss == PUMAS_MODE_CSDA)) {
                        do_stepping = 0;
                }
        }

        if (do_stepping) {
                /* Transport with a detailed stepping. */
                e = transport_with_stepping(physics, context, state, &medium,
                    &locals, step_max_medium, step_max_type, step_max_locals,
                    error_);
        } else {
                /* This is a purely deterministic case. */
                e = transport_with_csda(
                    physics, context, state, medium, &locals, error_);
        }

        if (event != NULL) *event = e;
        if (media != NULL) media[1] = medium;
        return ERROR_RAISE();
}

/* Public library function: transported particle info. */
enum pumas_return pumas_physics_particle(const struct pumas_physics * physics,
    enum pumas_particle * particle, double * lifetime, double * mass)
{
        ERROR_INITIALISE(pumas_physics_particle);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }
        if (particle != NULL) *particle = physics->particle;
        if (lifetime != NULL) *lifetime = physics->ctau;
        if (mass != NULL) *mass = physics->mass;

        return PUMAS_RETURN_SUCCESS;
}

/* Public library functions: elements handling. */
enum pumas_return pumas_physics_element_name(
    const struct pumas_physics * physics, int index, const char ** element)
{
        ERROR_INITIALISE(pumas_physics_element_name);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }
        if ((index < 0) || (index >= physics->n_elements)) {
                return ERROR_INVALID_ELEMENT(index);
        }
        *element = physics->element[index]->name;

        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_physics_element_index(
    const struct pumas_physics * physics, const char * element, int * index)
{
        ERROR_INITIALISE(pumas_physics_element_index);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }

        const int i = element_index(physics, element);
        if (i < 0) {
                return ERROR_FORMAT(PUMAS_RETURN_UNKNOWN_ELEMENT,
                    "unknown element `%s'", element);
        } else {
                *index = i;
        }

        return PUMAS_RETURN_SUCCESS;
}

int pumas_physics_element_length(const struct pumas_physics * physics)
{
        if (physics == NULL) return 0;
        return physics->n_elements;
}

enum pumas_return pumas_physics_element_properties(
    const struct pumas_physics * physics, int index, double * Z, double * A,
    double * I)
{
        ERROR_INITIALISE(pumas_physics_element_properties);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }
        if ((index < 0) || (index >= physics->n_elements)) {
                return ERROR_INVALID_ELEMENT(index);
        }

        if (Z != NULL) *Z = physics->element[index]->Z;
        if (A != NULL) *A = physics->element[index]->A;
        if (I != NULL) *I = physics->element[index]->I;

        return PUMAS_RETURN_SUCCESS;
}

/* Public library functions: materials handling. */
enum pumas_return pumas_physics_material_name(
    const struct pumas_physics * physics, int index, const char ** material)
{
        ERROR_INITIALISE(pumas_physics_material_name);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }
        if ((index < 0) || (index >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(index);
        }
        *material = physics->material_name[index];

        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_physics_material_index(
    const struct pumas_physics * physics, const char * material, int * index)
{
        ERROR_INITIALISE(pumas_physics_material_index);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }
        material_index(physics, material, index, error_);

        return ERROR_RAISE();
}

enum pumas_return pumas_physics_material_properties(
    const struct pumas_physics * physics, int material, int * length,
    double * density, double * I,
    struct pumas_physics_density_effect * density_effect, int * components,
    double * fractions)
{
        ERROR_INITIALISE(pumas_physics_material_properties);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }

        int i0 = physics->n_materials;
        if ((I != NULL) || (density_effect != NULL)) {
                i0 -= physics->n_composites;
        }

        if ((material < 0) || (material > i0 - 1)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        if (length != NULL) *length = physics->elements_in[material];
        if (density != NULL) *density = physics->material_density[material];
        if (I != NULL) *I = physics->material_I[material];
        if (density_effect != NULL) {
                struct pumas_physics_density_effect * const d =
                    physics->material_density_effect + material;
                memcpy(density_effect, d, sizeof(*density_effect));
        }
        int i;
        for (i = 0; i < physics->elements_in[material]; i++) {
                const struct material_component * component =
                    &physics->composition[material][i];
                if (components != NULL) components[i] = component->element;
                if (fractions != NULL) fractions[i] = component->fraction;
        }

        return PUMAS_RETURN_SUCCESS;
}

int pumas_physics_material_length(const struct pumas_physics * physics)
{
        if (physics == NULL) return 0;
        return physics->n_materials;
}

/* Public library functions: accessing and modifying composite materials. */
int pumas_physics_composite_length(const struct pumas_physics * physics)
{
        if (physics == NULL) return 0;
        return physics->n_composites;
}

enum pumas_return pumas_physics_composite_update(struct pumas_physics * physics,
    int material, const double * fractions)
{
        ERROR_INITIALISE(pumas_physics_composite_update);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }

        if (fractions == NULL) {
                return ERROR_MESSAGE(PUMAS_RETURN_VALUE_ERROR,
                    "NULL pointer for fractions");
        }

        const int i0 = physics->n_materials - physics->n_composites;
        if ((material < i0) || (material > physics->n_materials - 1)) {
                return ERROR_FORMAT(PUMAS_RETURN_INDEX_ERROR,
                    "invalid index for composite material [%d]", material);
        }

        const int icomp =
            material - physics->n_materials + physics->n_composites;
        int i;
        for (i = 0; i < physics->composite[icomp]->n_components; i++) {
                struct composite_component * component =
                    physics->composite[icomp]->component + i;
                const double d = fractions[i];
                component->fraction = (d > 0.) ? d : 0.;
        }

        const enum pumas_return rc =
            compute_composite_density(physics, material, error_);
        if ((rc != PUMAS_RETURN_SUCCESS) || (fractions == NULL))
                goto clean_and_exit;
        compute_composite(physics, material, error_);

clean_and_exit:
        /* Free temporary workspace and return. */
        compute_coulomb_parameters(physics, -1, -1, error_);
        return ERROR_RAISE();
}

enum pumas_return pumas_physics_composite_properties(
    const struct pumas_physics * physics, int material, int * length,
    int * components, double * fractions)
{
        ERROR_INITIALISE(pumas_physics_composite_properties);

        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }

        const int i0 = physics->n_materials - physics->n_composites;
        if ((material < i0) || (material > physics->n_materials - 1)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        const int icomp =
            material - physics->n_materials + physics->n_composites;
        if (length != NULL) *length = physics->composite[icomp]->n_components;
        int i;
        for (i = 0; i < physics->composite[icomp]->n_components; i++) {
                struct composite_component * component =
                    physics->composite[icomp]->component + i;
                if (components != NULL) components[i] = component->material;
                if (fractions != NULL) fractions[i] = component->fraction;
        }

        return PUMAS_RETURN_SUCCESS;
}

/* Public library functions: info on tabulations. */
enum pumas_return pumas_physics_table_value(
    const struct pumas_physics * physics, enum pumas_property property,
    enum pumas_mode scheme, int material, int row, double * value)
{
        ERROR_INITIALISE(pumas_physics_table_value);

        /* Check the input parameters. */
        *value = 0.;
        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }

        if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        } else if ((row < 0) || (row >= physics->n_energies)) {
                return ERROR_FORMAT(
                    PUMAS_RETURN_INDEX_ERROR, "invalid `row' index [%d]", row);
        }

        if (property == PUMAS_PROPERTY_KINETIC_ENERGY) {
                *value = *table_get_K(physics, row);
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_GRAMMAGE) {
                if ((scheme <= PUMAS_MODE_VIRTUAL) ||
                    (scheme >= PUMAS_MODE_DETAILED)) {
                        return ERROR_INVALID_SCHEME(scheme);
                }
                *value = *table_get_X(physics, scheme, material, row);
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_PROPER_TIME) {
                if ((scheme <= PUMAS_MODE_VIRTUAL) ||
                    (scheme >= PUMAS_MODE_DETAILED)) {
                        return ERROR_INVALID_SCHEME(scheme);
                }
                *value = *table_get_T(physics, scheme, material, row);
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_ENERGY_LOSS) {
                if ((scheme <= PUMAS_MODE_VIRTUAL) ||
                    (scheme >= PUMAS_MODE_DETAILED)) {
                        return ERROR_INVALID_SCHEME(scheme);
                }
                *value = *table_get_dE(physics, scheme, material, row);
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_MAGNETIC_ROTATION) {
                const double factor = LARMOR_FACTOR / physics->mass;
                *value =
                    *table_get_T(physics, PUMAS_MODE_CSDA, material, row) *
                    factor;
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_CROSS_SECTION) {
                *value = *table_get_CS(physics, material, row);
                return PUMAS_RETURN_SUCCESS;
        } else {
                return ERROR_FORMAT(PUMAS_RETURN_INDEX_ERROR,
                    "invalid `property' index [%d]", property);
        }
}

int pumas_physics_table_length(const struct pumas_physics * physics)
{
        if (physics == NULL) return 0;
        return physics->n_energies;
}

enum pumas_return pumas_physics_table_index(
    const struct pumas_physics * physics, enum pumas_property property,
    enum pumas_mode scheme, int material, double value, int * index)
{
        ERROR_INITIALISE(pumas_physics_table_index);

        /* Check some input parameters. */
        *index = 0;
        if (physics == NULL) {
                return ERROR_NOT_INITIALISED();
        }
        if ((material < 0) || (material >= physics->n_materials)) {
                return ERROR_INVALID_MATERIAL(material);
        }

        /* Get the tabulated value's index. */
        const double * table;
        if (property == PUMAS_PROPERTY_KINETIC_ENERGY)
                table = table_get_K(physics, 0);
        else if (property == PUMAS_PROPERTY_GRAMMAGE) {
                if ((scheme <= PUMAS_MODE_VIRTUAL) ||
                    (scheme >= PUMAS_MODE_DETAILED)) {
                        return ERROR_INVALID_SCHEME(scheme);
                }
                table = table_get_X(physics, scheme, material, 0);
        } else if (property == PUMAS_PROPERTY_PROPER_TIME) {
                if ((scheme <= PUMAS_MODE_VIRTUAL) ||
                    (scheme >= PUMAS_MODE_DETAILED)) {
                        return ERROR_INVALID_SCHEME(scheme);
                }
                table = table_get_T(physics, scheme, material, 0);
        } else if (property == PUMAS_PROPERTY_MAGNETIC_ROTATION) {
                value *= physics->mass / LARMOR_FACTOR;
                table = table_get_T(physics, PUMAS_MODE_CSDA, material, 0);
        } else {
                return ERROR_FORMAT(PUMAS_RETURN_INDEX_ERROR,
                    "invalid `property' index [%d]", property);
        }

        const int imax = physics->n_energies - 1;
        if (value < table[0]) {
                return ERROR_FORMAT(PUMAS_RETURN_VALUE_ERROR,
                    "out of range `value' [%.5lE < %.5lE]", value, table[0]);
        } else if (value == table[imax]) {
                *index = imax;
                return PUMAS_RETURN_SUCCESS;
        } else if (value > table[imax]) {
                *index = imax;
                return ERROR_FORMAT(PUMAS_RETURN_VALUE_ERROR,
                    "out of range `value' [%.5lE > %.5lE]", value, table[imax]);
        }

        int i1 = 0, i2 = imax;
        table_bracket(table, value, &i1, &i2);
        *index = i1;

        return PUMAS_RETURN_SUCCESS;
}

/* Low level routines: encapsulation of the tabulated CEL and DEL properties as
 * functions.
 *
 * Note that whenever a simulation context is required it is actualy optionnal
 * and can be `NULL`. The context is used for speeding up the interpolation
 * using a memory of previous indices.
 */
/**
 * Total grammage for a deterministic CEL as function of initial kinetic energy.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the total grammage in kg/m^2 otherwise `0`.
 */
double cel_grammage(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double kinetic)
{
        const int imax = physics->n_energies - 1;
        if (kinetic < *table_get_K(physics, 0)) return 0.;

        if (kinetic >= *table_get_K(physics, imax)) {
                /* Constant energy loss model. */
                const double a = *table_get_a_max(physics, material);
                const double b = *table_get_b_max(physics, scheme, material);
                const double K0 = *table_get_K(physics, imax);
                const double K1 = a / b + physics->mass;

                return *table_get_X(physics, scheme, material, imax) +
                    1. / b * log((kinetic + K1) / (K0 + K1));
        }

        /* Interpolation. */
        return table_interpolate(physics, context, table_get_K(physics, 0),
            table_get_X(physics, scheme, material, 0), kinetic);
}

/**
 * Total grammage for a deterministic CEL as function of total proper time.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param time     The total proper time variation.
 * @return On success, the total grammage in kg/m^2 otherwise `0`.
 */
double cel_grammage_as_time(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double time)
{
        const int imax = physics->n_energies - 1;
        if (time < *table_get_T(physics, scheme, material, 0)) return 0.;

        if (time >= *table_get_T(physics, scheme, material, imax)) {
                /* Constant energy loss model. */
                const double a = *table_get_a_max(physics, material);
                const double b = *table_get_b_max(physics, scheme, material);
                const double E0 = *table_get_K(physics, imax) + physics->mass;
                const double tmax =
                    *table_get_T(physics, scheme, material, imax);
                const double r =
                    E0 / (a + b * E0) * exp(a * (time - tmax) / physics->mass);
                const double kinetic = a * r / (1. - b * r) - physics->mass;
                const double K1 = a / b + physics->mass;
                return *table_get_X(physics, scheme, material, imax) +
                    1. / b * log((kinetic + K1) / (E0 - physics->mass + K1));
        }

        /* Interpolation. */
        return table_interpolate(physics, context,
            table_get_T(physics, scheme, material, 0),
            table_get_X(physics, scheme, material, 0), time);
}

/**
 * Total proper time for a deterministic CEL in a homogeneous medium.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the normalised proper time in kg/m^2 otherwise `0`.
 */
double cel_proper_time(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double kinetic)
{
        const int imax = physics->n_energies - 1;
        if (kinetic < *table_get_K(physics, 0)) return 0.;

        if (kinetic >= *table_get_K(physics, imax)) {
                /* Constant energy loss model. */
                const double a = *table_get_a_max(physics, material);
                const double b = *table_get_b_max(physics, scheme, material);
                const double E0 = *table_get_K(physics, imax) + physics->mass;
                const double E1 = kinetic + physics->mass;

                return *table_get_T(physics, scheme, material, imax) +
                    physics->mass / a *
                    log((E1 / E0) * (a + b * E0) / (a + b * E1));
        }

        /* Interpolation. */
        return table_interpolate(physics, context, table_get_K(physics, 0),
            table_get_T(physics, scheme, material, 0), kinetic);
}

/**
 * The initial kinetic energy for a given total grammage assuming a determistic
 * CEL.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param grammage The total grammage depth.
 * @return On success, the initial kinetic energy in GeV otherwise `0`.
 */
double cel_kinetic_energy(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double grammage)
{
        const int imax = physics->n_energies - 1;
        if (grammage < *table_get_X(physics, scheme, material, 0)) return 0.;

        if (grammage >= *table_get_X(physics, scheme, material, imax)) {
                /* Constant energy loss model. */
                const double a = *table_get_a_max(physics, material);
                const double b = *table_get_b_max(physics, scheme, material);
                const double K0 = *table_get_K(physics, imax);
                const double K1 = a / b + physics->mass;

                return (K0 + K1) * exp(b * (grammage -
                                               *table_get_X(physics, scheme,
                                                   material, imax))) -
                    K1;
        }

        /* Interpolation. */
        return table_interpolate(physics, context,
            table_get_X(physics, scheme, material, 0), table_get_K(physics, 0),
            grammage);
}

/**
 * The average CEL.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the CEL in GeV/(kg/m^2) otherwise `0`.
 */
double cel_energy_loss(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double kinetic)
{
        const int imax = physics->n_energies - 1;
        if (kinetic < *table_get_K(physics, 0)) return 0.;

        if (kinetic >= *table_get_K(physics, imax)) {
                /* Constants energy loss model. */
                return *table_get_a_max(physics, material) +
                    *table_get_b_max(physics, scheme, material) *
                    (kinetic + physics->mass);
        }

        /* Interpolation. */
        return table_interpolate(physics, context, table_get_K(physics, 0),
            table_get_dE(physics, scheme, material, 0), kinetic);
}

/**
 * The normalised magnetic rotation angle for a deterministic CEL.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The index of the material in which the
 *                 particle travels.
 * @param kinetic  The initial kinetic energy, in GeV.
 * @return The normalised rotation angle in kg/m^2/T.
 *
 * The magnetic rotation angle is proportional to the proper time integral.
 * Therefore it is computed from the proper time table.
 */
double cel_magnetic_rotation(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic)
{
        const int imax = physics->n_energies - 1;
        const double factor = LARMOR_FACTOR / physics->mass;
        double * const T = table_get_T(physics, PUMAS_MODE_CSDA, material, 0);
        if (kinetic <= *table_get_K(physics, 0)) return T[imax] * factor;

        if (kinetic >= *table_get_K(physics, imax)) {
                /*
                 * Neglect any magnetic rotation above the max tabulated value
                 * of the kinetic energy.
                 */
                return 0.;
        }

        /* Interpolation. */
        return (T[imax] - table_interpolate(physics, context,
                              table_get_K(physics, 0), T, kinetic)) *
            factor;
}

/**
 * The total cross section for inelastic DELs.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the DEL total cross section, otherwise `0`.
 */
double del_cross_section(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic)
{
        const int imax = physics->n_energies - 1;
        if (kinetic < *table_get_K(physics, 0)) return 0.;

        if (kinetic >= *table_get_K(physics, imax)) {
                /* Constant cross section model. */
                return *table_get_CS(physics, material, imax);
        }

        /* Interpolation. */
        return table_interpolate(physics, context, table_get_K(physics, 0),
            table_get_CS(physics, material, 0), kinetic);
}

/**
 * Total number of interaction lengths for a deterministic CEL.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the number of interaction lenths, otherwise `0`.
 */
double del_interaction_length(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic)
{
        const int imax = physics->n_energies - 1;
        if (kinetic < *table_get_K(physics, 0)) return 0.;

        if (kinetic >= *table_get_K(physics, imax)) {
                /* constant loss model. */
                const double k0 = *table_get_K(physics, imax);
                const double a0 = *table_get_a_max(physics, material);
                const double b0 =
                    *table_get_b_max(physics, PUMAS_MODE_HYBRID, material);
                const double cs = *table_get_CS(physics, material, imax);
                const double dZ =
                    cs / b0 * log((a0 + b0 * (kinetic + physics->mass)) /
                                  (a0 + b0 * (k0 + physics->mass)));
                return *table_get_NI_in(physics, material, imax) + dZ;
        }

        /* Interpolation. */
        return table_interpolate(physics, context, table_get_K(physics, 0),
            table_get_NI_in(physics, material, 0), kinetic);
}

/**
 * Initial kinetic energy for a given number of interaction lengths.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param nI       The number of interaction lengths.
 * @return On success, the initial kinetic energy, otherwise `0`.
 */
double del_kinetic_from_interaction_length(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double nI)
{
        const int imax = physics->n_energies - 1;
        if (nI < *table_get_NI_in(physics, material, 0)) return 0.;

        if (nI >= *table_get_NI_in(physics, material, imax)) {
                /* constant loss model. */
                const double k0 = *table_get_K(physics, imax);
                const double a0 = *table_get_a_max(physics, material);
                const double b0 =
                    *table_get_b_max(physics, PUMAS_MODE_HYBRID, material);
                const double cs = *table_get_CS(physics, material, imax);
                const double nI0 = *table_get_NI_in(physics, material, imax);
                return ((a0 + b0 * (k0 + physics->mass)) *
                               exp(b0 * (nI - nI0) / cs) -
                           a0) /
                    b0 -
                    physics->mass;
        }

        /* Interpolation. */
        return table_interpolate(physics, context,
            table_get_NI_in(physics, material, 0), table_get_K(physics, 0), nI);
}

/**
 * Total number of EHS interaction lengths for a deterministic CEL.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the number of interaction lengths, otherwise `0`.
 */
double ehs_interaction_length(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double kinetic)
{
        const int imax = physics->n_energies - 1;
        if (kinetic < *table_get_K(physics, 0)) return 0.;

        if (kinetic >= *table_get_K(physics, imax)) {
                /* linear extrapolation. */
                const double k0 = *table_get_K(physics, imax - 1);
                const double k1 = *table_get_K(physics, imax);
                const double nI0 =
                    *table_get_NI_el(physics, scheme, material, imax - 1);
                const double nI1 =
                    *table_get_NI_el(physics, scheme, material, imax);
                return nI1 + (kinetic - k1) * (nI1 - nI0) / (k1 - k0);
        }

        /* Interpolation. */
        return table_interpolate(physics, context, table_get_K(physics, 0),
            table_get_NI_el(physics, scheme, material, 0), kinetic);
}

/**
 * Initial kinetic energy for a total number of EHS interaction lengths.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param nI       The number of interaction lengths.
 * @return On success, the initial kinetic energy, otherwise `0`.
 */
double ehs_kinetic_from_interaction_length(const struct pumas_physics * physics,
    struct pumas_context * context, enum pumas_mode scheme, int material,
    double nI)
{
        const int imax = physics->n_energies - 1;
        if (nI < *table_get_NI_el(physics, scheme, material, 0)) return 0.;

        if (nI >= *table_get_NI_el(physics, scheme, material, imax)) {
                /* linear extrapolation. */
                const double k0 = *table_get_K(physics, imax - 1);
                const double k1 = *table_get_K(physics, imax);
                const double nI0 =
                    *table_get_NI_el(physics, scheme, material, imax - 1);
                const double nI1 =
                    *table_get_NI_el(physics, scheme, material, imax);
                return k1 + (nI - nI1) * (k1 - k0) / (nI1 - nI0);
        }

        /* Interpolation. */
        return table_interpolate(physics, context,
            table_get_NI_el(physics, scheme, material, 0),
            table_get_K(physics, 0), nI);
}

/*
 * Low level routines: generic and specific interpolations of various
 * tabulated data.
 */
/**
 * Interpolation of the Multiple SCattering (MSC) parameters.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The propagation material.
 * @param kinetic  The kinetic energy.
 * @param mu0      The interpolated angular cuttof for Coulomb scattering.
 * @param invlb1   The interpolated 1st transport inverse grammage.
 *
 * Interpolate the cutt-off angle for Coulomb scattering and the total 1st
 * transport path length for MSC.
 */
void table_get_msc(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic, double * mu0,
    double * invlb1)
{
        const int imax = physics->n_energies - 1;
        if (kinetic < *table_get_K(physics, 1)) {
                *mu0 = *table_get_Mu0(physics, material, 1);
                /* Use asymptotic limit as lb1 ~ sqrt(kinetic). */
                *invlb1 = *table_get_Ms1(physics, material, 1) *
                    sqrt((*table_get_K(physics, 1)) / kinetic);
        } else if (kinetic >= *table_get_K(physics, imax)) {
                *mu0 = *table_get_Mu0(physics, material, imax);
                /* Use asymptotic limit as lb1 ~ kinetic. */
                *invlb1 = *table_get_Ms1(physics, material, imax) *
                    (*table_get_K(physics, imax)) / kinetic;
        } else {
                const int i1 = table_index(
                    physics, context, table_get_K(physics, 0), kinetic);
                const int i2 = i1 + 1;
                double h = (kinetic - *table_get_K(physics, i1)) /
                    (*table_get_K(physics, i2) - *table_get_K(physics, i1));
                *mu0 = *table_get_Mu0(physics, material, i1) +
                    h * (*table_get_Mu0(physics, material, i2) -
                            *table_get_Mu0(physics, material, i1));
                *invlb1 = *table_get_Ms1(physics, material, i1) +
                    h * (*table_get_Ms1(physics, material, i2) -
                            *table_get_Ms1(physics, material, i1));
        }
}

/**
 * Interpolate an arbitrary tabulated property.
 *
 * @param Physics Handle for physics tables.
 * @param context The simulation context.
 * @param table_X Table of x values.
 * @param table_Y Table of y values.
 * @param x       Point at which the interpolant is evaluated.
 * @return The interpolated value.
 *
 * Compute a linear interpolant of property Y at x given table_Y values as
 * table_X values. Note that table_X must be strictly monotone but not
 * necessarily with a constant stepping. The `context` parameter is used for
 * mapping the index of x in table_X from a memory. It can be `NULL`, though
 * providing a `context` leads to a speed up on average for Monte-Carlo
 * stepping.
 * **Warning** : there is no bound check. x must be in the range of table_X.
 */
double table_interpolate(const struct pumas_physics * physics,
    struct pumas_context * context, const double * table_X,
    const double * table_Y, double x)
{
        const int i1 = table_index(physics, context, table_X, x);
        const int i2 = i1 + 1;
        double h = (x - table_X[i1]) / (table_X[i2] - table_X[i1]);
        return table_Y[i1] + h * (table_Y[i2] - table_Y[i1]);
}

/**
 * Find the index closest to `value`, from below.
 *
 * @param Physics Handle for physics tables.
 * @param context The simulation context.
 * @param table   The tabulated values.
 * @param value   The value to bracket.
 * @return In case of success the closest index from below is returned,
 * otherwise -1 in case of underflow or `pumas_physics::n_energies-1` in
 * case of overflow.
 *
 * Compute the table index for the given entry `value` using a dichotomy
 * algorithm. If a `context` is not `NULL`, `value` is checked against the
 * last used indices in the table, before doing the dichotomy search.
 */
int table_index(const struct pumas_physics * physics,
    struct pumas_context * context, const double * table, double value)
{
        int * last = NULL;
        if (context != NULL) {
                /* Check if the last used indices are still relevant. */
                struct simulation_context * const ctx =
                    (struct simulation_context * const)context;
                if (table == table_get_K(physics, 0))
                        last = ctx->index_K_last;
                else
                        last = ctx->index_X_last;

                if ((value >= table[last[0]]) && (value < table[last[0] + 1]))
                        return last[0];
                if ((value >= table[last[1]]) && (value < table[last[1] + 1]))
                        return last[1];
        }

        /* Check the boundaries. */
        const int imax = physics->n_energies - 1;
        if (value < table[0]) return -1;
        if (value >= table[imax]) return imax;

        /* Bracket the value. */
        int i1 = 0, i2 = imax;
        table_bracket(table, value, &i1, &i2);

        if (context != NULL) {
                /* Update the last used indices. */
                if (i1 != last[0])
                        last[0] = i1;
                else if (i1 != last[1])
                        last[1] = i1;
        }
        return i1;
}

/**
 * Recursive bracketing of a table value.
 *
 * @param table The tabulated values.
 * @param value The value to bracket.
 * @param p1 The lower bracketing index.
 * @param p2 The upper bracketing index.
 * @return At final return p1 points to the closest from below index
 * to `value` and p2 to the closest from above.
 *
 * Refine the bracketing indices [p1, p2] of `value` in table using a
 * recursive procedure.
 */
void table_bracket(const double * table, double value, int * p1, int * p2)
{
        int i3 = (*p1 + *p2) / 2;
        if (value >= table[i3])
                *p1 = i3;
        else
                *p2 = i3;
        if (*p2 - *p1 >= 2) table_bracket(table, value, p1, p2);
}

/*
 * Low level routines: inlined encapsulations of table accesses.
 *
 * The tables are stored linearly in the shared data segment, allocated on the
 * heap. The following routines provide pointers on the corresponding table
 * elements.
 */
/**
 * Encapsulation of the kinetic energy table.
 *
 * @param Physics  Handle for physics tables.
 * @param row The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_K(const struct pumas_physics * physics, int row)
{
        return physics->table_K + row;
}

/**
 * Encapsulation of the total grammage table.
 *
 * @param Physics  Handle for physics tables.
 * @param scheme   The energy loss scheme.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_X(
    const struct pumas_physics * physics, int scheme, int material, int row)
{
        scheme = (scheme > PUMAS_MODE_HYBRID) ? PUMAS_MODE_HYBRID : scheme;
        return physics->table_X +
            (scheme * physics->n_materials + material) * physics->n_energies +
            row;
}

/**
 * Encapsulation of the proper time table.
 *
 * @param Physics  Handle for physics tables.
 * @param scheme   The energy loss scheme.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_T(
    const struct pumas_physics * physics, int scheme, int material, int row)
{
        scheme = (scheme > PUMAS_MODE_HYBRID) ? PUMAS_MODE_HYBRID : scheme;
        return physics->table_T +
            (scheme * physics->n_materials + material) * physics->n_energies +
            row;
}

/**
 * Encapsulation of the dE/dX table.
 *
 * @param Physics  Handle for physics tables.
 * @param scheme   The energy loss scheme.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_dE(
    const struct pumas_physics * physics, int scheme, int material, int row)
{
        scheme = (scheme > PUMAS_MODE_HYBRID) ? PUMAS_MODE_HYBRID : scheme;
        return physics->table_dE +
            (scheme * physics->n_materials + material) * physics->n_energies +
            row;
}

/**
 * Encapsulation of the number of EHS interaction lengths.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_NI_el(
    const struct pumas_physics * physics, int scheme, int material, int row)
{
        scheme = (scheme >= PUMAS_MODE_HYBRID) ? PUMAS_MODE_HYBRID :
                                                   PUMAS_MODE_CSDA;
        return physics->table_NI_el +
            (scheme * physics->n_materials + material) * physics->n_energies +
            row;
}

/**
 * Encapsulation of the number of interaction lengths for inelastic DELs.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_NI_in(
    const struct pumas_physics * physics, int material, int row)
{
        return physics->table_NI_in + material * physics->n_energies + row;
}

/**
 * Encapsulation of the cross section table, for inelastic DELs.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_CS(
    const struct pumas_physics * physics, int material, int row)
{
        return physics->table_CS + material * physics->n_energies + row;
}

/**
 * Encapsulation of the fractional cross-sections table, for inelastic DELs.
 *
 * @param Physics   Handle for physics tables.
 * @param process   The process index.
 * @param component The component index.
 * @param row       The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_CSf(
    const struct pumas_physics * physics, int process, int component, int row)
{
        return physics->table_CSf +
            (process * physics->n_components + component) *
            physics->n_energies +
            row;
}

/**
 * Encapsulation of the cross-sections normalisation table.
 *
 * @param Physics Handle for physics tables.
 * @param process The process index.
 * @param element The element index.
 * @param row     The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_CSn(
    const struct pumas_physics * physics, int process, int element, int row)
{
        return physics->table_CSn +
            (process * physics->n_elements + element) * physics->n_energies +
            row;
}

/**
 * Encapsulation of the fractional lower threshold for inelastic DELs.
 *
 * @param Physics Handle for physics tables.
 * @param process The process index.
 * @param element The element index.
 * @param row     The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_Xt(
    const struct pumas_physics * physics, int process, int element, int row)
{
        return physics->table_Xt +
            (process * physics->n_elements + element) * physics->n_energies +
            row;
}

/**
 * Encapsulation of the lower kinetic threshold for inelatic DELs.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @return A pointer to the table element.
 */
double * table_get_Kt(const struct pumas_physics * physics, int material)
{
        return physics->table_Kt + material;
}

/**
 * Encapsulation of the temporary average CEL table, element wise.
 *
 * @param Physics  Handle for physics tables.
 * @param property The index of the tabulated property.
 * @param process  The process index.
 * @param element  The element index.
 * @param row      The kinetic energy row index.
 * @param table    The temporary table.
 * @return A pointer to the table element.
 *
 * The `property` index controls the tabulated property as following.
 * `property = 0` is the total CEL and `property = 1` is the restricted
 * CEL, splitted with DELs.
 */
double * table_get_cel(const struct pumas_physics * physics, int process,
    int element, int row, double * table)
{
        return table +
            (process * physics->n_elements + element) * physics->n_energies +
            row;
}

/**
 * Encapsulation of the magnetic deflection table, when using CSDA in a
 * homogeneous medium of infinite extension.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @param order    The order of the Taylor expansion.
 * @param row      The kinetic row index.
 * @return A pointer to the table element.
 */
double * table_get_Li(
    const struct pumas_physics * physics, int material, int order, int row)
{
        return physics->table_Li +
            (material * N_LARMOR_ORDERS + order) * physics->n_energies + row;
}

/**
 * Encapsulation of the last tabulated ionisation dE/dX.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @return A pointer to the table element.
 */
double * table_get_a_max(const struct pumas_physics * physics, int material)
{
        return physics->table_a_max + material;
}

/**
 * Encapsulation of the last tabulated radiative dE/dX parameter.
 *
 * @param Physics  Handle for physics tables.
 * @param scheme   The energy loss scheme.
 * @param material The material index.
 * @return A pointer to the table element.
 */
double * table_get_b_max(
    const struct pumas_physics * physics, int scheme, int material)
{
        scheme = (scheme > PUMAS_MODE_HYBRID) ? PUMAS_MODE_HYBRID : scheme;
        return physics->table_b_max + scheme * physics->n_materials + material;
}

/**
 * Encapsulation of the angular cutoff for Coulomb scattering.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_Mu0(
    const struct pumas_physics * physics, int material, int row)
{
        return physics->table_Mu0 + material * physics->n_energies + row;
}

/**
 * Encapsulation of the interaction length for EHS events.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_Lb(
    const struct pumas_physics * physics, int material, int row)
{
        return physics->table_Lb + material * physics->n_energies + row;
}

/*!
 * Encapsulation of the Multiple SCattering (MSC) 1st transport path length.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_Ms1(
    const struct pumas_physics * physics, int material, int row)
{
        return physics->table_Ms1 + material * physics->n_energies + row;
}

/*!
 * Encapsulation of the temporary MSC 1st transport path length, element wise.
 *
 * @param Physics  Handle for physics tables.
 * @param element  The material index.
 * @param row      The kinetic energy row index.
 * @param table    The temporary table.
 * @return A pointer to the table element.
 */
double * table_get_ms1(
    const struct pumas_physics * physics, int element, int row, double * table)
{
        return table + element * physics->n_energies + row;
}

/**
 * Encapsulation of the polynomial coefficients of the DCS model.
 *
 * @param Physics Handle for physics tables.
 * @param element The element index.
 * @param process The process index.
 * @param row     The kinetic energy row index.
 * @return A pointer to the table element.
 */
float * table_get_dcs_coeff(const struct pumas_physics * physics,
    const struct atomic_element * element, int process, int kinetic)
{
        const int n =
            DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1 + DCS_SAMPLING_N;
        return element->dcs_data +
            (process * (physics->n_energies - physics->dcs_model_offset) +
                kinetic) *
            n;
}

/**
 * Encapsulation of the tabulated DCS values.
 *
 * @param Physics Handle for physics tables.
 * @param element The element index.
 * @param process The process index.
 * @param row     The kinetic energy row index.
 * @return A pointer to the table element.
 */
float * table_get_dcs_value(const struct pumas_physics * physics,
    const struct atomic_element * element, int process, int kinetic)
{
        const int n = DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1;
        return element->dcs_data +
            (process * (physics->n_energies - physics->dcs_model_offset) +
                kinetic) *
            (n + DCS_SAMPLING_N) +
            n;
}

/*
 * Low level routines: propagation.
 */
/**
 * CSDA propagation routine for a uniform and infinite medium.
 *
 * @param Physics      Handle for physics tables.
 * @param context      The simulation context.
 * @param state        The initial/final state.
 * @param medium_index The index of the propagation medium.
 * @param locals       Handle for the local properties of the uniform medium.
 * @param error_       The error data.
 * @return The end condition event.
 *
 * Transport the particle through a uniform medium using the CSDA. At output
 * the final kinetic energy, the effective distance traveled and the
 * proper time spent are updated.
 *
 * Backward propagation allows to compute the required minimum kinetic
 * energy for crossing the medium. The corresponding distance and
 * proper time are also given.
 * **Warning** : the initial state must have been initialized before the call.
 */
enum pumas_event transport_with_csda(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium * medium, struct medium_locals * locals,
    struct error_context * error_)
{
        /* Unpack and backup the initial state. */
        const int material = medium->material;
        const double density = locals->api.density;
        const double ki = state->energy;
        const double di = state->distance;
        const double ti = state->time;
        const double xi =
            cel_grammage(physics, context, PUMAS_MODE_CSDA, material, ki);
        const double Ti =
            cel_proper_time(physics, context, PUMAS_MODE_CSDA, material, ki);

        /* Register the start of the the track, if recording. */
        enum pumas_event event = PUMAS_EVENT_NONE;
        struct pumas_recorder * recorder = context->recorder;
        const int record = (recorder != NULL);
        if (record)
                record_state(context, medium, event | PUMAS_EVENT_START, state);

        /* Get the end point with CSDA. */
        double xB;
        if (context->event & PUMAS_EVENT_LIMIT_GRAMMAGE) {
                xB = context->limit.grammage - state->grammage;
                event = PUMAS_EVENT_LIMIT_GRAMMAGE;
        } else
                xB = DBL_MAX;
        if (context->event & PUMAS_EVENT_LIMIT_DISTANCE) {
                const double xD = density * (context->limit.distance - di);
                if (xD < xB) {
                        xB = xD;
                        event = PUMAS_EVENT_LIMIT_DISTANCE;
                }
        }
        int decayed = 0;
        double time_max = (context->event & PUMAS_EVENT_LIMIT_TIME) ?
            context->limit.time : 0.;
        if (context->mode.decay == PUMAS_MODE_DECAY) {
                struct simulation_context * c =
                    (struct simulation_context *)context;
                if ((time_max <= 0.) || (c->lifetime < time_max)) {
                        time_max = c->lifetime;
                        decayed = 1;
                }
        }
        if (time_max > 0.) {
                const double dt = time_max - ti;
                if (dt <= 0.) {
                        xB = 0.;
                        event = PUMAS_EVENT_LIMIT_TIME;
                } else {
                        const double sgn =
                            (context->mode.direction == PUMAS_MODE_FORWARD) ?
                            1. : -1.;
                        const double Tf = Ti - sgn * dt * density;
                        if (Tf > 0.) {
                                const double xT = fabs(
                                    xi - cel_grammage_as_time(physics, context,
                                             PUMAS_MODE_CSDA, material, Tf));
                                if (xT < xB) {
                                        xB = xT;
                                        event = PUMAS_EVENT_LIMIT_TIME;
                                }
                        }
                }
        }
        if (xB <= 0.) return event;

        double xf, kf;
        if (context->mode.direction == PUMAS_MODE_FORWARD) {
                if (xi > xB) {
                        xf = xi - xB;
                        kf = cel_kinetic_energy(
                            physics, context, PUMAS_MODE_CSDA, material, xf);
                } else {
                        xf = kf = 0.;
                        event = PUMAS_EVENT_LIMIT_ENERGY;
                }

                if (context->event & PUMAS_EVENT_LIMIT_ENERGY) {
                        if (ki <= context->limit.energy)
                                return PUMAS_EVENT_LIMIT_ENERGY;
                        if (kf < context->limit.energy) {
                                kf = context->limit.energy;
                                xf = cel_grammage(physics, context,
                                    PUMAS_MODE_CSDA, material, kf);
                                event = PUMAS_EVENT_LIMIT_ENERGY;
                        }
                }
        } else {
                xf = xB + xi;
                kf = cel_kinetic_energy(
                    physics, context, PUMAS_MODE_CSDA, material, xf);
                if (context->event & PUMAS_EVENT_LIMIT_ENERGY) {
                        if (ki >= context->limit.energy)
                                return PUMAS_EVENT_LIMIT_ENERGY;
                        if (kf > context->limit.energy) {
                                kf = context->limit.energy;
                                xf = cel_grammage(physics, context,
                                    PUMAS_MODE_CSDA, material, kf);
                                event = PUMAS_EVENT_LIMIT_ENERGY;
                        }
                }
        }

        /* Update the end point statistics. */
        const double distance = fabs(xf - xi) / density;
        state->energy = kf;
        if (event == PUMAS_EVENT_LIMIT_DISTANCE)
                state->distance = context->limit.distance;
        else
                state->distance += distance;
        if (event == PUMAS_EVENT_LIMIT_GRAMMAGE)
                state->grammage = context->limit.grammage;
        else
                state->grammage += distance * density;
        if (event == PUMAS_EVENT_LIMIT_TIME) {
                state->time = time_max;
                state->decayed = decayed;
                if (decayed) event = PUMAS_EVENT_VERTEX_DECAY;
        } else
                state->time += fabs(Ti - cel_proper_time(physics, context,
                                             PUMAS_MODE_CSDA, material, kf)) /
                    density;
        if (context->mode.direction == PUMAS_MODE_BACKWARD)
                state->weight *= cel_energy_loss(physics, context,
                                     PUMAS_MODE_CSDA, material, kf) /
                    cel_energy_loss(physics, context, PUMAS_MODE_CSDA,
                                     material, ki);
        if (context->mode.decay == PUMAS_MODE_WEIGHT)
                state->weight *= exp(-fabs(ti - state->time) / physics->ctau);

        /* Update the position and direction. */
        if ((locals->magnetized != 0)) {
                if (transport_csda_deflect(physics, context, state, medium,
                        locals, ki, distance, error_) != PUMAS_RETURN_SUCCESS)
                        return event;
        } else {
                double path;
                if (context->mode.direction == PUMAS_MODE_FORWARD)
                        path = distance;
                else
                        path = -distance;
                state->position[0] += path * state->direction[0];
                state->position[1] += path * state->direction[1];
                state->position[2] += path * state->direction[2];
        }

        /* Register the end of the track, if recording. */
        if (record)
                record_state(context, medium, event | PUMAS_EVENT_STOP, state);

        return event;
}

/**
 * Apply the magnetic deflection in CSDA scheme and uniform medium.
 *
 * @param Physics      Handle for physics tables.
 * @param context      The simulation context.
 * @param state        The final state.
 * @param medium       The propagation medium.
 * @param locals       Handle for the local properties of the uniform medium.
 * @param ki           The initial kinetic energy.
 * @param distance     The travelled distance.
 * @param error_       The error data.
 * @return #PUMAS_RETURN_SUCCESS on success, #PUMAS_ERROR otherwise.
 *
 * Compute the magnetic deflection between two kinetic energies using the CSDA.
 * The final kinetic energy is read from the state. At return the final state
 * position and direction are updated.
 */
enum pumas_return transport_csda_deflect(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium * medium, struct medium_locals * locals, double ki,
    double distance, struct error_context * error_)
{
        /* Unpack arguments */
        const double charge = state->charge;
        const double kf = state->energy;
        double * const position = state->position;
        double * const direction = state->direction;
        const int material = medium->material;
        const double density = locals->api.density;
        const double * const magnet = locals->api.magnet;

        /* Compute the local basis and the transverse magnetic field
         * amplitude. Note that ex, ey and ez are not normed but already
         * multiplied by the cos(theta) and sin(theta) factors.
         */
        double b0, ex[3], ey[3], ez[3];
        {
                b0 = sqrt(magnet[0] * magnet[0] + magnet[1] * magnet[1] +
                    magnet[2] * magnet[2]);
                double b0_i = 1.0 / b0;

                /* Negative charge convention for particles. */
                ey[0] = (direction[2] * magnet[1] - direction[1] * magnet[2]) *
                    b0_i;
                ey[1] = (direction[0] * magnet[2] - direction[2] * magnet[0]) *
                    b0_i;
                ey[2] = (direction[1] * magnet[0] - direction[0] * magnet[1]) *
                    b0_i;

                {
                        double d = b0_i * b0_i * (magnet[0] * direction[0] +
                                                     magnet[1] * direction[1] +
                                                     magnet[2] * direction[2]);
                        ez[0] = magnet[0] * d;
                        ez[1] = magnet[1] * d;
                        ez[2] = magnet[2] * d;
                }

                ex[0] = direction[0] - ez[0];
                ex[1] = direction[1] - ez[1];
                ex[2] = direction[2] - ez[2];
        }

        /* Update the position and direction. */
        double dx, dy, dz, dp;
        {
                /* Negative charge convention. */
                double a = -b0 * charge / density;
                double pi =
                    a * cel_magnetic_rotation(physics, context, material, ki);
                double pf =
                    a * cel_magnetic_rotation(physics, context, material, kf);
                dp = pf - pi;

                const double ps = fabs(pi + pf);
                if (ps == 0.) return PUMAS_RETURN_SUCCESS;
                const double dp0 = fabs(dp);
                if ((dp0 < 1E-02) && (dp0 / ps < 1E-02)) {
                        /* Use the leading order expressions for a small
                         * rotation varying linearly with the distance.
                         */
                        if (ki >= kf)
                                dz = distance;
                        else
                                dz = -distance;

                        const double sixth = 1 / 6.0;
                        dx = dz * (1.0 - sixth * dp * dp);
                        dy = 0.5 * dz * dp;
                } else {
                        double xi, yi, zi, xf, yf, zf;
                        if (csda_magnetic_transport(physics, context, material,
                                density, b0, charge, ki, pi, &xi, &yi, &zi,
                                error_) != PUMAS_RETURN_SUCCESS)
                                return error_->code;
                        if (csda_magnetic_transport(physics, context, material,
                                density, b0, charge, kf, pf, &xf, &yf, &zf,
                                error_) != PUMAS_RETURN_SUCCESS)
                                return error_->code;

                        /* Rotate back to the initial frame. */
                        {
                                const double si = sin(pi);
                                const double ci = cos(pi);
                                xf -= xi;
                                yf -= yi;
                                dx = ci * xf + si * yf;
                                dy = -si * xf + ci * yf;
                        }
                        dz = zf - zi;
                }
        }
        position[0] += dx * ex[0] + dy * ey[0] + dz * ez[0];
        position[1] += dx * ex[1] + dy * ey[1] + dz * ez[1];
        position[2] += dx * ex[2] + dy * ey[2] + dz * ez[2];

        const double s = sin(dp);
        const double c = cos(dp);
        direction[0] = c * ex[0] + s * ey[0] + ez[0];
        direction[1] = c * ex[1] + s * ey[1] + ez[1];
        direction[2] = c * ex[2] + s * ey[2] + ez[2];

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Compute the total magnetic transport within CSDA, for a uniform medium.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param density  The material density.
 * @param magnet   The transverse magnetic field amplitude.
 * @param charge   The particle charge.
 * @param kinetic  The kinetic energy.
 * @param phase    The magnetic angular rotation.
 * @param x        The total magnetic transport for x coordinate.
 * @param y        The total magnetic transport for y coordinate.
 * @param z        The total magnetic transport for z coordinate.
 * @param error_   The error data.
 * @return `PUMAS_RETURN_SUCCESS` on success or `PUMAS_ERROR` otherwise, i.e. if
 * an
 * out of bound error ooccured.
 */
enum pumas_return csda_magnetic_transport(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double density, double magnet,
    double charge, double kinetic, double phase, double * x, double * y,
    double * z, struct error_context * error_)
{
        const int imax = physics->n_energies - 1;
        double k0 = kinetic;
        int i1, i2;
        if (kinetic <= *table_get_K(physics, 0)) {
                k0 = *table_get_K(physics, 0);
                i1 = 0;
                i2 = 1;
        } else if (kinetic >= *table_get_K(physics, imax)) {
                /* Neglect deflection at very high energy. */
                *y = 0.0;
                *x = *z = -(cel_grammage(physics, context, PUMAS_MODE_CSDA,
                                material, kinetic) -
                    *table_get_X(physics, 0, material, imax));
                return PUMAS_RETURN_SUCCESS;
        } else {
                /* Interpolation of the step starting energy. */
                i1 = table_index(physics, context, table_get_K(physics, 0), k0);
                i2 = i1 + 1;
        }

        /* Check that the phase does not exceed the sine & cosine
         * interpolation range.
         */
        const double max_phi = 2.0 * M_PI;
        const double poly_x[N_LARMOR_ORDERS] = { 1.000000000e+000,
                0.000000000e+000, -5.000000000e-001, -1.048969541e-002,
                5.597396438e-002, -6.401353612e-003, -4.495849930e-004,
                6.583532135e-005 };
        const double poly_y[N_LARMOR_ORDERS] = { 0.000000000e+000,
                1.000000000e+000, -1.400755143e-002, -1.350173383e-001,
                -2.838778336e-002, 2.123237056e-002, -3.094290091e-003,
                1.409012754e-004 };

        if (fabs(phase) > max_phi)
                return ERROR_VREGISTER(PUMAS_RETURN_VALUE_ERROR,
                    "magnetic rotation is too strong [%.5lE > %.5lE]",
                    fabs(phase), max_phi);

        /* Compute the local x, y and z. */
        double x1, x2, y1, y2, z1, z2;
        {
                /* Negative charge convention. */
                double a = -magnet * charge / density;
                x1 = x2 = y1 = y2 = 0.;
                double u = 1.;
                int j;
                for (j = 0; j < N_LARMOR_ORDERS; j++) {
                        x1 += u * poly_x[j] *
                            (*table_get_Li(physics, material, j, i1));
                        y1 += u * poly_y[j] *
                            (*table_get_Li(physics, material, j, i1));
                        x2 += u * poly_x[j] *
                            (*table_get_Li(physics, material, j, i2));
                        y2 += u * poly_y[j] *
                            (*table_get_Li(physics, material, j, i2));
                        u *= a;
                }
                z1 = *table_get_Li(physics, material, 0, i1);
                z2 = *table_get_Li(physics, material, 0, i2);
        }

        /* Update the position. */
        {
                double h = (k0 - *table_get_K(physics, i1)) /
                    (*table_get_K(physics, i2) - *table_get_K(physics, i1));
                *x = (x1 + h * (x2 - x1)) / density;
                *y = (y1 + h * (y2 - y1)) / density;
                *z = (z1 + h * (z2 - z1)) / density;
        }

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Propagation in arbitrary media.
 *
 * @param Physics         Handle for physics tables.
 * @param context         The simulation context.
 * @param state           The initial/final state.
 * @param medium_ptr      The initial/final propagation medium.
 * @param locals          Handle for the local properties of the starting
 *                        medium.
 * @param error           The error data.
 * @param step_max_medium The step limitation from the medium.
 * @param step_max_type   The type of geometry step (exact or approximate).
 * @param step_max_locals The step limitation from the local properties.
 * @return The end condition event.
 *
 * Transport through a set of media described by a medium callback. At output
 * the final kinetic energy, the total distance travelled and the total proper
 * time spent are updated.
 *
 * **Warning** : The initial state must have been initialized before the call.
 */
enum pumas_event transport_with_stepping(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium ** medium_ptr, struct medium_locals * locals,
    double step_max_medium, enum pumas_step step_max_type,
    double step_max_locals, struct error_context * error_)
{
        /* Check the config */
        if ((context->random == NULL) &&
            ((context->mode.energy_loss > PUMAS_MODE_CSDA) ||
                (context->mode.scattering == PUMAS_MODE_FULL_SPACE) ||
                (context->mode.decay == PUMAS_MODE_DECAY))) {
                ERROR_REGISTER(
                    PUMAS_RETURN_MISSING_RANDOM, "no random engine provided");
                return PUMAS_EVENT_NONE;
        }

        /* Unpack data. */
        struct pumas_medium * medium = *medium_ptr;
        int material = medium->material;

        /* Check for a straight path in a uniform medium */
        const enum pumas_mode scheme = context->mode.energy_loss;
        int straight =
            ((context->mode.scattering == PUMAS_MODE_LONGITUDINAL) &&
            (scheme <= PUMAS_MODE_HYBRID) &&
            (step_max_locals <= 0.) && !locals->magnetized) ?
            1 :
            0;

        /* Register the start of the the track, if recording. */
        struct simulation_context * const context_ =
            (struct simulation_context *)context;
        context_->step_event = PUMAS_EVENT_NONE;
        struct pumas_recorder * recorder = context->recorder;
        const int record = (recorder != NULL);
        if (record)
                record_state(context, medium,
                    context_->step_event | PUMAS_EVENT_START, state);

        /* Initialise some temporary data for the propagation, weights, ect ...
         */
        double ti = state->time;
        double wi = state->weight;
        double Xi = state->grammage;
        double dei, Xf;
        if (scheme > PUMAS_MODE_VIRTUAL) {
                const double ki = state->energy;
                Xf = cel_grammage(physics, context, scheme, material, ki);
                dei = 1. /
                    cel_energy_loss(physics, context, scheme, material, ki);

        } else {
                Xf = dei = 0.;
        }

        /* Check for any initial violation of external limits. */
        if ((context->event & PUMAS_EVENT_LIMIT_DISTANCE) &&
            (state->distance >= context->limit.distance))
                return PUMAS_EVENT_LIMIT_DISTANCE;
        if ((context->event & PUMAS_EVENT_LIMIT_GRAMMAGE) &&
            (state->grammage >= context->limit.grammage))
                return PUMAS_EVENT_LIMIT_GRAMMAGE;
        if ((context->event & PUMAS_EVENT_LIMIT_TIME) &&
            (state->time >= context->limit.time))
                return PUMAS_EVENT_LIMIT_TIME;
        if (state->weight <= 0.) return PUMAS_EVENT_WEIGHT;

        /* Initialise the stepping data. */
        double grammage_max;
        context_->step_event = PUMAS_EVENT_NONE;
        context_->step_first = 1;
        context_->step_X_limit = (context->event & PUMAS_EVENT_LIMIT_ENERGY) ?
            cel_grammage(physics, context,
                (scheme > PUMAS_MODE_VIRTUAL) ? scheme : PUMAS_MODE_CSDA,
                material, context->limit.energy) :
            0.;
        context_->step_invlb1 = 0;
        context_->step_rLarmor = 0.;
        memset(context_->step_uT, 0x0, 3 * sizeof(*(context_->step_uT)));
        transport_limit(
            physics, context, state, material, Xi, Xf, &grammage_max);
        if (context_->step_event) return context_->step_event;

        /* Step through the media. */
        int step_index = 1;
        for (;;) {
                /* Do a transportation step. */
                if (step_index > 1) {
                        /* Update the geometric step length */
                        step_max_type = context->medium(
                            context, state, NULL, &step_max_medium);
                        if ((step_max_medium > 0.) &&
                            (step_max_type == PUMAS_STEP_CHECK))
                                step_max_medium += 0.5 * STEP_MIN;
                }

                struct pumas_medium * new_medium = NULL;
                if (step_transport(physics, context, state, straight, medium,
                        locals, grammage_max, step_max_medium, step_max_type,
                        &step_max_locals, &new_medium, error_)
                        != PUMAS_RETURN_SUCCESS)
                        return context_->step_event;
                step_index++;

                /* Check for any event. */
                int record_step;
                if ((record != 0) && (recorder->period > 0) &&
                    (step_index != 1)) {
                        record_step = (step_index % recorder->period) == 0;
                } else {
                        record_step = 0;
                }
                if (!context_->step_event && !record_step) continue;

                /* Update the weight if a boundary or hard energy loss
                 * occured.
                 */
                if (context->mode.energy_loss >= PUMAS_MODE_CSDA) {
                        const double w0 =
                            (context->mode.decay == PUMAS_MODE_WEIGHT) ?
                            wi * exp(-fabs(state->time - ti) / physics->ctau) :
                            wi;
                        state->weight =
                            (context->mode.direction == PUMAS_MODE_FORWARD) ?
                            w0 :
                            w0 * cel_energy_loss(physics, context, scheme,
                                     material, state->energy) *
                                dei;
                } else {
                        state->weight =
                            wi * exp(-fabs(state->time - ti) / physics->ctau);
                }

                /* Process the event. */
                if (!context_->step_event) {
                        /* Register the current state. */
                        record_state(
                            context, medium, context_->step_event, state);
                } else {
                        if (context_->step_event &
                            (PUMAS_EVENT_LIMIT | PUMAS_EVENT_WEIGHT |
                                PUMAS_EVENT_VERTEX_DECAY)) {
                                /* A boundary was reached. Let's stop the
                                 * simulation.
                                 */
                                break;
                        } else if (context_->step_event &
                            (PUMAS_EVENT_VERTEX_DEL |
                                       PUMAS_EVENT_VERTEX_COULOMB)) {
                                /* A discrete process occured. */
                                if (context_->step_event &
                                    PUMAS_EVENT_VERTEX_DEL) {
                                        /* Backup the pre step point if
                                         * recording.
                                         */
                                        const int rec =
                                            record && (recorder->period > 0);
                                        double ki, ui[3];
                                        if (rec != 0) {
                                                ki = state->energy;
                                                memcpy(ui, state->direction,
                                                    sizeof(ui));
                                        }

                                        /* Apply the inelastic DEL. */
                                        transport_do_del(
                                            physics, context, state, material);

                                        /* Record the pre step point. */
                                        if (rec != 0) {
                                                const double kf =
                                                    state->energy;
                                                double uf[3];
                                                memcpy(uf, state->direction,
                                                    sizeof(uf));
                                                state->energy = ki;
                                                memcpy(state->direction, ui,
                                                    sizeof(state->direction));
                                                record_state(context, medium,
                                                    context_->step_event,
                                                    state);
                                                state->energy = kf;
                                                memcpy(state->direction, uf,
                                                    sizeof(state->direction));
                                        }

                                        /* Check for any stop condition */
                                        if (context->event &
                                            context_->step_event)
                                                break;

                                        /* Record the post step point. */
                                        if (rec != 0)
                                                record_state(context, medium,
                                                    PUMAS_EVENT_NONE, state);

                                        /* Reset the stepping data memory since
                                         * the kinetic energy has changed.
                                         */
                                        context_->step_first = 1;
                                        context_->step_invlb1 = 0;
                                        context_->step_rLarmor = 0.;
                                        memset(context_->step_uT, 0x0,
                                            3 * sizeof(*(context_->step_uT)));
                                } else {
                                        /* An EHS event occured. */
                                        transport_do_ehs(
                                            physics, context, state, material);

                                        /* Check for any stop condition */
                                        if (context->event &
                                            context_->step_event)
                                                break;
                                }

                                /* Update the locals if needed. */
                                if (step_max_locals > 0.) {
                                        context->medium(
                                            context, state, NULL, NULL);
                                        step_max_locals = transport_set_locals(
                                            context, medium, state, locals);
                                        if (locals->api.density <= 0.) {
                                                ERROR_REGISTER_NEGATIVE_DENSITY(
                                                    physics->material_name
                                                        [medium->material]);
                                                return context_->step_event;
                                        }
                                }
                        } else if (context_->step_event & PUMAS_EVENT_MEDIUM) {
                                /* A medium change occured. Let's update the
                                 * medium.
                                 */
                                medium = new_medium;
                                if ((medium == NULL) ||
                                    (context->event & PUMAS_EVENT_MEDIUM))
                                        break;
                                material = medium->material;
                                context->medium(context, state, NULL, NULL);
                                memset(&locals->api, 0x0, sizeof(locals->api));
                                locals->magnetized = 0;
                                step_max_locals = transport_set_locals(
                                    context, medium, state, locals);
                                if (locals->api.density <= 0.) {
                                        ERROR_REGISTER_NEGATIVE_DENSITY(
                                            physics->material_name
                                                [medium->material]);
                                        return context_->step_event;
                                }
                                straight = ((context->mode.scattering ==
                                    PUMAS_MODE_LONGITUDINAL) &&
                                    (scheme <= PUMAS_MODE_HYBRID) &&
                                    (step_max_locals <= 0.) &&
                                    !locals->magnetized) ?
                                    1 :
                                    0;

                                /* Update the kinetic limit converted
                                 * to grammage for this material.
                                 */
                                enum pumas_mode tmp_scheme =
                                    scheme > PUMAS_MODE_VIRTUAL ?
                                    scheme :
                                    PUMAS_MODE_CSDA;
                                context_->step_X_limit =
                                    (context->event &
                                        PUMAS_EVENT_LIMIT_ENERGY) ?
                                    cel_grammage(physics, context, tmp_scheme,
                                        material, context->limit.energy) :
                                    0.;

                                /* Reset the stepping data memory. */
                                context_->step_first = 1;
                                context_->step_invlb1 = 0;
                                context_->step_rLarmor = 0.;
                                memset(context_->step_uT, 0x0,
                                    3 * sizeof(*(context_->step_uT)));

                                /* Record the change of medium. */
                                if ((record) &&
                                    !(context->event & PUMAS_EVENT_MEDIUM))
                                        record_state(context, medium,
                                            context_->step_event, state);
                        } else {
                                /*  This should not happen. */
                                assert(0);
                        }

                        /* Update the initial conditions and the tracking of
                         * stepping events.
                         */
                        ti = state->time;
                        wi = state->weight;
                        Xi = state->grammage;
                        if (context->mode.energy_loss >= PUMAS_MODE_CSDA) {
                                const double ki = state->energy;
                                Xf = cel_grammage(
                                    physics, context, scheme, material, ki);
                                dei = 1. / cel_energy_loss(physics, context,
                                               scheme, material, ki);
                        }
                        transport_limit(physics, context, state, material, Xi,
                            Xf, &grammage_max);
                        if (context_->step_event) break;
                }
        }

        /* Protect final kinetic energy value against rounding errors. */
        if (context->mode.direction == PUMAS_MODE_FORWARD) {
                const double kinetic_min =
                    (context->event & PUMAS_EVENT_LIMIT_ENERGY) ?
                    context->limit.energy :
                    0.;
                if (fabs(state->energy - kinetic_min) < FLT_EPSILON)
                        state->energy = kinetic_min;
        } else {
                if ((context->event & PUMAS_EVENT_LIMIT_ENERGY) &&
                    (fabs(state->energy - context->limit.energy) <
                        FLT_EPSILON))
                        state->energy = context->limit.energy;
        }

        /* Register the end of the track, if recording. */
        if (record)
                record_state(context, medium,
                    context_->step_event | PUMAS_EVENT_STOP, state);

        *medium_ptr = medium;
        return context_->step_event;
}

/**
 * Set the local properties of a medium.
 *
 * @param context The simulation context.
 * @param state   The Monte-Carlo state.
 * @param locals  The local properties.
 * @return The proposed max step length is returned. A null or negative value
 * indicates a uniform and infinite medium.
 *
 * This is an encapsulation of the API locals callback where internal data are
 * also initialised.
 */
double transport_set_locals(const struct pumas_context * context,
    struct pumas_medium * medium, struct pumas_state * state,
    struct medium_locals * locals)
{
        struct pumas_locals * loc = (struct pumas_locals *)locals;
        if (medium->locals == NULL) {
                loc->density =
                    locals->physics->material_density[medium->material];
                memset(loc->magnet, 0x0, sizeof(loc->magnet));
                locals->magnetized = 0;

                return 0;
        } else {
                const double step_max = medium->locals(medium, state, loc);
                if (loc->density <= 0) {
                        loc->density =
                            locals->physics->material_density[medium->material];
                }

                const double * const b = loc->magnet;
                locals->magnetized =
                    ((b[0] != 0.) || (b[1] != 0.) || (b[2] != 0.)) ? 1 : 0;
                return step_max * context->accuracy;
        }
}

/**
 * Prepare the various limits for a MC propagation.
 *
 * @param Physics      Handle for physics tables.
 * @param context      The simulation context.
 * @param material     The index of the propagation material.
 * @param ki           The kinetic energy.
 * @param Xi           The total travelled grammage.
 * @param Xf           The total grammage over which the particle can travel.
 * @param grammage_max The total grammage at which a limit is reached.
 *
 * Compute the limit, translated into grammage, for Monte-Carlo steps.
 * At return *distance* is overwriten with the new limit and the context event
 * flags are updated.
 */
void transport_limit(const struct pumas_physics * physics,
    struct pumas_context * context, const struct pumas_state * state,
    int material, double Xi, double Xf, double * grammage_max)
{
        /* Initialise the stepping event flags. */
        struct simulation_context * const context_ =
            (struct simulation_context *)context;
        context_->step_event = context_->step_foreseen = PUMAS_EVENT_NONE;

        /* Check the Monte-Carlo weight. */
        if (state->weight <= 0.) {
                context_->step_event = PUMAS_EVENT_WEIGHT;
                return;
        }

        /* Check the kinetic limits. */
        if (context->mode.direction == PUMAS_MODE_FORWARD) {
                const double kinetic_min =
                    (context->event & PUMAS_EVENT_LIMIT_ENERGY) ?
                    context->limit.energy :
                    0.;
                if (state->energy <= kinetic_min) {
                        context_->step_event = PUMAS_EVENT_LIMIT_ENERGY;
                        return;
                }
        } else if ((context->event & PUMAS_EVENT_LIMIT_ENERGY) &&
            (state->energy >= context->limit.energy)) {
                context_->step_event = PUMAS_EVENT_LIMIT_ENERGY;
                return;
        };

        /* Initialise with the context grammage limit. */
        if (context->event & PUMAS_EVENT_LIMIT_GRAMMAGE) {
                *grammage_max = context->limit.grammage;
                context_->step_foreseen = PUMAS_EVENT_LIMIT_GRAMMAGE;

        } else {
                *grammage_max = 0.;
        }

        /* Check the NO LOSS case. */
        const enum pumas_mode scheme = context->mode.energy_loss;
        if (scheme == PUMAS_MODE_VIRTUAL) {
                if (context->mode.scattering == PUMAS_MODE_FULL_SPACE) {
                        const double X = Xi -
                            coulomb_ehs_length(
                                physics, context, material, state->energy) *
                                log(context->random(context));
                        if ((*grammage_max == 0.) || (X < *grammage_max)) {
                                *grammage_max = X;
                                context_->step_foreseen =
                                    PUMAS_EVENT_VERTEX_COULOMB;
                                return;
                        }
                }
                return;
        }

        /* Check for an inelastic DEL. */
        const double sgn =
            (context->mode.direction == PUMAS_MODE_FORWARD) ? 1. : -1.;
        enum pumas_event foreseen = PUMAS_EVENT_NONE;
        double kinetic_limit = 0.;
        if (scheme == PUMAS_MODE_HYBRID) {
                const double nI = del_interaction_length(physics, context,
                                      material, state->energy) +
                    sgn * log(context->random(context));
                if (nI > 0.) {
                        const double k = del_kinetic_from_interaction_length(
                            physics, context, material, nI);
                        if ((context->mode.direction == PUMAS_MODE_BACKWARD) ||
                            (k > *table_get_Kt(physics, material))) {
                                kinetic_limit = k;
                                foreseen = PUMAS_EVENT_VERTEX_DEL;
                        }
                }
        }

        /* Check for an EHS event. */
        if ((scheme < PUMAS_MODE_DETAILED) &&
            (context->mode.scattering == PUMAS_MODE_FULL_SPACE)) {
                const double nI = ehs_interaction_length(physics, context,
                                      scheme, material, state->energy) +
                    sgn * log(context->random(context));
                if (nI > 0.) {
                        const double k = ehs_kinetic_from_interaction_length(
                            physics, context, scheme, material, nI);
                        if ((kinetic_limit <= 0.) ||
                            ((context->mode.direction == PUMAS_MODE_FORWARD) &&
                             (k > kinetic_limit)) ||
                            ((context->mode.direction == PUMAS_MODE_BACKWARD) &&
                             (k < kinetic_limit))) {
                                kinetic_limit = k;
                                foreseen = PUMAS_EVENT_VERTEX_COULOMB;
                        }
                }
        }

        /* Return if no discrete event might occur. */
        if (foreseen == PUMAS_EVENT_NONE) return;

        /* Convert the kinetic limit to a grammage one and update. */
        const double X = Xi +
            sgn * (Xf - cel_grammage(
                            physics, context, scheme, material, kinetic_limit));
        if ((*grammage_max <= 0) || (X < *grammage_max)) {
                *grammage_max = X;
                context_->step_foreseen = foreseen;
        }
}

/**
 * Apply an inelastic DEL during the Monte-Carlo propagation.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param state    The initial/final state.
 */
void transport_do_del(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material)
{
        /* Update the energy. */
        double ki, kf;
        int process = -1;
        polar_function_t * polar_func;
        if (context->mode.direction == PUMAS_MODE_FORWARD) {
                ki = state->energy;
                polar_func = del_randomise_forward(
                    physics, context, state, material, &process);
                kf = state->energy;
        } else {
                kf = state->energy;
                polar_func = del_randomise_reverse(
                    physics, context, state, material, &process);
                ki = state->energy;
        }

        /* Update the event flag */
        struct simulation_context * context_ =
            (struct simulation_context *)context;
        if (process < 0) {
                context_->step_event = PUMAS_EVENT_NONE;
        } else {
                enum pumas_event event_for_[N_DEL_PROCESSES] = {
                        PUMAS_EVENT_VERTEX_BREMSSTRAHLUNG,
                        PUMAS_EVENT_VERTEX_PAIR_CREATION,
                        PUMAS_EVENT_VERTEX_PHOTONUCLEAR,
                        PUMAS_EVENT_VERTEX_DELTA_RAY
                };
                context_->step_event = event_for_[process];
        }

        /* Update the direction. */
        if ((context->mode.scattering == PUMAS_MODE_FULL_SPACE) &&
            (polar_func != NULL)) {
                const double ct = polar_func(physics, context, ki, kf);
                step_rotate_direction(context, state, ct);
        }
}

/*!
 * Apply an Elastic Hard Scattering (EHS) event during the Monte-Carlo
 * propagation.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param state    The initial/final state.
 */
void transport_do_ehs(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material)
{
        /* Unpack data. */
        struct simulation_context * ctx = (struct simulation_context *)context;
        struct coulomb_workspace * workspace = ctx->workspace;
        const double kinetic = state->energy;

        /* Get the cut-off angle and the soft scattering 1st moment. */
        double mu0 = 0., invlb1 = 0.;
        table_get_msc(physics, context, material, kinetic, &mu0, &invlb1);

        /* Compute the scattering parameters and the the total Coulomb
         * cross section.
         */
        int i;
        struct coulomb_data * data;
        double cs_tot = 0.;
        for (i = 0, data = workspace->data; i < physics->elements_in[material];
             i++, data++) {
                /* Compute the scattering parameters. */
                const struct material_component * const component =
                    &physics->composition[material][i];
                const struct atomic_element * const element =
                    physics->element[component->element];
                double kinetic0;
                coulomb_frame_parameters(
                    physics, kinetic, element->A, &kinetic0, data->fCM);
                data->fspin = coulomb_spin_factor(physics, kinetic0);
                coulomb_screening_parameters(physics, context, kinetic0,
                    component->element, data->screening);
                coulomb_pole_decomposition(data->screening, data->a, data->b);
                data->invlambda = component->fraction /
                    coulomb_wentzel_path(physics, kinetic0, element->Z,
                                      element->A, data->screening[0]);

                /* Compute the restricted Coulomb cross-section in the
                 * CM frame.
                 */
                data->cs_hard =
                    data->invlambda * coulomb_restricted_cs(mu0, data->fspin,
                                          data->screening, data->a, data->b);
                cs_tot += data->cs_hard;
        }

        /* Randomise the hard scatterer element. */
        const double cc0 = context->random(context) * cs_tot;
        int ihard = physics->elements_in[material] - 1;
        double cs = 0.;
        for (i = 0, data = workspace->data; i < physics->elements_in[material];
             i++, data++) {
                cs += data->cs_hard;
                if (cs >= cc0) {
                        ihard = i;
                        break;
                }
        }

        /* Compute the asymptotic hard angular parameter. */
        data = workspace->data + ihard;
        const double A = data->screening[0];
        const double mu2 = 1. - mu0;
        const double zeta = context->random(context);
        double mu1 = mu0 + (A + mu0) * zeta * mu2 / (1. + A - zeta * mu2);

        /* Configure for the root solver. */
        workspace->material = material;
        workspace->ihard = ihard;
        workspace->cs_h = data->cs_hard * (1. - zeta);

        /* Check the asymptotic value. */
        const double f1 =
            transport_hard_coulomb_objective(physics, mu1, workspace);
        if ((f1 < 0.) && (-f1 / workspace->cs_h >= 1E-06)) {
                /* We hit a form factor suppression. Let's call the
                 * root solver.
                 */
                const double f0 = data->cs_hard * zeta;
                math_find_root(transport_hard_coulomb_objective, physics, mu0,
                    mu1, &f0, &f1, 1E-06 * mu0, 1E-06, 100, workspace, &mu1);
        }

        /* Transform from the CM frame to the Lab frame. */
        const double gamma = data->fCM[0];
        const double tau = data->fCM[1];
        const double a = gamma * (tau + 1. - 2. * mu1);
        const double ct_h = a / sqrt(4. * mu1 * (1. - mu1) + a * a);

        /* Apply the rotation. */
        step_rotate_direction(context, state, ct_h);
}

/**
 * Objective function for solving the Coulomb scattering angle.
 *
 * @param Physics    Handle for physics tables.
 * @param mu         The proposed Coulomb scattering angular parameter.
 * @param parameters A pointer to the scattering workspace.
 * @return The difference between the current and expected restricted cross-
 * section.
 *
 * This is a wrapper for the root solver. It provides the objective function
 * to resolve.
 */
double transport_hard_coulomb_objective(
    const struct pumas_physics * physics, double mu, void * parameters)
{
        /* Unpack the workspace. */
        struct coulomb_workspace * workspace = parameters;
        struct coulomb_data * data = workspace->data + workspace->ihard;

        /* Compute the restricted cross section. */
        double cs_exp =
            data->invlambda * coulomb_restricted_cs(mu, data->fspin,
                                  data->screening, data->a, data->b);

        /* Return the difference with the expectation. */
        return cs_exp - workspace->cs_h;
}

/**
 * Randomise an inelastic DEL in forward MC.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param state    The initial/final state.
 * @param material The index of the propagation material.
 * @param process  The index opf the randomised process.
 * @return The polar function for the randomisation of the corresponding TT or
 * `NULL` if none.
 *
 * Below 10 GeV the DEL is randomised from a power law bias PDF. Above,
 * a ziggurat algorithm is used.
 */

polar_function_t * del_randomise_forward(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material,
    int * process)
{
        /* Check for a *do nothing* process. */
        if (state->energy <= *table_get_Kt(physics, material)) return NULL;

        /* Randomise the target element and the DEL sub-process. */
        struct del_info info;
        del_randomise_target(physics, context, state, material, &info);
        dcs_function_t * const dcs_func = dcs_get(info.process);
        const struct atomic_element * element = physics->element[info.element];

        if (info.process == 3) {
                const double m1 = physics->mass - ELECTRON_MASS;
                if (state->energy <= 0.5 * m1 * m1 / ELECTRON_MASS) {
                        state->energy = dcs_ionisation_randomise(physics,
                            context, element, state->energy, physics->cutoff);
                        *process = info.process;
                        return polar_get(info.process);
                }
        }

        /* Interpolate the lower fractional threshold. */
        int row;
        double xmin, hK;
        if (state->energy >= *table_get_K(physics, physics->n_energies - 1)) {
                row = physics->n_energies - 1;
                hK = 0.;
                xmin = *table_get_Xt(physics, info.process, info.element, row);
        } else {
                row = table_index(
                    physics, context, table_get_K(physics, 0), state->energy);
                if (row <= 0) {
                        row = 0;
                        hK = 0.;
                        xmin = *table_get_Xt(
                            physics, info.process, info.element, row);
                } else {
                        const double * const k = table_get_K(physics, row);
                        const double * const y = table_get_Xt(
                            physics, info.process, info.element, row);
                        hK = (state->energy - k[0]) / (k[1] - k[0]);
                        xmin = hK * y[1] + (1. - hK) * y[0];
                }
        }

        /* Set the upper fractionnal threshold. */
        double xmax = 1.;
        if (info.process == 3) {
                /* Ionisation case with upper kinematic threshold. */
                const double P2 =
                    state->energy * (state->energy + 2. * physics->mass);
                const double E = state->energy + physics->mass;
                const double Wmax =
                    2. * ELECTRON_MASS * P2 /
                    (physics->mass * physics->mass +
                        ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
                xmax = Wmax / state->energy;
                if (xmax > 1.) xmax = 1.;
        }
        if (xmin >= xmax) return NULL;

        if (state->energy <= 1.E+01) {
                /* Randomise the final energy with a bias model. */
                const double alpha[N_DEL_PROCESSES] = { 1.5, 3., 1.6, 2.3 };
                double r, w_bias;
                del_randomise_power_law(
                    context, alpha[info.process], xmin, xmax, &r, &w_bias);

                /* Update the kinetic energy and the Monte-Carlo weight. */
                const double d = dcs_evaluate(physics, context, dcs_func,
                    element, state->energy, state->energy * (1 - r));
                state->energy *= r;
                state->weight *= info.reverse.weight * d * w_bias;
        } else {
                /* Above 10 GeV rely on a direct sampling. */
                float tmp[DCS_SAMPLING_N], *dcs_samples = NULL;
                if ((info.process < N_DEL_PROCESSES - 1) &&
                    (state->energy >= DCS_MODEL_MIN_KINETIC)) {
                        /*  Interpolate the tabulated values and pass them
                         * to the sampling routine.
                         */
                        dcs_samples = tmp;
                        const int n = DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q +
                            1 + DCS_SAMPLING_N;
                        const float * y = table_get_dcs_value(physics, element,
                            info.process, row - physics->dcs_model_offset);
                        memcpy(tmp, y, DCS_SAMPLING_N * sizeof(float));
                        if (hK > 0.) {
                                y += n;
                                int i;
                                for (i = 0; i < DCS_SAMPLING_N; i++)
                                        tmp[i] = tmp[i] * (float)(1. - hK) +
                                            (float)(hK * y[i]);
                        }
                }
                del_randomise_ziggurat(physics, context, state, dcs_func,
                    element, xmin, xmax, dcs_samples);
        }

        *process = info.process;
        return polar_get(info.process);
}

/**
 * Randomise an inelastic DEL in backward MC.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param state    The initial/final state.
 * @param material The index of the propagation material.
 * @param process  The index of the randomised process.
 * @return The polar function for the randomisation of the corresponding TT or
 * `NULL` if none.
 *
 * A mixture PDF with a weight function is used for the randomisation
 * of the ancestor's state. First we check for a *do nothing* event, then the
 * DEL is processed by randomising the initial kinetic energy over a power law
 * distribution and then randomising the target element.
 */
polar_function_t * del_randomise_reverse(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material,
    int * process)
{
        /* Check for a pure CEL event. */
        const double lnq0 = -log(physics->cutoff);
        const double kt = *table_get_Kt(physics, material);
        const double kf = state->energy;
        const double pCEL =
            (kf >= kt) ? 0. : lnq0 / (lnq0 + log(kt / (kt - kf)));
        if (pCEL && (context->random(context) < pCEL)) {
                state->weight /= pCEL;
                return NULL;
        }
        const double xf =
            (kf >= kt * (1. - physics->cutoff)) ? physics->cutoff : 1. - kf / kt;

        /* Randomise the initial energy with a bias model. */
        double xmax, alpha;
        const double m1 = physics->mass - ELECTRON_MASS;
        if (state->energy < 0.5 * m1 * m1 / ELECTRON_MASS) {
                const double m2 = physics->mass + ELECTRON_MASS;
                xmax = 2. * ELECTRON_MASS *
                    (state->energy + 2. * physics->mass) / (m2 * m2);
                if (xmax < xf) return NULL;
                alpha = RMC_ALPHA_LOW;
        } else {
                alpha = RMC_ALPHA_HIGH;
                xmax = 1.;
        }

        double r, w_bias;
        del_randomise_power_law(context, alpha, xf, xmax, &r, &w_bias);
        w_bias /= r; /* Jacobian factor. */
        state->energy /= r;

        /* Randomise the target element and the SEL process. */
        struct del_info info;
        info.reverse.Q = state->energy - kf;
        del_randomise_target(physics, context, state, material, &info);
        dcs_function_t * const dcs_func = dcs_get(info.process);
        const struct atomic_element * element = physics->element[info.element];

        /* Update the MC weight. */
        const double f = dcs_evaluate(physics, context, dcs_func, element,
                             state->energy, state->energy - kf) *
            info.reverse.weight;
        state->weight *= w_bias * f *
            del_cross_section(physics, context, material, state->energy) /
            (del_cross_section(physics, context, material, kf) * (1. - pCEL));

        *process = info.process;
        return polar_get(info.process);
}

/**
 * Randomise the DEL over a generic power law PDF.
 *
 * @param context The simulation context.
 * @param alpha The power law exponent.
 * @param xmin The minimum fractional energy transfer.
 * @param xmax The maximum fractional energy transfer.
 * @param p_r The fractional final energy.
 * @param p_w The biasing weight factor.
 *
 * Randomise the initial energy with a power law bias model. Due to rounding
 * it can happen that the randomised energy is negative instead of very large.
 * If this occurs we try another step.
 */
void del_randomise_power_law(struct pumas_context * context, double alpha,
    double xmin, double xmax, double * p_r, double * p_w)
{
        double r = 0., w_bias;
        if (alpha == 1.) {
                const double lnq = log(xmax / xmin);
                for (;;) {
                        const double z = context->random(context);
                        r = xmin * exp(z * lnq);
                        if ((r < xmin) || (r >= xmax)) continue;
                        w_bias = lnq * r;
                        break;
                }
        } else {
                const double a1 = 1. - alpha;
                const double x0 = pow(xmin, a1);
                const double x1 = pow(xmax, a1);
                for (;;) {
                        const double z = context->random(context);
                        const double tmp = x0 + z * (x1 - x0);
                        r = pow(tmp, 1. / a1);
                        if ((r < xmin) || (r >= xmax)) continue;
                        w_bias = (x1 - x0) * r / (a1 * tmp);
                        break;
                }
        }

        *p_r = 1. - r;
        *p_w = w_bias;
}

/**
 * Randomise the DEL using a ziggurat like algorithm.
 *
 * @param Physics Handle for physics tables.
 * @param context The simulation context.
 * @param alpha   The power law exponent.
 * @param xmin    The minimum fractional energy transfer.
 * @param xmax    The maximum fractional energy transfer.
 *
 * Randomise the initial energy using a rejection sampling method. The DCS must
 * be a decreasing function of the energy loss. It is bounded by a piecewise
 * uniform pdf.
 */
void del_randomise_ziggurat(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state,
    dcs_function_t * dcs_func, const struct atomic_element * element,
    double xmin, double xmax, float * dcs_sampling)
{
        double x_b[DCS_SAMPLING_N + 1], dcs_b[DCS_SAMPLING_N];
        int n_b = DCS_SAMPLING_N, ib;
        if (dcs_sampling == NULL) {
                /* The DCS is uniformly sampled. */
                const double dx =
                    state->energy * (xmax - xmin) / DCS_SAMPLING_N;
                x_b[0] = xmin * state->energy;
                dcs_b[0] = dcs_evaluate(physics, context, dcs_func, element,
                    state->energy, x_b[0]);
                for (ib = 1; ib < DCS_SAMPLING_N; ib++) {
                        x_b[ib] = x_b[ib - 1] + dx;
                        dcs_b[ib] = dcs_evaluate(physics, context, dcs_func,
                            element, state->energy, x_b[ib]);
                        if ((dcs_b[ib] <= 0.) || (dcs_b[ib] > dcs_b[ib - 1])) {
                                /* Protect against numeric errors. */
                                x_b[ib] = xmax * state->energy;
                                dcs_b[ib] = 0.;
                                n_b = ib;
                                break;
                        }
                }
                if (n_b == DCS_SAMPLING_N)
                        x_b[DCS_SAMPLING_N] = x_b[DCS_SAMPLING_N - 1] + dx;
        } else {
                /* Use the pre-computed values if provided. */
                const double dnu = (1. - xmin) / DCS_SAMPLING_N;
                double nu;
                for (ib = 0, nu = xmin; ib < DCS_SAMPLING_N; ib++, nu += dnu) {
                        x_b[ib] = nu * state->energy;
                        dcs_b[ib] = (double)dcs_sampling[ib];
                }
                x_b[DCS_SAMPLING_N] =
                    x_b[DCS_SAMPLING_N - 1] + dnu * state->energy;
        }

        /* Build the CDF of the piecewise enveloppe. */
        double cdf_b[DCS_SAMPLING_N];
        cdf_b[0] = dcs_b[0] * (x_b[1] - x_b[0]);
        for (ib = 1; ib < n_b; ib++)
                cdf_b[ib] = cdf_b[ib - 1] + dcs_b[ib] * (x_b[ib + 1] - x_b[ib]);

        /* Randomise the energy transfer. */
        double xsol;
        for (;;) {
                const double zeta = context->random(context) * cdf_b[n_b - 1];
                int i;
                for (i = 0; i < n_b - 1; i++)
                        if (zeta <= cdf_b[i]) break;
                xsol =
                    (x_b[i + 1] - x_b[i]) * context->random(context) + x_b[i];
                const double dd = dcs_evaluate(
                    physics, context, dcs_func, element, state->energy, xsol);
                if (context->random(context) * dcs_b[i] <= dd) break;
        }

        /* Update the kinetic energy. */
        state->energy -= xsol;
}

/**
 * Randomise the target element and the sub-process for a DEL.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param state    The initial/final state.
 * @param  info    Temporary data for the DEL.
 *
 * In forward mode the target is selected according to the total cross-section.
 * In backward mode, since both the initial and final kinetic energy are known,
 * the target is selected according to the DCS. At output the *info* structure
 * is filled with the target element data and the MC weight, for backward mode.
 */
void del_randomise_target(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material,
    struct del_info * info)
{
        /* Interpolate the kinetic table. */
        int i1, i2;
        double h;
        i1 = table_index(
            physics, context, table_get_K(physics, 0), state->energy);
        if (i1 < 0) {
                i1 = i2 = 0;
                h = 0.;
        } else if (i1 >= physics->n_energies - 1) {
                i1 = i2 = physics->n_energies - 1;
                h = 0.;
        } else {
                i2 = i1 + 1;
                const double K1 = *table_get_K(physics, i1);
                h = (state->energy - K1) / (*table_get_K(physics, i2) - K1);
        }

        /* Randomise the target element and the DEL process. */
        const struct material_component * component;
        int ic, ic0 = 0, ip;
        for (ic = 0; ic < material; ic++) ic0 += physics->elements_in[ic];
        if (context->mode.direction == PUMAS_MODE_FORWARD) {
                /* Randomise according to the total cross section. */
                double zeta = context->random(context);
                component = physics->composition[material];
                for (ic = ic0; ic < ic0 + physics->elements_in[material];
                     ic++, component++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                const double f1 =
                                    *table_get_CSf(physics, ip, ic, i1);
                                const double f2 =
                                    *table_get_CSf(physics, ip, ic, i2);
                                const double csf = (f2 - f1) * h + f1;
                                if (!(zeta > csf)) {
                                        double csn;
                                        double csn1 = *table_get_CSn(physics,
                                            ip, component->element, i1);
                                        if (csn1 == 0.) {
                                                /* Linear interpolation. */
                                                const double csn2 =
                                                    *table_get_CSn(physics, ip,
                                                        component->element, i2);
                                                csn = (csn2 - csn1) * h + csn1;
                                        } else {
                                                /* Log interpolation. */
                                                csn1 = log(csn1);
                                                const double csn2 =
                                                    log(*table_get_CSn(physics,
                                                        ip, component->element,
                                                        i2));
                                                csn = exp(
                                                    (csn2 - csn1) * h + csn1);
                                        }
                                        info->reverse.weight = 1. / csn;
                                        goto target_found;
                                }
                        }
                assert(0); /* We should never reach this point ... */
        } else {
                /* Randomise according to the differential cross section. */
                double stot = 0.;
                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                        component = physics->composition[material];
                        for (ic = ic0;
                             ic < ic0 + physics->elements_in[material];
                             ic++, component++) {
                                const struct atomic_element * element =
                                    physics->element[component->element];
                                const double d = dcs_evaluate(physics, context,
                                    dcs_get(ip), element, state->energy,
                                    info->reverse.Q);
                                stot += d * component->fraction;
                        }
                }

                double zeta;
                for (;;) {
                        /* Prevent rounding errors. */
                        zeta = context->random(context) * stot;
                        if ((zeta > 0.) && (zeta < stot)) break;
                }
                double s = 0., csf_last = 0.;
                component = physics->composition[material];
                for (ic = ic0; ic < ic0 + physics->elements_in[material];
                     ic++, component++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                const int iel = component->element;
                                const struct atomic_element * element =
                                    physics->element[iel];
                                const double si =
                                    dcs_evaluate(physics, context, dcs_get(ip),
                                        element, state->energy,
                                        info->reverse.Q) *
                                    component->fraction;
                                s += si;
                                const double csf1 =
                                    *table_get_CSf(physics, ip, ic, i1);
                                const double csf2 =
                                    *table_get_CSf(physics, ip, ic, i2);
                                const double csf = (csf2 - csf1) * h + csf1;
                                if (!(zeta > s)) {
                                        double csn;
                                        double csn1 = *table_get_CSn(
                                            physics, ip, iel, i1);
                                        if (csn1 == 0.) {
                                                /* Linear interpolation. */
                                                const double csn2 =
                                                    *table_get_CSn(
                                                        physics, ip, iel, i2);
                                                csn = (csn2 - csn1) * h + csn1;
                                        } else {
                                                /* Log interpolation. */
                                                csn1 = log(csn1);
                                                const double csn2 =
                                                    log(*table_get_CSn(
                                                        physics, ip, iel, i2));
                                                csn = exp(
                                                    (csn2 - csn1) * h + csn1);
                                        }
                                        info->reverse.weight =
                                            (csf - csf_last) * stot /
                                            (si * csn);
                                        goto target_found;
                                }
                                csf_last = csf;
                        }
                assert(0); /* We should never reach this point ... */
        }

target_found:
        info->element = component->element;
        info->process = ip;
}

/*
 * Low level routines: MC stepping.
 */
/**
 * Perform a Monte-Carlo transport step in a single medium.
 *
 * @param Physics            Handle for physics tables.
 * @param context            The simulation context.
 * @param state              The initial/final state.
 * @param straight           Flag for a straight step.
 * @param medium_index       The index of the propagation medium.
 * @param locals             Handle for the local properties of the medium.
 * @param grammage_max       The maximum grammage until a limit is reached.
 * @param step_max_medium    The stepping limitation from medium boundaries.
 * @param step_max_type      The type of geometry step (approximate or exact).
 * @param step_max_locals    The stepping limitation from a non uniform medium.
 * @param out_index          The index of the end step medium or a negative
 *                           value if a boundary condition was reached.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise
 * `PUMAS_ERROR`.
 *
 * At return the state kinetic energy, position, direction and distance are
 * updated. In case of a stochastic CEL the proper time is also updated.
 */
enum pumas_return step_transport(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int straight,
    struct pumas_medium * medium, struct medium_locals * locals,
    double grammage_max, double step_max_medium,
    enum pumas_step step_max_type, double * step_max_locals,
    struct pumas_medium ** out_medium, struct error_context * error_)
{
        /* Unpack the data. */
        struct simulation_context * const context_ =
            (struct simulation_context *)context;
        const enum pumas_mode scheme = context->mode.energy_loss;
        double * const direction = state->direction;
        double * const position = state->position;
        const int material = medium->material;
        const double density = locals->api.density;
        const double density_i = 1. / density;
        const double momentum =
            sqrt(state->energy * (state->energy + 2. * physics->mass));

        /* Total grammage for the initial kinetic energy.  */
        const int tmp_scheme =
            (scheme == PUMAS_MODE_VIRTUAL) ? PUMAS_MODE_CSDA : scheme;
        const double Xtot = cel_grammage(
            physics, context, tmp_scheme, material, state->energy);

        /* Compute the local step length. */
        double step_loc, rLarmor = 0., uT[3];
        double invlb1 = 0.;
        if (!straight) {
                /* Compute the kinetic step length. */
                double r = context->accuracy;
                const double k_threshold = 1E+09;
                if ((scheme == PUMAS_MODE_DETAILED) &&
                    (state->energy > k_threshold)) {
                        /* In detailed mode, at very high energies shorter
                         * steps are needed.
                         */
                        double f = k_threshold / state->energy;
                        if (f < 0.1) f = 0.1;
                        r *= f;
                }
                step_loc = r * density_i * Xtot;

                if (context->mode.scattering == PUMAS_MODE_FULL_SPACE) {
                        /* Compute the soft scattering path length. */
                        if (context_->step_first != 0) {
                                double mu0;
                                table_get_msc(physics, context, material,
                                    state->energy, &mu0, &invlb1);
                                invlb1 *= density;
                        } else
                                invlb1 = context_->step_invlb1;
                        const double stepT = context->accuracy / invlb1;
                        if (stepT < step_loc) step_loc = stepT;
                }

                /* Compute the Larmor radius and the magnetic deflection
                 * direction.
                 */
                if (context_->step_first != 0) {
                        if (locals->magnetized == 0) {
                                rLarmor = 0.;
                        } else {
                                const double * const B = locals->api.magnet;
                                uT[0] =
                                    direction[1] * B[2] - direction[2] * B[1];
                                uT[1] =
                                    direction[2] * B[0] - direction[0] * B[2];
                                uT[2] =
                                    direction[0] * B[1] - direction[1] * B[0];
                                const double BT = sqrt(uT[0] * uT[0] +
                                    uT[1] * uT[1] + uT[2] * uT[2]);
                                if (BT == 0.) {
                                        rLarmor = 0.;
                                } else {
                                        const double iBT = 1. / BT;
                                        rLarmor =
                                            momentum * iBT / LARMOR_FACTOR;
                                        uT[0] *= iBT;
                                        uT[1] *= iBT;
                                        uT[2] *= iBT;
                                }
                        }
                } else {
                        rLarmor = context_->step_rLarmor;
                        memcpy(uT, context_->step_uT, 3 * sizeof(*uT));
                }

                /* Update the step length. */
                if (rLarmor > 0.) {
                        const double stepM = context->accuracy * rLarmor;
                        if (stepM < step_loc) step_loc = stepM;
                }
        } else {
                /* This is a straight step. */
                if (grammage_max <= 0.)
                        step_loc = 1E+09 * density_i;
                else
                        step_loc = grammage_max * density_i;
        }

        /* Update the `first step` flag. */
        context_->step_first = 0;

        /* Check the spatial resolution and the geometry step length. */
        if (step_loc < STEP_MIN) step_loc = STEP_MIN;
        double step = (step_max_medium <= 0) ?
            step_loc :
            (step_max_medium < step_loc ? step_max_medium : step_loc);
        if ((*step_max_locals > 0.) && (step > *step_max_locals))
                step = *step_max_locals;

        /* Check the total distance limitation. */
        *out_medium = medium;
        enum pumas_event event = PUMAS_EVENT_NONE;
        if (context->event & PUMAS_EVENT_LIMIT_DISTANCE) {
                const double d = context->limit.distance - state->distance;
                if (d <= step) {
                        step = d;
                        event = PUMAS_EVENT_LIMIT_DISTANCE;
                }
        }

        /* Update the position. */
        const double sgn =
            (context->mode.direction == PUMAS_MODE_FORWARD) ? 1. : -1.;
        position[0] += step * sgn * direction[0];
        position[1] += step * sgn * direction[1];
        position[2] += step * sgn * direction[2];

        /* Check for a change of medium. */
        struct pumas_medium * end_medium = NULL;
        context->medium(context, state, &end_medium, NULL);
        double end_position[3] = { position[0], position[1], position[2] };
        if (end_medium != medium) {
                if (step_max_type == PUMAS_STEP_CHECK) {
                        /* Check for an exact boundary. */
                        const double step_min =
                            (step < STEP_MIN) ? step : STEP_MIN;
                        double pi[3];
                        memcpy(pi, position, sizeof(pi));
                        double s1 = 0., s2 = -step_min;
                        position[0] = pi[0] + s2 * sgn * direction[0];
                        position[1] = pi[1] + s2 * sgn * direction[1];
                        position[2] = pi[2] + s2 * sgn * direction[2];
                        struct pumas_medium * tmp_medium = NULL;
                        context->medium(context, state, &tmp_medium, NULL);
                        if (tmp_medium != medium) {
                                /* Locate the medium change by dichotomy. */
                                if (tmp_medium != end_medium)
                                        end_medium = tmp_medium;
                                s1 = s2;
                                s2 = -step;
                                while (fabs(s1 - s2) > STEP_MIN) {
                                        double s3 = 0.5 * (s1 + s2);
                                        position[0] =
                                            pi[0] + s3 * sgn * direction[0];
                                        position[1] =
                                            pi[1] + s3 * sgn * direction[1];
                                        position[2] =
                                            pi[2] + s3 * sgn * direction[2];
                                        tmp_medium = NULL;
                                        context->medium(
                                            context, state, &tmp_medium, NULL);
                                        if (tmp_medium == medium) {
                                                s2 = s3;
                                        } else {
                                                s1 = s3;
                                                /* Update the end medium if
                                                 * required.
                                                 */
                                                if (tmp_medium != end_medium)
                                                        end_medium = tmp_medium;
                                        }
                                }
                                position[0] = pi[0] + s2 * sgn * direction[0];
                                position[1] = pi[1] + s2 * sgn * direction[1];
                                position[2] = pi[2] + s2 * sgn * direction[2];
                                step += s1;
                                end_position[0] =
                                    pi[0] + s1 * sgn * direction[0];
                                end_position[1] =
                                    pi[1] + s1 * sgn * direction[1];
                                end_position[2] =
                                    pi[2] + s1 * sgn * direction[2];

                                /* Force the last medium call to occur at the
                                 * final position. */
                                tmp_medium = NULL;
                                context->medium(
                                    context, state, &tmp_medium, NULL);
                        }
                }
                event = PUMAS_EVENT_MEDIUM;
                *out_medium = end_medium;
        }

        /*  Get the end step locals. */
        double Bi[3] = { 0., 0., 0. };
        if (locals->magnetized != 0) {
                Bi[0] = locals->api.magnet[0];
                Bi[1] = locals->api.magnet[1];
                Bi[2] = locals->api.magnet[2];
        }
        if ((*step_max_locals > 0.) && ((step_max_type != PUMAS_STEP_RAW) ||
            (event != PUMAS_EVENT_MEDIUM))) {
                /* Update the locals. */
                *step_max_locals = transport_set_locals(
                    context, medium, state, locals);
                if (locals->api.density <= 0.) {
                        ERROR_REGISTER_NEGATIVE_DENSITY(
                            physics->material_name[medium->material]);
                        return PUMAS_RETURN_DENSITY_ERROR;
                }
        }

        /* Offset the end step position for a boundary crossing. */
        if (event & PUMAS_EVENT_MEDIUM) {
                position[0] = end_position[0];
                position[1] = end_position[1];
                position[2] = end_position[2];
        }

        /* Set the end step kinetic energy. */
        double k1 = state->energy, dk = 0.;
        const double dX = 0.5 * step * (density + locals->api.density);
        if ((scheme >= PUMAS_MODE_CSDA) && (scheme <= PUMAS_MODE_HYBRID)) {
                /* Deterministic CEL with check for any kinetic limit. */
                const double X = Xtot - sgn * dX;
                if ((context->mode.direction == PUMAS_MODE_FORWARD) &&
                    (X <= context_->step_X_limit)) {
                        k1 = (context->event & PUMAS_EVENT_LIMIT_ENERGY) ?
                            context->limit.energy :
                            0.;
                        const double grammage =
                            state->grammage + Xtot - context_->step_X_limit;
                        if ((grammage_max <= 0.) || (grammage < grammage_max)) {
                                grammage_max = grammage;
                                event = PUMAS_EVENT_LIMIT_ENERGY;
                        }
                } else if ((context->mode.direction == PUMAS_MODE_BACKWARD) &&
                    (context_->step_X_limit > 0.) &&
                    (X > context_->step_X_limit)) {
                        k1 = (context->event & PUMAS_EVENT_LIMIT_ENERGY) ?
                            context->limit.energy :
                            0.;
                        const double grammage =
                            state->grammage + context_->step_X_limit - Xtot;
                        if ((grammage_max <= 0.) || (grammage < grammage_max)) {
                                grammage_max = grammage;
                                event = PUMAS_EVENT_LIMIT_ENERGY;
                        }
                } else
                        k1 = cel_kinetic_energy(
                            physics, context, scheme, material, X);
        } else if (scheme == PUMAS_MODE_DETAILED) {
                /* Fluctuate the CEL around its average value. */
                step_fluctuate(
                    physics, context, state, material, Xtot, dX, &k1, &dk);

                /* Check for a kinetic limit. */
                double kinetic_limit = -1.;
                if (context->mode.direction == PUMAS_MODE_FORWARD) {
                        const double kinetic_min =
                            (context->event & PUMAS_EVENT_LIMIT_ENERGY) ?
                            context->limit.energy :
                            0.;
                        if (k1 <= kinetic_min) kinetic_limit = kinetic_min;
                } else {
                        if ((context->event & PUMAS_EVENT_LIMIT_ENERGY) &&
                            (k1 >= context->limit.energy))
                                kinetic_limit = context->limit.energy;
                }
                if (kinetic_limit >= 0.) {
                        const double grammage = state->grammage +
                            fabs(state->energy - kinetic_limit) / dk * dX;
                        k1 = kinetic_limit;
                        if ((grammage_max <= 0.) || (grammage < grammage_max)) {
                                grammage_max = grammage;
                                event = PUMAS_EVENT_LIMIT_ENERGY;
                        }
                }

                /* Check for discrete events. */
                double Xmax = dX;
                if (grammage_max > 0.) {
                        const double dX1 = grammage_max - state->grammage;
                        if (Xmax > dX1) Xmax = dX1;
                }
                double k_x = -1.;

                /* Randomise an inelastic DEL. */
                double xs_del = del_cross_section(
                    physics, context, material, state->energy);
                const double tmp_del =
                    del_cross_section(physics, context, material, k1);
                if (tmp_del > xs_del) xs_del = tmp_del;
                if (xs_del <= 0.) goto no_del_event;
                const double X_del = -log(context->random(context)) / xs_del;
                if (X_del >= Xmax) goto no_del_event;
                const double k_del = state->energy - sgn * X_del * dk / dX;
                if (k_del <= 0.) goto no_del_event;
                const double r =
                    del_cross_section(physics, context, material, k_del) /
                    xs_del;
                if (context->random(context) > r) goto no_del_event;
                k_x = k_del;
                Xmax = X_del;
                event = PUMAS_EVENT_VERTEX_DEL;
        no_del_event:

                if (context->mode.scattering == PUMAS_MODE_FULL_SPACE) {
                        /* Randomise an EHS. */
                        double kmin, kmax;
                        if (k1 <= state->energy)
                                kmin = k1, kmax = state->energy;
                        else
                                kmax = k1, kmin = state->energy;
                        if (kmin < physics->table_K[1]) {
                                const double tmp = kmin;
                                kmin = kmax, kmax = tmp;
                        }
                        double lb_ehs = coulomb_ehs_length(
                            physics, context, material, kmin);
                        if (lb_ehs <= 0.) goto no_ehs_event;
                        const double X_ehs =
                            -log(context->random(context)) * lb_ehs;
                        if (X_ehs >= Xmax) goto no_ehs_event;
                        const double k_ehs =
                            state->energy - sgn * X_ehs * dk / dX;
                        if (k_ehs <= 0.) goto no_ehs_event;
                        const double r = coulomb_ehs_length(physics, context,
                                             material, k_ehs) /
                            lb_ehs;
                        if ((r <= 0.) || (context->random(context) > 1. / r))
                                goto no_ehs_event;
                        k_x = k_ehs;
                        Xmax = X_ehs;
                        event = PUMAS_EVENT_VERTEX_COULOMB;
                no_ehs_event:;
                }

                /* Apply the discrete event, if any. */
                if (k_x >= 0.) {
                        k1 = k_x;
                        grammage_max = state->grammage + Xmax;
                }
        }

        /* Check and update the grammage. */
        const double Xi = state->grammage;
        state->grammage += dX;
        const double sf0 = step;
        double h_int = 0.;
        if ((grammage_max > 0.) && (state->grammage >= grammage_max)) {
                const double dX_ = grammage_max - Xi;
                if ((fabs(density - locals->api.density) <= FLT_EPSILON) ||
                    (sf0 <= FLT_EPSILON)) {
                        step = dX_ * density_i;
                        h_int = dX_ / (state->grammage - Xi);
                } else {
                        const double drho = locals->api.density - density;
                        double tmp =
                            density * density + 2. * dX_ * drho / sf0;
                        tmp = (tmp > FLT_EPSILON) ? sqrt(tmp) : 0;
                        h_int = (tmp - density) / drho;
                        step *= h_int;
                }
                state->grammage = grammage_max;
                if (!(event & PUMAS_EVENT_LIMIT_ENERGY)) {
                        /*  Update the kinetic energy. */
                        if (scheme <= PUMAS_MODE_HYBRID) {
                                if (scheme != PUMAS_MODE_VIRTUAL)
                                        k1 = cel_kinetic_energy(physics,
                                            context, scheme, material, Xtot -
                                                sgn * (state->grammage - Xi));
                                event = context_->step_foreseen;
                        } else if (!(event & (PUMAS_EVENT_VERTEX_COULOMB |
                                                 PUMAS_EVENT_VERTEX_DEL))) {
                                /*
                                 * This is a grammage limit in the detailed
                                 * scheme.
                                 */
                                k1 = state->energy -
                                    sgn * (state->grammage - Xi) * dk / dX;
                                if (k1 < 0.) k1 = 0.;
                                event = PUMAS_EVENT_LIMIT_GRAMMAGE;
                        }
                }

                /* Update the position. */
                const double ds_ = step - sf0;
                position[0] += ds_ * sgn * direction[0];
                position[1] += ds_ * sgn * direction[1];
                position[2] += ds_ * sgn * direction[2];
        }

        /* Check the proper time limit. */
        int decayed = 0;
        double time_max = (context->event & PUMAS_EVENT_LIMIT_TIME) ?
            context->limit.time : 0.;
        if (context->mode.decay == PUMAS_MODE_DECAY) {
                if ((time_max <= 0.) || (context_->lifetime < time_max)) {
                        time_max = context_->lifetime;
                        decayed = 1;
                }
        }

        const double sf1 = step;
        if (straight && (scheme != PUMAS_MODE_VIRTUAL)) {
                const double Ti = cel_proper_time(
                    physics, context, scheme, material, state->energy);
                if (time_max > 0.) {
                        const double Tf =
                            Ti - sgn * (time_max - state->time) * density;
                        if (Tf > 0.) {
                                const double dxT = fabs(
                                    Xtot - cel_grammage_as_time(physics,
                                               context, scheme, material, Tf));
                                if (Xi + dxT < state->grammage) {
                                        /* A proper time limit is reached. */
                                        event = PUMAS_EVENT_LIMIT_TIME;
                                        state->time = time_max;
                                        state->decayed = decayed;
                                        step = dxT * density_i;
                                        if (step > sf1) step = sf1;
                                        state->grammage = Xi + dxT;
                                        const double xf = Xtot - sgn * dxT;
                                        k1 = cel_kinetic_energy(physics,
                                            context, scheme, material, xf);
                                }
                        }
                }

                if (event != PUMAS_EVENT_LIMIT_TIME) {
                        const double Tf = cel_proper_time(
                            physics, context, scheme, material, k1);
                        state->time += fabs(Tf - Ti) * density_i;
                }
        } else {
                const double p_f =
                    (k1 <= 0) ? momentum : sqrt(k1 * (k1 + 2. * physics->mass));
                const double ti = state->time;
                state->time +=
                    0.5 * step * physics->mass * (1. / momentum + 1. / p_f);
                if ((time_max > 0.) && (state->time >= time_max)) {
                        /* A proper time limit is reached. */
                        event = PUMAS_EVENT_LIMIT_TIME;
                        state->time = time_max;
                        state->decayed = decayed;

                        /* Interpolate the step length. */
                        const double a = (momentum / p_f - 1.) / sf1;
                        const double c =
                            -2. * (time_max - ti) * momentum / physics->mass;
                        if (a != 0.) {
                                double delta = 1. - a * c;
                                if (delta < 0.) delta = 0.;
                                step = 1. / a * (sqrt(delta) - 1.);
                        } else {
                                step = -0.5 * c;
                        }
                        if (step < 0.) step = 0.;

                        /*  Update the kinetic energy. */
                        if (scheme != PUMAS_MODE_VIRTUAL) {
                                const double p1 = momentum / (1 + a * step);
                                k1 = sqrt(p1 * p1 -
                                         physics->mass * physics->mass) -
                                    physics->mass;
                                if (k1 < 0.) k1 = 0.;
                        }

                        /* Correct the grammage. */
                        const double Xf = Xi + dX;
                        if (density == locals->api.density)
                                state->grammage = Xi + step * density;
                        else
                                state->grammage = Xi +
                                    step * (density +
                                               0.5 * step / sf0 *
                                                   (locals->api.density -
                                                       density));
                        if (state->grammage > Xf) state->grammage = Xf;
                }
        }

        if (event == PUMAS_EVENT_LIMIT_TIME) {
                /* Correct the position. */
                const double ds_ = step - sf1;
                position[0] += ds_ * sgn * direction[0];
                position[1] += ds_ * sgn * direction[1];
                position[2] += ds_ * sgn * direction[2];
        }

        /* Update the event flag. */
        context_->step_event = event;

        /* Update the kinetic energy. */
        state->energy = k1;

        /* Update the travelled distance. */
        state->distance += step;

        /* Compute the multiple scattering path length. */
        if (context->mode.scattering == PUMAS_MODE_FULL_SPACE) {
                double mu0, invlb1_;
                if (state->energy <= 0.) {
                        context_->step_invlb1 = invlb1;
                } else {
                        table_get_msc(physics, context, material,
                            state->energy, &mu0, &invlb1_);
                        context_->step_invlb1 = locals->api.density * invlb1_;
                }
        }

        /* Update the end step properties. */
        if (locals->magnetized) {
                /* Interpolate the end point magnetic field if needed. */
                double B[3] = { locals->api.magnet[0], locals->api.magnet[1],
                        locals->api.magnet[2] };
                int magnetized = locals->magnetized;
                if (h_int > 0.) {
                        B[0] = Bi[0] + h_int * (B[0] - Bi[0]);
                        B[1] = Bi[1] + h_int * (B[1] - Bi[1]);
                        B[2] = Bi[2] + h_int * (B[2] - Bi[2]);
                        magnetized =
                            (B[0] * B[0] + B[1] * B[1] + B[2] * B[2] == 0.);
                }

                /* Compute the Larmor radius and the magnetic deflection
                 * direction.
                 */
                if (magnetized == 0) {
                        context_->step_rLarmor = 0.;
                        memset(context_->step_uT, 0x0, sizeof(uT));
                } else {
                        double * const uT_ = context_->step_uT;
                        uT_[0] = direction[1] * B[2] - direction[2] * B[1];
                        uT_[1] = direction[2] * B[0] - direction[0] * B[2];
                        uT_[2] = direction[0] * B[1] - direction[1] * B[0];
                        const double BT = sqrt(uT_[0] * uT_[0] +
                            uT_[1] * uT_[1] + uT_[2] * uT_[2]);
                        if (BT == 0.) {
                                rLarmor = 0.;
                        } else {
                                const double p = (state->energy <= 0) ?
                                    momentum :
                                    sqrt(state->energy *
                                        (state->energy + 2. * physics->mass));
                                const double iBT = 1. / BT;
                                context_->step_rLarmor =
                                    p * iBT / LARMOR_FACTOR;
                                uT_[0] *= iBT;
                                uT_[1] *= iBT;
                                uT_[2] *= iBT;
                        }
                }
        } else {
                context_->step_rLarmor = rLarmor;
                memcpy(context_->step_uT, uT, sizeof(uT));
        }

        /* Apply the magnetic deflection. */
        if ((rLarmor > 0.) || (context_->step_rLarmor > 0.)) {
                const double u[3] = { direction[0], direction[1],
                        direction[2] };
                double theta_i =
                    (rLarmor > 0.) ? state->charge * step / rLarmor : 0.;
                double theta_f = (context_->step_rLarmor > 0.) ?
                    state->charge * step / context_->step_rLarmor :
                    0.;
                if (context->mode.direction == PUMAS_MODE_BACKWARD) {
                        theta_i = -theta_i;
                        theta_f = -theta_f;
                }
                const double ci = cos(theta_i);
                const double si = sin(theta_i);
                const double cf = cos(theta_f);
                const double sf = sin(theta_f);
                const double c = 0.5 * (ci + cf);
                direction[0] =
                    c * u[0] + 0.5 * (si * uT[0] + sf * context_->step_uT[0]);
                direction[1] =
                    c * u[1] + 0.5 * (si * uT[1] + sf * context_->step_uT[1]);
                direction[2] =
                    c * u[2] + 0.5 * (si * uT[2] + sf * context_->step_uT[2]);
                const double norm = 1. / sqrt(direction[0] * direction[0] +
                                             direction[1] * direction[1] +
                                             direction[2] * direction[2]);
                direction[0] *= norm;
                direction[1] *= norm;
                direction[2] *= norm;
        }

        /* Apply the multiple scattering. */
        if ((invlb1 > 0.) || (context_->step_invlb1 > 0.)) {
                double ilb1 = 0.5 * step * (invlb1 + context_->step_invlb1);
                if (ilb1 > 1.) ilb1 = 1.;
                double ct;
                do
                        ct = 1. + ilb1 * log(context->random(context));
                while (ct < -1.);
                step_rotate_direction(context, state, ct);
        }

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Apply a stochastic CEL on a MC step.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param state    The particle Monte-Carlo state.
 * @param material The target material.
 * @param Xtot     The initial CSDA total grammage.
 * @param dX       The step grammage length.
 * @param kf       The expected final state kinetic energy.
 * @param dE       The total energy loss for path dX.
 *
 * PENELOPE's model is used for the soft energy loss randomisation. Note that
 * `dE` differs from `ki-kf` if the kinetic energy falls to `0` before the
 * end of the step.
 */
static void step_fluctuate(const struct pumas_physics * physics,
    struct pumas_context * context, struct pumas_state * state, int material,
    double Xtot, double dX, double * kf, double * dE)
{
        const enum pumas_mode scheme = context->mode.energy_loss;
        const double sgn =
            (context->mode.direction == PUMAS_MODE_FORWARD) ? 1. : -1.;
        double k1, dk = 0.;
        k1 = cel_kinetic_energy(
            physics, context, scheme, material, Xtot - sgn * dX);
        if (k1 > 0.) {
                const double dk1 =
                    sqrt(0.5 * dX *
                        (step_fluctuations2(physics, material, state->energy) +
                             step_fluctuations2(physics, material, k1)));
                const double dk0 = fabs(state->energy - k1);
                if (dk0 >= 3. * dk1) {
                        double u;
                        do
                                u = step_randn(context);
                        while (fabs(u) > 3.);
                        u /= 1.015387;
                        k1 += u * dk1;
                } else if (dk0 >= 1.7320508 * dk1) {
                        const double u =
                            1.7320508 * (1. - 2. * context->random(context));
                        k1 += u * dk1;
                } else {
                        const double dk32 = 3. * dk1 * dk1;
                        const double dk02 = dk0 * dk0;
                        const double a =
                            1. - (dk32 - dk02) / (dk32 + 3. * dk02);
                        if (context->random(context) <= a) {
                                const double b = 0.5 * (dk32 + 3. * dk02) / dk0;
                                const double u = context->random(context);
                                k1 = state->energy - sgn * b * u;
                        } else {
                                k1 = state->energy;
                        }
                }
                dk = fabs(state->energy - k1);
        }
        if (dk == 0.)
                dk = cel_energy_loss(
                         physics, context, scheme, material, state->energy) *
                    dX;

        /* Copy back the result. */
        *kf = k1;
        *dE = dk;
}

/**
 * Squared standard deviation of the stochastic CEL.
 *
 * @param Physics  Handle for physics tables.
 * @param material The target material.
 * @param kinetic  The projectile kinetic energy.
 * @return The squared standard deviation.
 *
 * The Gaussian thick aborber approximation of ICRU Report 49 is assumed.
 */
static double step_fluctuations2(
    const struct pumas_physics * physics, int material, double kinetic)
{
        const double r = ELECTRON_MASS / physics->mass;
        const double g = 1. + kinetic / physics->mass;
        const double b2 = 1. - 1. / (g * g);
        double qmax = 2. * ELECTRON_MASS * b2 * g * g / (1. + r * (2. * g + r));
        const double qcut = physics->cutoff * kinetic;
        if (qmax > qcut) qmax = qcut;
        return 1.535375E-05 * physics->material_ZoA[material] *
            (1. / b2 - 0.5) * qmax;
}

/**
 * Gaussian normal random number.
 *
 * @param context The simulation context.
 * @return a random number distributed according to a normal distribution.
 *
 * The Box-Muller algorithm is used. The random variates are generated in pairs.
 * The *context* is used as local storage.
 */
static double step_randn(struct pumas_context * context)
{
        struct simulation_context * const context_ =
            (struct simulation_context *)context;
        context_->randn_done = !context_->randn_done;
        if (!context_->randn_done) return context_->randn_next;

        const double r = sqrt(-2. * log(context->random(context)));
        const double phi = 2. * M_PI * context->random(context);
        const double c = cos(phi);
        const double s = sin(phi);
        context_->randn_next = r * c;
        return r * s;
}

/**
 * Rotate the direction randomly but constraining the polar angle.
 *
 * @param context   The simulation context.
 * @param state     The initial/final state.
 * @param cos_theta The cosine of the polar rotation angle.
 *
 * The direction is randomly rotated around the initial direction with the
 * constraint that the cosine of the angle between both directions is
 * *cos_theta*.
 */
void step_rotate_direction(struct pumas_context * context,
    struct pumas_state * state, double cos_theta)
{
        /* Unpack data. */
        double * const direction = state->direction;

        /* Check the numerical sine. */
        const double stsq = 1. - cos_theta * cos_theta;
        if (stsq <= 0.) return;
        const double st = sqrt(stsq);

        /* select the co-vectors for the local basis. */
        double u0x = 0., u0y = 0., u0z = 0.;
        const double a0 = fabs(direction[0]);
        const double a1 = fabs(direction[1]);
        const double a2 = fabs(direction[2]);
        if (a0 > a1) {
                if (a0 > a2) {
                        const double nrm =
                            1. / sqrt(direction[0] * direction[0] +
                                     direction[2] * direction[2]);
                        u0x = -direction[2] * nrm, u0z = direction[0] * nrm;
                } else {
                        const double nrm =
                            1. / sqrt(direction[1] * direction[1] +
                                     direction[2] * direction[2]);
                        u0y = direction[2] * nrm, u0z = -direction[1] * nrm;
                }
        } else {
                if (a1 > a2) {
                        const double nrm =
                            1. / sqrt(direction[0] * direction[0] +
                                     direction[1] * direction[1]);
                        u0x = direction[1] * nrm, u0y = -direction[0] * nrm;
                } else {
                        const double nrm =
                            1. / sqrt(direction[1] * direction[1] +
                                     direction[2] * direction[2]);
                        u0y = direction[2] * nrm, u0z = -direction[1] * nrm;
                }
        }
        const double u1x = u0y * direction[2] - u0z * direction[1];
        const double u1y = u0z * direction[0] - u0x * direction[2];
        const double u1z = u0x * direction[1] - u0y * direction[0];

        /* Apply the rotation. */
        const double phi = M_PI * (1. - 2. * context->random(context));
        const double cp = cos(phi);
        const double sp = sin(phi);
        direction[0] = cos_theta * direction[0] + st * (cp * u0x + sp * u1x);
        direction[1] = cos_theta * direction[1] + st * (cp * u0y + sp * u1y);
        direction[2] = cos_theta * direction[2] + st * (cp * u0z + sp * u1z);
}

/*
 * Low level routine: indexing of atomic elements.
 */
/**
 * Find the table index of an element given its name.
 *
 * @param Physics  Handle for physics tables.
 * @param name     The element name.
 * @return The element index or `-1` if it wasn't found.
 */
int element_index(const struct pumas_physics * physics, const char * name)
{
        int index = 0;
        for (; index < physics->n_elements; index++)
                if (strcmp(physics->element[index]->name, name) == 0)
                        return index;
        return -1;
}

/*
 * Low level routine: indexing of materials.
 */
/**
 * Find the table index of a material given its name.
 *
 * @param Physics  Handle for physics tables.
 * @param name     The material name.
 * @param error_   The error data.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned.
 */
static enum pumas_return material_index(const struct pumas_physics * physics,
    const char * material, int * index, struct error_context * error_)
{
        int i;
        for (i = 0; i < physics->n_materials; i++) {
                if (strcmp(physics->material_name[i], material) == 0) {
                        *index = i;
                        return PUMAS_RETURN_SUCCESS;
                }
        }

        return ERROR_VREGISTER(
            PUMAS_RETURN_UNKNOWN_MATERIAL, "unknown material `%s'", material);
}

/*
 * Low level routines: various functions related to Coulomb scattering.
 */
/**
 * Compute the atomic and nuclear screening parameters.
 *
 * @param Physics   Handle for physics tables.
 * @param context   The simulation context.
 * @param kinetic   The kinetic energy.
 * @param element   The scatterer atomic element.
 * @param screening The computed screening parameters.
 *
 * The atomic screening parameter is computed according to Kuraev's
 * parameterisation for relativistic particles or using Moliere's one at low
 * energies.
 */
void coulomb_screening_parameters(const struct pumas_physics * physics,
    struct pumas_context * context, double kinetic, int element,
    double * screening)
{
        /* Nuclear screening. */
        const double third = 1. / 3;
        const double A13 = pow(physics->element[element]->A, third);
        const double R1 = 1.02934 * A13 + 0.435;
        const double R2 = 2.;
        const double p2 = kinetic * (kinetic + 2. * physics->mass);
        const double d = 5.8406E-02 / p2;
        screening[1] = d / (R1 * R1);
        screening[2] = d / (R2 * R2);

        /* Atomic Moliere screening with Coulomb correction from Kuraev et al.
         * Phys. Rev. D 89, 116016 (2014). Valid for ultra-relativistic
         * particles only.
         */
        struct atomic_element * e = physics->element[element];
        const double etot = kinetic + physics->mass;
        const double ZE = e->Z * etot;
        const double zeta2 = 5.3251346E-05 * (ZE * ZE) / p2;
        double cK;
        if (zeta2 > 1.) {
                /* Let's perform the serie computation. */
                int i, n = 10 + (int)(e->Z);
                double f = 0.;
                for (i = 1; i <= n; i++) f += zeta2 / (i * (i * i + zeta2));
                cK = exp(f);
        } else {
                /* Let's use Kuraev's approximate expression. */
                cK = exp(1. - 1. / (1. + zeta2) +
                    zeta2 * (0.2021 + zeta2 * (0.0083 * zeta2 - 0.0369)));
        }

        /* Original Moliere's atomic screening, considered as a reference
         * value at low energies.
         */
        const double cM = 1. + 3.34 * zeta2;

        /* Atomic screening interpolation. */
        double r = kinetic / etot;
        r *= r;
        const double c = r * cK + (1. - r) * cM;
        screening[0] = 5.179587126E-12 * pow(e->Z, 2. / 3.) * c / p2;
}

/**
 * Compute the Coulomb mean free path assuming a Wentzel DCS with only atomic
 * screening.
 *
 * @param Physics   Handle for physics tables.
 * @param kinetic   The kinetic energy.
 * @param Z         The atomic number of the target element.
 * @param A         The atomic mass of the target element.
 * @param screening The atomic screening factor.
 * @return The mean free path in kg/m^2.
 */
double coulomb_wentzel_path(const struct pumas_physics * physics,
    double kinetic, double Z, double A, double screening)
{
        const double d = kinetic * (kinetic + 2. * physics->mass) /
            (Z * (kinetic + physics->mass));
        return A * 2.54910918E+08 * screening * (1. + screening) * d * d;
}

/**
 * Encapsulation of the interaction length for EHS.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param material The propagation material.
 * @param kinetic  The kinetic energy.
 * @return The grammage length for hard scattering events.
 */
double coulomb_ehs_length(const struct pumas_physics * physics,
    struct pumas_context * context, int material, double kinetic)
{
        const double p2 = kinetic * (kinetic + 2. * physics->mass);
        const int imax = physics->n_energies - 1;
        if (kinetic < *table_get_K(physics, 1)) {
                return *table_get_Lb(physics, material, 1) / p2;
        } else if (kinetic >= *table_get_K(physics, imax)) {
                return *table_get_Lb(physics, material, imax) / p2;
        } else {
                const int i1 = table_index(
                    physics, context, table_get_K(physics, 0), kinetic);
                const int i2 = i1 + 1;
                double h = (kinetic - *table_get_K(physics, i1)) /
                    (*table_get_K(physics, i2) - *table_get_K(physics, i1));
                return (*table_get_Lb(physics, material, i1) +
                           h * (*table_get_Lb(physics, material, i2) -
                                   *table_get_Lb(physics, material, i1))) /
                    p2;
        }
}

/**
 * Compute the spin factor to the Coulomb DCS.
 *
 * @param Physics Handle for physics tables.
 * @param kinetic The kinetic energy.
 * @return The spin factor.
 */
double coulomb_spin_factor(const struct pumas_physics * physics, double kinetic)
{
        const double e = kinetic + physics->mass;
        return kinetic * (e + physics->mass) / (e * e);
}

/**
 * Compute the parameters of the CM to the Lab frame transform for the
 * Coulomb DCS.
 *
 * @param Physics    Handle for physics tables.
 * @param kinetic    The kinetic energy.
 * @param Ma         The target mass, in atomic unit.
 * @param kinetic0   The CM kinetic energy.
 * @param parameters The vector of Lorentz parameters (gamma, tau).
 */
void coulomb_frame_parameters(const struct pumas_physics * physics,
    double kinetic, double Ma, double * kinetic0, double * parameters)
{
        Ma *= 0.931494; /* Convert from atomic unit to GeV. */
        double M2 = physics->mass + Ma;
        M2 *= M2;
        const double sCM12i = 1. / sqrt(M2 + 2. * Ma * kinetic);
        parameters[0] = (kinetic + physics->mass + Ma) * sCM12i;
        *kinetic0 =
            (kinetic * Ma + physics->mass * (physics->mass + Ma)) * sCM12i -
            physics->mass;
        if (*kinetic0 < 1E-09) *kinetic0 = 1E-09;
        const double etot = kinetic + physics->mass + Ma;
        const double betaCM2 =
            kinetic * (kinetic + 2. * physics->mass) / (etot * etot);
        double rM2 = physics->mass / Ma;
        rM2 *= rM2;
        parameters[1] = sqrt(rM2 * (1. - betaCM2) + betaCM2);
}

/**
 * Compute the coefficients of the Coulomb DCS pole reduction.
 *
 * @param screening The screening parameters: A, N1 and N2.
 * @param a         The vector of coefficients for 1st order poles.
 * @param b         The vector of coefficients for 2nd order poles.
 *
 * The Coulomb DCS goes as 1/(A+mu)^2*1/(N1+mu)^2*1/(N2+mu)^2, with A the atomic
 * screening angle and Ni the nuclear ones. It is reduced to a sum over poles
 * as: a0/(A+mu)+b0/(A+mu)^2+...+a2/(N2+mu)+b2/(N2+mu)^2.
 */
void coulomb_pole_decomposition(double * screening, double * a, double * b)
{
        const double d01 = 1. / (screening[0] - screening[1]);
        const double d02 = 1. / (screening[0] - screening[2]);
        const double d12 = 1. / (screening[1] - screening[2]);
        b[0] = d01 * d01 * d02 * d02;
        b[1] = d01 * d01 * d12 * d12;
        b[2] = d12 * d12 * d02 * d02;
        a[0] = 2. * b[0] * (d01 + d02);
        a[1] = 2. * b[1] * (d12 - d01);
        a[2] = -2. * b[2] * (d12 + d02);
}

/**
 * Compute the restricted EHS cross-section in the CM frame.
 *
 * @param mu        The angular cut-off value.
 * @param fspin     The spin correction factors.
 * @param screening The screening parameters.
 * @param a         The 1st order pole coefficients.
 * @param b         The 2nd order pole coefficients.
 *
 * The restricted cross-section is integrated from *mu* to `1`.
 */
double coulomb_restricted_cs(
    double mu, double fspin, double * screening, double * a, double * b)
{
        if (mu >= 1.) return 0.;

        const double nuclear_screening =
            (screening[1] < screening[2]) ? screening[1] : screening[2];
        if (mu < 1E-08 * nuclear_screening) {
                /* We neglect the nucleus finite size. */
                const double L = log((screening[0] + 1.) / (screening[0] + mu));
                const double r =
                    (1. - mu) / ((screening[0] + mu) * (screening[0] + 1.));
                const double k = screening[0] * (1. + screening[0]);
                return k * (r - fspin * (L - screening[0] * r));
        } else {
                /* We need to take all factors into account using a pole
                 * reduction. */
                double I0[3], I1[3], J0[3], J1[3];
                int i;
                for (i = 0; i < 3; i++) {
                        const double L =
                            log((screening[i] + 1.) / (screening[i] + mu));
                        const double r = (1. - mu) /
                            ((screening[i] + mu) * (screening[i] + 1.));
                        I0[i] = r;
                        J0[i] = L;
                        I1[i] = L - screening[i] * r;
                        J1[i] = mu - screening[i] * L;
                }

                const double k = screening[0] * (1. + screening[0]) *
                    screening[1] * screening[1] * screening[2] * screening[2];
                double cs = 0.;
                for (i = 0; i < 3; i++) {
                        cs += a[i] * (J0[i] - fspin * J1[i]) +
                            b[i] * (I0[i] - fspin * I1[i]);
                }
                return k * cs;
        }
}

/**
 * Compute the order 0 and 1 transport coefficients.
 *
 * @param mu          The angular cut-off.
 * @param fspin       The spin factor.
 * @param screening   The screening factors.
 * @param a           The poles 1st order coefficients.
 * @param b           The poles 2nd order coefficients.
 * @param coefficient The computed transport coefficients.
 */
void coulomb_transport_coefficients(double mu, double fspin, double * screening,
    double * a, double * b, double * coefficient)
{
        const double nuclear_screening =
            (screening[1] < screening[2]) ? screening[1] : screening[2];
        if (mu < 1E-08 * nuclear_screening) {
                /* We neglect the nucleus finite size. */
                const double L = log(1. + mu / screening[0]);
                const double r = mu / (mu + screening[0]);
                const double k = screening[0] * (1. + screening[0]);
                coefficient[0] = k * (r / screening[0] - fspin * (L - r));
                const double I2 = mu - screening[0] * (r - 2. * L);
                coefficient[1] = 2. * k * (L - r - fspin * I2);
        } else {
                /* We need to take all factors into account using a pole
                 * reduction. */
                double I0[3], I1[3], I2[3], J0[3], J1[3], J2[3];
                int i;
                double mu2 = 0.5 * mu * mu;
                for (i = 0; i < 3; i++) {
                        double r = mu / (mu + screening[i]);
                        double L = log(1. + mu / screening[i]);
                        double mu1 = mu;
                        I0[i] = r / screening[i];
                        J0[i] = L;
                        I1[i] = L - r;
                        r *= screening[i];
                        L *= screening[i];
                        J1[i] = mu1 - L;
                        I2[i] = mu1 - 2. * L + r;
                        L *= screening[i];
                        mu1 *= screening[i];
                        J2[i] = mu2 + L - mu1;
                }

                const double k = screening[0] * (1. + screening[0]) *
                    screening[1] * screening[1] * screening[2] * screening[2];
                coefficient[0] = coefficient[1] = 0.;
                for (i = 0; i < 3; i++) {
                        coefficient[0] += a[i] * (J0[i] - fspin * J1[i]) +
                            b[i] * (I0[i] - fspin * I1[i]);
                        coefficient[1] += a[i] * (J1[i] - fspin * J2[i]) +
                            b[i] * (I1[i] - fspin * I2[i]);
                }
                coefficient[0] *= k;
                coefficient[1] *= 2. * k;
        }
}

/**
 * The 1st transport coefficient for multiple scattering on electronic
 * shells.
 *
 * @param Physics Handle for physics tables.
 * @param element The target atomic element.
 * @param kinetic The projectile initiale kinetic energy.
 * @return The inverse of the electronic 1st transport path length in kg/m^2.
 *
 * The contribution from atomic electronic shells is computed following
 * Salvat et al., NIMB316 (2013) 144-159, considering only close interactions
 * and approximating the electronic structure by a single shell of energy I
 * with occupancy Z.
 */
double transverse_transport_ionisation(const struct pumas_physics * physics,
    const struct atomic_element * element, double kinetic)
{
        /* Soft close interactions, restricted to physics->cutoff. */
        const double momentum2 = kinetic * (kinetic + 2. * physics->mass);
        const double E = kinetic + physics->mass;
        const double Wmax = 2. * ELECTRON_MASS * momentum2 /
            (physics->mass * physics->mass +
                                ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
        const double W0 = 2. * momentum2 / ELECTRON_MASS;
        const double mu_max = Wmax / W0;
        double mu3 = kinetic * physics->cutoff / W0;
        if (mu3 > mu_max) mu3 = mu_max;
        const double mu2 = 0.62 * element->I / W0;
        if (mu2 >= mu3) return 0.;
        const double a0 = 0.5 * W0 / momentum2;
        const double a1 = -1. / Wmax;
        const double a2 = E * E / (W0 * momentum2);
        const double cs0 = 1.535336E-05 / element->A; /* m^2/kg/GeV. */
        return 2. * cs0 * element->Z *
            (0.5 * a0 * (mu3 * mu3 - mu2 * mu2) + a1 * (mu3 - mu2) +
                   a2 * log(mu3 / mu2));
}

/**
 * The 1st transport coefficient for multiple scattering in soft photonuclear
 * events.
 *
 * @param Physics Handle for physics tables.
 * @param element The target atomic element.
 * @param kinetic The projectile initial kinetic energy.
 * @return The inverse of the photonuclear 1st transport path length in kg/m^2.
 *
 * The doubly differential cross section is computed assuming a Q^2 dependency
 * as given by GEANT4 Physics Reference Manual and further approximating 1-nu
 * as 1. Following: Q^2 = mu*p^2, and the integration over mu can be done
 * analyticaly. The integration over nu is done numericaly with a Gaussian
 * quadrature.
 */
double transverse_transport_photonuclear(const struct pumas_physics * physics,
    const struct atomic_element * element, double kinetic)
{
        /* Integration over the kinetic transfer, q, done with a log sampling.
         */
        const double E = kinetic + physics->mass;
        double x0 = log(1E-06);
        double x1 = 0.;
        math_gauss_quad(100, &x0, &x1); /* Initialisation. */

        double xi, wi, lbipn = 0.;
        while (math_gauss_quad(0, &xi, &wi) == 0) { /* Stepping. */
                const double nu = physics->cutoff * exp(xi);
                const double q = nu * kinetic;

                /* Analytical integration over mu. */
                const double m02 = 0.4;
                const double q2 = q * q;
                const double tmax = 1.876544 * q;
                const double tmin =
                    q2 * physics->mass * physics->mass / (E * (E - q));
                const double b1 = 1. / (1. - q2 / m02);
                const double c1 = 1. / (1. - m02 / q2);
                double L1 = b1 * log((q2 + tmax) / (q2 + tmin));
                double L2 = c1 * log((m02 + tmax) / (m02 + tmin));
                const double I0 = log(tmax / tmin) - L1 - L2;
                L1 *= q2;
                L2 *= m02;
                const double I1 = L1 + L2;
                L1 *= q2;
                L2 *= m02;
                const double I2 =
                    (tmax - tmin) * (b1 * q2 + c1 * m02) - L1 - L2;
                const double ratio =
                    (I1 * tmax - I2) / ((I0 * tmax - I1) * kinetic *
                                           (kinetic + 2. * physics->mass));

                /* Update the double integral value.  */
                lbipn += ratio * nu *
                    dcs_photonuclear(physics, element, kinetic, nu * kinetic) *
                    wi;
        }
        return 2. * lbipn;
}

/* Low level routine: helper function for recording a MC state. */
/**
 * Register a Monte-Carlo state.
 *
 * @param recorder     The recorder handle.
 * @param medium       The medium in which the particle is located.
 * @param event        The current stepping event.
 * @param state        The Monte-Carlo state to record.
 *
 * This routine adds the given state to the recorder's stack.
 */
void record_state(struct pumas_context * context, struct pumas_medium * medium,
    enum pumas_event event, struct pumas_state * state)
{
        /* Check for a user supplied recorder */
        struct pumas_recorder * recorder = context->recorder;
        if (recorder->record != NULL) {
                recorder->record(context, state, medium, event);
                return;
        }

        /* Do the default recording ... */
        struct frame_recorder * const rec =
            (struct frame_recorder * const)recorder;
        struct frame_stack * stack = rec->stack;
        struct pumas_frame * frame = NULL;

        if ((stack == NULL) || (stack->size < sizeof(*frame))) {
                /* Allocate a new memory segment. */
                const int size = 4096;
                stack = allocate(size);
                if (stack == NULL) return;
                stack->size = size - sizeof(*stack);
                stack->frame = stack->frames;
                stack->next = rec->stack;
                rec->stack = stack;
        }

        /* Allocate the new frame in the stack. */
        frame = stack->frame++;
        stack->size -= sizeof(*frame);

        /* Link the new frame. */
        if (recorder->first == NULL)
                recorder->first = frame;
        else
                rec->last->next = frame;
        rec->last = frame;
        frame->next = NULL;
        frame->medium = medium;
        frame->event = event;
        memcpy(&(frame->state), state, sizeof(*state));
        recorder->length++;
}

/* Low level routine: utility function for memory alignment. */
/**
 * Compute the padded memory size.
 *
 * @param size     The requested memory size.
 * @param pad_size The memory padding size.
 * @return The padded memory size.
 *
 * The padded memory size is the smallest integer multiple of *pad_size* and
 * greater or equal to *size*. It allows to align memory addresses on multiples
 * of pad_size.
 */
int memory_padded_size(int size, int pad_size)
{
        int i = size / pad_size;
        if ((size % pad_size) != 0) i++;
        return i * pad_size;
}

/* Low level routine: error handling. */

/**
 * Utility function for formating errors.
 *
 * @param error_     The error data.
 * @param rc         The return code.
 * @param caller     The calling function from which to return.
 * @param file       The file where the error occured.
 * @param line       The line where the error occured.
 * @param message    A brief description of the error.
 * @return The return code is forwarded.
 */
static enum pumas_return error_format(struct error_context * error_,
    enum pumas_return rc, const char * file, int line, const char * format, ...)
{
        if (error_ == NULL) return rc;

        error_->code = rc;
        if ((s_error.handler == NULL) || (rc == PUMAS_RETURN_SUCCESS))
                return rc;

        /* Format the error message */
        const int n =
            snprintf(error_->message, ERROR_MSG_LENGTH, "{ %s [#%d], %s:%d } ",
                pumas_error_function(error_->function), rc, file, line);
        if (n < ERROR_MSG_LENGTH - 1) {
                va_list ap;
                va_start(ap, format);
                vsnprintf(
                    error_->message + n, ERROR_MSG_LENGTH - n, format, ap);
                va_end(ap);
        }

        return rc;
}

/**
 * Utility function for handling errors.
 *
 * @param ezrror_ The error data.
 * @return The return code is forwarded.
 */
static enum pumas_return error_raise(struct error_context * error_)
{
        if ((s_error.handler == NULL) || (error_->code == PUMAS_RETURN_SUCCESS))
                return error_->code;

        if (s_error.catch) {
                if (s_error.catch_error.code == PUMAS_RETURN_SUCCESS) {
                        memcpy(&s_error.catch_error, error_,
                            sizeof(s_error.catch_error));
                }
                return error_->code;
        }
        s_error.handler(error_->code, error_->function, error_->message);

        return error_->code;
}

/*
 * Low level routines: I/O and parsing.
 */
/**
 * Parse a dE/dX file.
 *
 * @param Physics     Handle for physics tables.
 * @param fid         The file handle.
 * @param material    The material index.
 * @param filename    The name of the curent file.
 * @param error_      The error data.
 * @return On succces `PUMAS_RETURN_SUCCESS`, otherwise `PUMAS_ERROR`.
 *
 * Parse a dE/dX data table in PDG text file format.
 */
enum pumas_return io_parse_dedx_file(struct pumas_physics * physics, FILE * fid,
    int material, const char * filename, struct error_context * error_)
{
        char * buffer = NULL;

        /* Parse the Stenheimer coefficients from the header lines. */
        int line = 0, read_coef = 0;
        int i;
        for (i = 0; i < physics->n_energy_loss_header; i++) {
                io_read_line(fid, &buffer, filename, line, error_);
                line++;
                if (error_->code != PUMAS_RETURN_SUCCESS) return error_->code;

                if (read_coef) {
                        double * const d =
                            (double * const)(
                                &physics->material_density_effect[material]);
                        double * const I =
                            (double * const)(
                                &physics->material_I[material]);
                        if (sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf",
                            d, d + 1, d + 2, d + 3, I, d + 4, d + 5) != 7) {
                                return ERROR_REGISTER(
                                    PUMAS_RETURN_FORMAT_ERROR,
                                    "could not read Sternheimer coef");
                        }
                        physics->material_I[material] *= 1E-09;
                        read_coef = 0;
                } else if (strstr(buffer, "Sternheimer coef") != NULL) {
                        read_coef = 1;
                }
        }

        /* Initialise the new table. */
        int row = 0;
        *table_get_T(physics, PUMAS_MODE_CSDA, material, row) = 0.;
        *table_get_T(physics, PUMAS_MODE_HYBRID, material, row) = 0.;
        *table_get_K(physics, row) = 0.;
        *table_get_dE(physics, PUMAS_MODE_CSDA, material, row) = 0.;
        *table_get_dE(physics, PUMAS_MODE_HYBRID, material, row) = 0.;
        *table_get_NI_el(physics, PUMAS_MODE_CSDA, material, row) = 0.;
        *table_get_NI_el(physics, PUMAS_MODE_HYBRID, material, row) = 0.;
        *table_get_NI_in(physics, material, row) = 0.;
        *table_get_CS(physics, material, row) = 0.;
        int ip;
        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                *table_get_CSf(physics, ip, material, row) = 0.;
                *table_get_CSn(physics, ip, material, row) = 0.;
        }
        row++;

        /* Scan the new table. */
        while (error_->code == PUMAS_RETURN_SUCCESS) {
                io_read_line(fid, &buffer, filename, line, error_);
                line++;
                if (error_->code != PUMAS_RETURN_SUCCESS) break;
                io_parse_dedx_row(
                    physics, buffer, material, &row, filename, line, error_);
        }

        if (error_->code != PUMAS_RETURN_SUCCESS) {
                if ((error_->code == PUMAS_RETURN_END_OF_FILE) || feof(fid))
                        error_->code = PUMAS_RETURN_SUCCESS;
        }
        if (error_->code != PUMAS_RETURN_SUCCESS) return error_->code;

        if (row != physics->n_energies)
                return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                    "inconsistent number of rows in energy loss tables [%s, %d "
                    "!= %d]",
                    filename, row, physics->n_energies);

        compute_regularise_del(physics, material);

        return error_->code;
}

/**
 * Parse a row of a dE/dX file.
 *
 * @param Physics  Handle for physics tables.
 * @param buffer   The read buffer containing the row data.
 * @param material The index of the material.
 * @param row      The index of the parsed row.
 * @param filename The name of the current file.
 * @param line     The line being processed.
 * @param error_   The error data.
 * @return On succees `PUMAS_RETURN_SUCCESS`, or `PUMAS_ERROR` otherwise.
 *
 * Parse a row of a dE/dX data table formated in PDG text file format.
 */
enum pumas_return io_parse_dedx_row(struct pumas_physics * physics,
    char * buffer, int material, int * row, const char * filename, int line,
    struct error_context * error_)
{
        /*
         * Skip the peculiar values since they differ from material to
         * material.
         */
        if ((strstr(buffer, "Minimum ionization") != NULL) ||
            (strstr(buffer, "critical energy") != NULL))
                return PUMAS_RETURN_SUCCESS;

        /* parse the new data line */
        double k, a, be, de, brems, pair, photo;
        int count = sscanf(buffer, "%lf %*f %lf %lf %lf %lf %lf %lf", &k, &a,
            &brems, &pair, &photo, &be, &de);
        if ((count != 7) || (*row >= physics->n_energies))
                return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                    "invalid data line [@%s:%d]", filename, line);

        /* Change units. MeV to GeV and MeV cm^2/g to GeV m^2/kg. */
        k *= 1E-03;
        a *= 1E-04;
        be *= 1E-04;
        de *= 1E-04;
        brems *= 1E-04;
        pair *= 1E-04;
        photo *= 1E-04;

        /* Check the consistency of kinetic values. */
        if (material == 0) {
                *table_get_K(physics, *row) = k;
        } else if (*table_get_K(physics, *row) != k)
                return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                    "inconsistent kinetic energy value [@%s:%d, %.5lE != "
                    "%.5lE]",
                    filename, line, *table_get_K(physics, *row), k);

        /* Compute the fractional contributions of the energy loss
         * processes to the CEL and to DELs.
         */
        static double * cel_table = NULL;
        if (material == 0) {
                /* Precompute the per element terms. */
                if ((cel_table = compute_cel_and_del(physics, *row)) == NULL)
                        return ERROR_REGISTER_MEMORY();
        }

        struct material_component * component = physics->composition[material];
        double frct_cel[] = { 0., 0., 0., 0. };
        double frct_cs[] = { 0., 0., 0., 0. };
        int ic, ic0 = 0;
        for (ic = 0; ic < material; ic++) ic0 += physics->elements_in[ic];
        for (ic = ic0; ic < ic0 + physics->elements_in[material];
             ic++, component++) {
                int iel = component->element;
                const double w = component->fraction;
                /* Loop over processes. */
                int ip;
                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                        const double f =
                            (*table_get_CSn(physics, ip, iel, *row)) * w;
                        *table_get_CSf(physics, ip, ic, *row) = f;
                        frct_cs[ip] += f;
                        frct_cel[ip] +=
                            *table_get_cel(physics, ip, iel, *row, cel_table) *
                            w;
                }
        }

        const double cel_max[] = { brems, pair, photo, a };
        double be_cel = 0., frct_cs_del = 0.;
        int ip;
        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                be_cel +=
                    (frct_cel[ip] < cel_max[ip]) ? frct_cel[ip] : cel_max[ip];
                frct_cs_del += frct_cs[ip];
        }

        /* Normalise the fractional cross-sections terms. */
        if (frct_cs_del <= 0.) {
                for (ic = ic0; ic < ic0 + physics->elements_in[material]; ic++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++)
                                *table_get_CSf(physics, ip, ic, *row) = 0.;
        } else {
                double sum_tot = 0.;
                for (ic = ic0; ic < ic0 + physics->elements_in[material]; ic++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++)
                                sum_tot +=
                                    *table_get_CSf(physics, ip, ic, *row);
                double sum = 0.;
                for (ic = ic0; ic < ic0 + physics->elements_in[material]; ic++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                sum += *table_get_CSf(physics, ip, ic, *row);
                                *table_get_CSf(physics, ip, ic, *row) =
                                    sum / sum_tot;
                        }

                /* Protect against rounding errors. */
                ic = ic0 + physics->elements_in[material] - 1;
                ip = N_DEL_PROCESSES - 1;
                *table_get_CSf(physics, ip, ic, *row) = 1.;
        }

        /* Update the table values */
        const double de_cel = de - be_cel;

        /* If this is the last entry, save the energy loss values. */
        if (*row == physics->n_energies - 1) {
                const double etot = k + physics->mass;
                *table_get_a_max(physics, material) = a;
                *table_get_b_max(physics, PUMAS_MODE_CSDA, material) =
                    be / etot;
                *table_get_b_max(physics, PUMAS_MODE_HYBRID, material) =
                    (be - be_cel) / etot;
        }

        /* End point statistics */
        *table_get_dE(physics, PUMAS_MODE_CSDA, material, *row) = de;
        *table_get_dE(physics, PUMAS_MODE_HYBRID, material, *row) = de_cel;
        *table_get_CS(physics, material, *row) = frct_cs_del;

        /* Weighted integrands */
        const double dei = 1. / de_cel;
        *table_get_X(physics, PUMAS_MODE_HYBRID, material, *row) = dei;
        *table_get_NI_in(physics, material, *row) = frct_cs_del * dei;

        (*row)++;
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Read a line from a file.
 *
 * @param fid      The file handle.
 * @param buf      A pointer to the new line or `NULL` in case of faillure.
 * @param filename The file currently processed.
 * @param line     The current line in the file.
 * @param error_   The error data.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise an error code.
 *
 * This routine manages a dynamic buffer. If *fid* is not `NULL`, this routine
 * reads a new line and fills in a pointer to a string holding the read data.
 * Otherwise, if *fid* is `NULL` any memory previously allocated by the routine
 * is released.
 */
enum pumas_return io_read_line(FILE * fid, char ** buf, const char * filename,
    int line, struct error_context * error_)
{
        static int size = 0;
        static char * buffer = NULL;

        if (buf != NULL) *buf = NULL;
        if (fid == NULL) {
                /* Release the buffer memory. */
                deallocate(buffer);
                buffer = NULL;
                size = 0;
                return PUMAS_RETURN_SUCCESS;
        }

        if (buffer == NULL) {
                /* Allocate the buffer if not already done. */
                size = 2048;
                buffer = allocate(size * sizeof(*buffer));
                if (buffer == NULL) return ERROR_REGISTER_MEMORY();
        }

        /* Get a full line. */
        char * s = buffer;
        int to_read = size;
        for (;;) {
                const char endline1 = '\n';
                const char endline2 = '\r';
                const char * r = fgets(s, to_read, fid);
                if (r == NULL) return ERROR_REGISTER_EOF(filename);
                int n_read = strlen(s);
                if ((n_read >= to_read - 1) && (s[to_read - 2] != endline1) &&
                    (s[to_read - 2] != endline2)) {
                        size += 2048;
                        char * new_buffer =
                            reallocate(buffer, size * sizeof(*buffer));
                        if (new_buffer == NULL) {
                                return ERROR_REGISTER_MEMORY();
                        }
                        buffer = new_buffer;
                        s += to_read - 1;
                        to_read = 2049;
                        continue;
                }
                break;
        }

        *buf = buffer;
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Parse the global settings from a MDF.
 *
 * @param Physics      Handle for physics tables.
 * @param mdf          The MDF handle.
 * @param dedx_path    The path to the energy loss table(s).
 * @param error_       The error data.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise an error code.
 */
enum pumas_return mdf_parse_settings(const struct pumas_physics * physics,
    struct mdf_buffer * mdf, const char * dedx_path,
    struct error_context * error_)
{
        /* Initialisation of settings. */
        mdf->n_energies = 0;
        mdf->dcs_model_offset = 0;
        mdf->n_energy_loss_header = -1;
        mdf->n_materials = 0;
        mdf->n_composites = 0;
        mdf->n_elements = 0;
        mdf->n_components = 0;
        mdf->max_components = 0;
        mdf->size_composite = 0;
        mdf->size_dedx_path = 0;
        mdf->size_elements_names = 0;
        mdf->size_materials_names = 0;

        /* Prepare the path to the first dedx file. */
        char * full_path = NULL;
        char * filename = NULL;
        int offset_dir, size_path, size_name = 0;
        if (!mdf->dry_mode) {
                mdf->size_dedx_path = strlen(dedx_path) + 1;
                if (mdf_format_path(dedx_path, mdf->mdf_path, &full_path,
                        &offset_dir, &size_path,
                        error_) != PUMAS_RETURN_SUCCESS)
                        goto clean_and_exit;
        }

        /* Prepare the index tables. */
        if (mdf_settings_index(MDF_INDEX_INITIALISE, 0, error_) !=
            PUMAS_RETURN_SUCCESS)
                goto clean_and_exit;

        /* Loop on the XML nodes. */
        const int pad_size = sizeof(*(physics->data));
        mdf->left = 0;
        mdf->line = 1;
        mdf->depth = MDF_DEPTH_EXTERN;
        struct mdf_node node;
        for (;;) {
                /* Get the next node. */
                mdf_get_node(mdf, &node, error_);
                /* Check the termination. */
                if (error_->code == PUMAS_RETURN_END_OF_FILE) {
                        if (mdf->depth == MDF_DEPTH_EXTERN) {
                                /* This is a normal terminations. */
                                error_->code = PUMAS_RETURN_SUCCESS;
                                break;
                        } else {
                                goto clean_and_exit;
                        }
                } else if (error_->code != PUMAS_RETURN_SUCCESS)
                        goto clean_and_exit;
                if (node.key == MDF_KEY_OTHER) continue;

                /* Process material nodes. */
                if (node.key == MDF_KEY_MATERIAL) {
                        if (node.head == 1) {
                                /* This a material head. */
                                mdf->n_materials++;
                                const int size = strlen(node.at1.name) + 1;
                                mdf->size_materials_names += size;

                                /* Copy the filename if first. */
                                if (filename == NULL) {
                                        size_name = strlen(node.at2.file);
                                        filename = allocate(size_name + 1);
                                        if (filename == NULL) {
                                                ERROR_REGISTER_MEMORY();
                                                goto clean_and_exit;
                                        }
                                        strcpy(filename, node.at2.file);
                                }

                                /* Update the names table and book a new index
                                 * table.
                                 */
                                if ((mdf_settings_name(size, 'M', node.at1.name,
                                         error_) != PUMAS_RETURN_SUCCESS) ||
                                    (mdf_settings_index(
                                         MDF_INDEX_WRITE_MATERIAL, 0, error_) !=
                                        PUMAS_RETURN_SUCCESS))
                                        goto clean_and_exit;
                        } else {
                                if (mdf->elements_in > mdf->max_components)
                                        mdf->max_components = mdf->elements_in;
                                /* Finalise the index table. */
                                if (mdf_settings_index(
                                        MDF_INDEX_FINALISE_MATERIAL,
                                        mdf->elements_in,
                                        error_) != PUMAS_RETURN_SUCCESS)
                                        goto clean_and_exit;
                        }
                }

                /* Process composite material nodes. */
                if (node.key == MDF_KEY_COMPOSITE) {
                        if (node.head == 1) {
                                /* This a material head. */
                                mdf->n_materials++;
                                mdf->n_composites++;
                                const int size = strlen(node.at1.name) + 1;
                                mdf->size_materials_names += size;

                                /* Book a new index table. */
                                if (mdf_settings_index(
                                        MDF_INDEX_INITIALISE_COMPOSITE,
                                        mdf->n_elements,
                                        error_) != PUMAS_RETURN_SUCCESS)
                                        goto clean_and_exit;
                        } else {
                                /* Analyse the index table. */
                                int elements_in = mdf_settings_index(
                                    MDF_INDEX_FINALISE_COMPOSITE, 0, error_);
                                mdf->n_components += elements_in;
                                if (elements_in > mdf->max_components)
                                        mdf->max_components = elements_in;
                                mdf->size_composite += memory_padded_size(
                                    sizeof(struct composite_material) +
                                        mdf->materials_in *
                                            sizeof(struct composite_component),
                                    pad_size);
                        }
                }

                /* Skip others pure closings. */
                if (node.head == 0) continue;

                /* Process the other head nodes. */
                else if (node.key == MDF_KEY_ELEMENT) {
                        mdf->n_elements++;
                        const int size = strlen(node.at1.name) + 1;
                        mdf->size_elements_names +=
                            memory_padded_size(size, pad_size);
                        if (mdf_settings_name(size, 'E', node.at1.name,
                                error_) != PUMAS_RETURN_SUCCESS)
                                goto clean_and_exit;
                } else if (node.key == MDF_KEY_ATOMIC_COMPONENT) {
                        mdf->n_components++;
                        int index =
                            mdf_settings_name(0, 'E', node.at1.name, error_);
                        if (index < 0) {
                                ERROR_VREGISTER(PUMAS_RETURN_UNKNOWN_ELEMENT,
                                    "unknown atomic element `%s' [@%s:%d]",
                                    node.at1.name, mdf->mdf_path, mdf->line);
                                goto clean_and_exit;
                        }
                        if (mdf_settings_index(MDF_INDEX_WRITE_MATERIAL, index,
                                error_) != PUMAS_RETURN_SUCCESS)
                                goto clean_and_exit;
                } else if (node.key == MDF_KEY_COMPOSITE_COMPONENT) {
                        /* Update the components count. */
                        mdf->materials_in++;

                        /* Update the index table. */
                        int index =
                            mdf_settings_name(0, 'M', node.at1.name, error_);
                        if (index < 0) {
                                ERROR_VREGISTER(PUMAS_RETURN_UNKNOWN_MATERIAL,
                                    "unknown material `%s' [@%s:%d]",
                                    node.at1.name, mdf->mdf_path, mdf->line);
                                goto clean_and_exit;
                        }
                        if (mdf_settings_index(MDF_INDEX_UPDATE_COMPOSITE,
                                index, error_) != PUMAS_RETURN_SUCCESS)
                                goto clean_and_exit;
                }
        }
        mdf->line = 0;

        /* Check the content .*/
        if (mdf->n_elements == 0) {
                /* There are no elements or materials. */
                ERROR_VREGISTER(PUMAS_RETURN_INCOMPLETE_FILE,
                    "no elements in MDF file `%s'", mdf->mdf_path);
                goto clean_and_exit;
        } else if (mdf->n_materials == 0) {
                /* There are no elements or materials. */
                ERROR_VREGISTER(PUMAS_RETURN_INCOMPLETE_FILE,
                    "no materials in MDF file `%s'", mdf->mdf_path);
                goto clean_and_exit;
        }

        /* Parse the kinetic data. */
        if (!mdf->dry_mode) {
                /* Format the full path to the 1st energy loss file. */
                const int size_new = offset_dir + size_name + 1;
                if (size_new > size_path) {
                        /* Get enough memory. */
                        char * new_name = reallocate(full_path, size_new);
                        if (new_name == NULL) {
                                ERROR_REGISTER_MEMORY();
                                goto clean_and_exit;
                        }
                        full_path = new_name;
                }
                strcpy(full_path + offset_dir, filename);

                mdf_parse_kinetic(mdf, full_path, error_);
        } else
                mdf->n_energies = 1;

clean_and_exit:
        /* Free the temporary memory and return. */
        mdf_settings_name(-1, 0x0, NULL, error_);
        mdf_settings_index(MDF_INDEX_FREE, 0, error_);
        deallocate(filename);
        deallocate(full_path);
        return error_->code;
}

/**
 * Manage a temporary mapping of material indices.
 *
 * @param operation The operation to perform.
 * @param value     An optional value for the operation.
 * @param error_    The error data.
 * @return The return value depends on the operation.
 */
int mdf_settings_index(int operation, int value, struct error_context * error_)
{
        const int chunk_size = 1024;
        static int * buffer = NULL;
        static int total_size = 0;
        static int free_size = 0;
        static int n_elements = 0;

        if (operation == MDF_INDEX_FREE) {
                /* Free the temporary memory. */
                deallocate(buffer);
                buffer = NULL;
                total_size = 0;
                free_size = 0;
                return PUMAS_RETURN_SUCCESS;
        } else if (operation == MDF_INDEX_FINALISE_MATERIAL) {
                /* Finalise the material index table. */
                int * header = buffer + (total_size - free_size - value - 1);
                *header = value;
                return PUMAS_RETURN_SUCCESS;
        } else if (operation == MDF_INDEX_UPDATE_COMPOSITE) {
                /* Update the composite index table. */
                int i, *table = buffer;
                for (i = 0; i < value; i++) {
                        table += (*table) + 1;
                }
                const int n = *(table++);
                int * has_element =
                    buffer + (total_size - free_size - n_elements);
                for (i = 0; i < n; i++) {
                        const int j = table[i];
                        if (has_element[j] == 0) has_element[j]++;
                }
                return PUMAS_RETURN_SUCCESS;
        } else if (operation == MDF_INDEX_FINALISE_COMPOSITE) {
                /* Finalise the composite index table. */
                int i, n = 0;
                int * table = buffer + (total_size - free_size - n_elements);
                for (i = 0; i < n_elements; i++) n += *(table++);
                free_size += n_elements;
                return n;
        }

        int size = (operation == MDF_INDEX_INITIALISE_COMPOSITE) ? value : 1;
        if (size > free_size) {
                /* Reserve memory for the next item. */
                const int delta_size =
                    ((size - free_size) / chunk_size + 1) * chunk_size;
                total_size += delta_size;
                int * new_buffer =
                    reallocate(buffer, total_size * sizeof(*buffer));
                if (new_buffer == NULL) {
                        deallocate(buffer);
                        buffer = NULL;
                        total_size = 0;
                        free_size = 0;
                        return ERROR_REGISTER_MEMORY();
                }
                buffer = new_buffer;
                free_size += delta_size;
        }

        if (operation == MDF_INDEX_INITIALISE) {
                return PUMAS_RETURN_SUCCESS;
        } else if (operation == MDF_INDEX_WRITE_MATERIAL) {
                /* Write a value to the material index table. */
                int * tail = buffer + (total_size - free_size);
                *tail = value;
                free_size--;
                return PUMAS_RETURN_SUCCESS;
        } else if (operation == MDF_INDEX_INITIALISE_COMPOSITE) {
                /* Book a new index table for a composite material. */
                n_elements = value;
                int * tail = buffer + (total_size - free_size);
                memset(tail, 0x0, n_elements * sizeof(int));
                free_size -= n_elements;
                return PUMAS_RETURN_SUCCESS;
        }

        return ERROR_VREGISTER(
            PUMAS_RETURN_INDEX_ERROR, "invalid operation `%s'", operation);
}

/**
 * Manage a temporary table of material names.
 *
 * @param size   The operation to perform or the size of the new name.
 * @param prefix A prefix for the category of the name (base, composite, ...).
 * @param name   The name to handle.
 * @param error_ The error data.
 * @return The return value depends on the operation performed.
 */
int mdf_settings_name(
    int size, char prefix, const char * name, struct error_context * error_)
{
        const int chunk_size = 4096;
        static char * buffer = NULL;
        static int total_size = 0;
        static int free_size = 0;

        if (size == 0) {
                /* Return the name index. */
                if (buffer == NULL) return -1;
                const char * ptr = buffer;
                int index = 0;
                while (*ptr != 0x0) {
                        if (*ptr == prefix) {
                                if (strcmp(ptr + 1, name) == 0) return index;
                                index++;
                        }
                        ptr += strlen(ptr + 1) + 2;
                }
                return -1;
        } else if (size < 0) {
                /* Free the temporary memory. */
                deallocate(buffer);
                buffer = NULL;
                total_size = 0;
                free_size = 0;
                return PUMAS_RETURN_SUCCESS;
        }

        size += 1;
        if (size + 1 > free_size) {
                /* Reserve memory for the next item. */
                const int delta_size =
                    ((size + 1 - free_size) / chunk_size + 1) * chunk_size;
                total_size += delta_size;
                char * new_buffer = reallocate(buffer, total_size);
                if (new_buffer == NULL) {
                        deallocate(buffer);
                        buffer = NULL;
                        total_size = 0;
                        free_size = 0;
                        return ERROR_REGISTER_MEMORY();
                }
                buffer = new_buffer;
                free_size += delta_size;
        }
        if (buffer == NULL) return -1;

        /* Dump the new item. */
        char * tail = buffer + (total_size - free_size);
        tail[0] = prefix;
        tail++;
        strcpy(tail, name);
        free_size -= size;
        *(tail + size - 1) = 0x0; /* Tag the end of the name list. */
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Parse the kinetic energy values from a dE/dX file.
 *
 * @param mdf    The MDF handle.
 * @param path   The full path to the dE/dX file.
 * @param error_ The error data.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
 */
enum pumas_return mdf_parse_kinetic(
    struct mdf_buffer * mdf, const char * path, struct error_context * error_)
{
        /* Initialise the settings. */
        mdf->n_energies = 0;
        mdf->dcs_model_offset = 0;
        mdf->line = 0;

        /* Open the dedx file. */
        FILE * fid = fopen(path, "r");
        if (fid == NULL) {
                return ERROR_VREGISTER(
                    PUMAS_RETURN_PATH_ERROR, "could not open file `%s'", path);
        }

        /* Skip the header lines. */
        char * buffer = NULL;
        for (mdf->n_energy_loss_header = 0;; mdf->n_energy_loss_header++) {
                io_read_line(fid, &buffer, path, mdf->line, error_);
                mdf->line++;
                if (error_->code != PUMAS_RETURN_SUCCESS) {
                        fclose(fid);
                        return error_->code;
                }
                const char * c;
                int i;
                for (i = 0, c = buffer; (*c == ' ') && (i < 4); i++) c++;
                if (isdigit(*c)) break;
        }

        /* Scan the table. */
        int nk = 1, offset = 1;
        for (;;) {
                /* Check for a comment. */
                if (strstr(buffer, "Minimum ionization") != NULL)
                        goto next_line;
                if (strstr(buffer, "critical energy") != NULL) goto next_line;

                /* parse the new data line. */
                double k;
                if (sscanf(buffer, "%lf", &k) <= 0) goto next_line;
                k *= 1E-03;
                nk++;
                if (k < DCS_MODEL_MIN_KINETIC) offset++;

        next_line:
                /* Check for a new line. */
                io_read_line(fid, &buffer, path, mdf->line, error_);
                mdf->line++;
                if (error_->code != PUMAS_RETURN_SUCCESS) {
                        if ((error_->code == PUMAS_RETURN_END_OF_FILE) ||
                            feof(fid)) {
                                error_->code = PUMAS_RETURN_SUCCESS;
                                break;
                        }
                }
        }
        io_read_line(NULL, NULL, NULL, 0, error_);

        /*  Update the settings. */
        if (error_->code == PUMAS_RETURN_SUCCESS) mdf->line = 0;
        mdf->n_energies = nk;
        mdf->dcs_model_offset = offset;

        fclose(fid);
        return error_->code;
}

/**
 * Parse the atomic elements from a MDF.
 *
 * @param Physics   Handle for physics tables.
 * @param mdf       The MDF handle.
 * @param error     The error data.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
 */
enum pumas_return mdf_parse_elements(const struct pumas_physics * physics,
    struct mdf_buffer * mdf, struct error_context * error_)
{
        /* Loop on the XML nodes. */
        const int m = physics->n_energies - physics->dcs_model_offset;
        const int n =
            DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1 + DCS_SAMPLING_N;
        const int pad_size = sizeof(*(physics->data));
        const int c_mem = memory_padded_size(
            m * n * (N_DEL_PROCESSES - 1) * sizeof(float), pad_size);
        rewind(mdf->fid);
        mdf->left = 0;
        mdf->line = 1;
        mdf->depth = MDF_DEPTH_EXTERN;
        struct mdf_node node;

        int iel = 0;
        for (;;) {
                /* Get the next node. */
                if (mdf_get_node(mdf, &node, error_) != PUMAS_RETURN_SUCCESS)
                        break;

                if (node.key != MDF_KEY_ELEMENT)
                        continue;
                else if ((node.head == 0) && (node.tail == 1)) {
                        if (iel >= physics->n_elements)
                                break;
                        else
                                continue;
                }

                /* Set the element data. */
                if (iel == 0) {
                        char * tmp = (char *)physics->element +
                            memory_padded_size(physics->n_elements *
                                             sizeof(physics->element[0]),
                                         pad_size);
                        physics->element[0] = (struct atomic_element *)tmp;
                } else {
                        const struct atomic_element * e =
                            physics->element[iel - 1];
                        physics->element[iel] =
                            (struct atomic_element *)((char *)(e) +
                                sizeof(struct atomic_element) + c_mem +
                                memory_padded_size(strlen(e->name) + 1,
                                                          pad_size));
                }
                struct atomic_element * e = physics->element[iel];
                e->dcs_data = (float *)e->data;
                e->name = (char *)(e->data) + c_mem;
                memset(e->dcs_data, 0x0, c_mem);
                strcpy(e->name, node.at1.name);
                if ((sscanf(node.at2.Z, "%lf", &(e->Z)) != 1) || (e->Z <= 0.)) {
                        ERROR_VREGISTER(PUMAS_RETURN_VALUE_ERROR,
                            "invalid atomic number `Z' [%s]", node.at2.Z);
                        break;
                }
                if ((sscanf(node.at3.A, "%lf", &(e->A)) != 1) || (e->A <= 0.)) {
                        ERROR_VREGISTER(PUMAS_RETURN_VALUE_ERROR,
                            "invalid atomic mass `A' [%s]", node.at3.A);
                        break;
                }
                if ((sscanf(node.at4.I, "%lf", &(e->I)) != 1) || (e->I <= 0.)) {
                        ERROR_VREGISTER(PUMAS_RETURN_VALUE_ERROR,
                            "invalid Mean Excitation Energy `I' [%s]",
                            node.at4.I);
                        break;
                }
                e->I *= 1E-09;

                /* Increment. */
                iel++;
                if ((node.tail == 1) && (iel >= physics->n_elements)) break;
        }

        return error_->code;
}

/**
 * Parse the base materials from a MDF.
 *
 * @param Physics   Handle for physics tables.
 * @param mdf       The MDF handle.
 * @param error_    The error data
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
 *
 * The materials dE/dX data are loaded according to the path provided
 * by the `<table>` XML node. If path is not absolute, it specifies a relative
 * path to the directory where the MDF is located.
 */
enum pumas_return mdf_parse_materials(struct pumas_physics * physics,
    struct mdf_buffer * mdf, struct error_context * error_)
{
        /* Format the path directory name. */
        char * filename = NULL;
        int offset_dir, size_name;
        if (!mdf->dry_mode) {
                if ((mdf_format_path(physics->dedx_path, physics->mdf_path,
                        &filename, &offset_dir, &size_name, error_)) !=
                    PUMAS_RETURN_SUCCESS)
                        return error_->code;
        }

        /* Loop on the XML nodes. */
        rewind(mdf->fid);
        mdf->left = 0;
        mdf->line = 1;
        mdf->depth = 0;
        struct mdf_node node;
        int imat = 0;
        const int pad_size = sizeof(*(physics->data));
        void * tmp_ptr = ((char *)physics->composition) +
            memory_padded_size(physics->n_materials *
                                 sizeof(struct material_component *),
                             pad_size);
        struct material_component * data = tmp_ptr;

        for (;;) {
                /* Get the next node. */
                mdf_get_node(mdf, &node, error_);
                if (error_->code != PUMAS_RETURN_SUCCESS) break;

                /* Set the material data. */
                if (node.key == MDF_KEY_MATERIAL) {
                        if (node.head == 0) {
                                /* This is a material closing. */
                                physics->elements_in[imat] = mdf->elements_in;

                                /* Compute the relative electron density. */
                                compute_ZoA(physics, imat);

                                /* Check for dry mode. */
                                if (mdf->dry_mode) goto update_count;

                                /* Read the energy loss data. */
                                FILE * fid = fopen(filename, "r");
                                if (fid == NULL) {
                                        ERROR_VREGISTER(PUMAS_RETURN_PATH_ERROR,
                                            "could not open file `%s'",
                                            filename);
                                        break;
                                }
                                io_parse_dedx_file(
                                    physics, fid, imat, filename, error_);
                                fclose(fid);
                                if (error_->code != PUMAS_RETURN_SUCCESS) break;

                        /* Update the material count. */
                        update_count:
                                imat++;
                                if (imat >= physics->n_materials -
                                        physics->n_composites)
                                        break;
                                continue;
                        }

                        /* We have a new material opening. */
                        if (imat == 0) {
                                physics->material_name[0] =
                                    (char *)(physics->material_name +
                                        physics->n_materials);
                        } else {
                                physics->material_name[imat] =
                                    physics->material_name[imat - 1] +
                                    strlen(physics->material_name[imat - 1]) +
                                    1;
                        }
                        strcpy(physics->material_name[imat], node.at1.name);
                        physics->composition[imat] = data;

                        /* Format the energy loss filename. */
                        if (!mdf->dry_mode) {
                                const int size_new =
                                    offset_dir + strlen(node.at2.file) + 1;
                                if (size_new > size_name) {
                                        /* Get enough memory. */
                                        char * new_name =
                                            reallocate(filename, size_new);
                                        if (new_name == NULL) {
                                                ERROR_REGISTER_MEMORY();
                                                break;
                                        }
                                        filename = new_name;
                                        size_name = size_new;
                                }
                                strcpy(filename + offset_dir, node.at2.file);
                        } else {
                                int n = strlen(node.at2.file) + 1;
                                physics->dedx_filename[imat] = allocate(n);
                                if (physics->dedx_filename[imat] == NULL) {
                                        ERROR_REGISTER_MEMORY();
                                        break;
                                }
                                memcpy(physics->dedx_filename[imat],
                                    node.at2.file, n);
                        }

                        /* Parse the default density. */
                        double rho;
                        if ((sscanf(node.at3.density, "%lf", &rho) != 1)
                            || (rho <= 0.)) {
                                ERROR_VREGISTER(
                                    PUMAS_RETURN_VALUE_ERROR,
                                    "invalid value for density [%s]",
                                   node.at3.density);
                                break;
                        }
                        rho *= 1E+03; /* g/cm^3 -> kg/m^3 */
                        physics->material_density[imat] = rho;
                }

                /* Skip other closings. */
                if (node.head == 0) continue;

                /* Set the composition data. */
                if (node.key == MDF_KEY_ATOMIC_COMPONENT) {
                        int i = mdf->elements_in - 1;
                        int iel = element_index(physics, node.at1.name);
                        if (iel < 0) return PUMAS_RETURN_UNKNOWN_ELEMENT;
                        physics->composition[imat][i].element = iel;

                        double f;
                        if ((sscanf(node.at2.fraction, "%lf", &f) != 1) ||
                            (f <= 0.)) {
                                ERROR_VREGISTER(PUMAS_RETURN_VALUE_ERROR,
                                    "invalid value for mass fraction [%s]",
                                    node.at2.fraction);
                                break;
                        }
                        physics->composition[imat][i].fraction = f;
                        data++;
                }
        }

        if (error_->code != PUMAS_RETURN_SUCCESS) {
                /* Clear the filenames, whenever allocated. */
                int i;
                for (i = 0; i < physics->n_materials - physics->n_composites;
                     i++) {
                        deallocate(physics->dedx_filename[imat]);
                        physics->dedx_filename[imat] = NULL;
                }
        }

        deallocate(filename);
        return error_->code;
}

/**
 * Parse the composite materials from a MDF.
 *
 * @param Physics   Handle for physics tables.
 * @param mdf       The MDF handle.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise an error code.
 *
 * The composite materials properties are given as a linear combination
 * of the base material ones.
 */
enum pumas_return mdf_parse_composites(struct pumas_physics * physics,
    struct mdf_buffer * mdf, struct error_context * error_)
{
        /* Loop on the XML nodes. */
        rewind(mdf->fid);
        mdf->left = 0;
        mdf->line = 1;
        mdf->depth = 0;
        struct mdf_node node;
        int icomp = 0;
        int elements_in = 0;
        const int imat = physics->n_materials - physics->n_composites - 1;
        struct material_component * data_el =
            (struct material_component *)((char *)(physics->composition[imat]) +
                physics->elements_in[imat] * sizeof(struct material_component));
        const int pad_size = sizeof(*(physics->data));
        void * tmp_ptr = ((char *)physics->composite) +
            memory_padded_size(physics->n_composites *
                                 sizeof(struct composite_material *),
                             pad_size);
        struct composite_material * data_ma = tmp_ptr;

        for (;;) {
                /* Get the next node. */
                if (mdf_get_node(mdf, &node, error_) != PUMAS_RETURN_SUCCESS)
                        break;

                /* Set the material data. */
                if (node.key == MDF_KEY_COMPOSITE) {
                        const int i0 = imat + icomp + 1;

                        if (node.head == 0) {
                                /* This is a composite closing. */
                                data_ma->n_components = mdf->materials_in;
                                physics->elements_in[i0] = elements_in;
                                data_ma =
                                    (struct composite_material *)((char *)
                                                                      data_ma +
                                        sizeof(struct composite_material) +
                                        mdf->materials_in *
                                            sizeof(struct composite_component));
                                data_el += elements_in;

                                /* update the composite material count. */
                                icomp++;
                                if (icomp >= physics->n_composites) break;
                                continue;
                        }

                        /* We have a new material opening. */
                        physics->material_name[i0] =
                            physics->material_name[i0 - 1] +
                            strlen(physics->material_name[i0 - 1]) + 1;
                        strcpy(physics->material_name[i0], node.at1.name);
                        physics->composition[i0] = data_el;
                        physics->composite[icomp] = data_ma;
                        elements_in = 0;
                }

                /* Skip other closings. */
                if (node.head == 0) continue;

                /* Set the composition data. */
                if (node.key == MDF_KEY_COMPOSITE_COMPONENT) {
                        int i = mdf->materials_in - 1;
                        int ima = 0;
                        if (material_index(physics, node.at1.name, &ima,
                                error_) != PUMAS_RETURN_SUCCESS)
                                break;
                        data_ma->component[i].material = ima;
                        double f;
                        if ((sscanf(node.at2.fraction, "%lf", &f) != 1) ||
                            (f < 0.)) {
                                ERROR_VREGISTER(PUMAS_RETURN_VALUE_ERROR,
                                    "invalid mass fraction [%s]",
                                    node.at2.fraction);
                                break;
                        }
                        data_ma->component[i].fraction = f;

                        const int i0 = imat + icomp + 1;
                        for (i = 0; i < physics->elements_in[ima]; i++) {
                                int iel = physics->composition[ima][i].element;
                                int j, n = elements_in;
                                int already_listed = 0;
                                for (j = 0; j < n; j++) {
                                        if (iel ==
                                            physics->composition[i0]
                                                                [j].element) {
                                                already_listed = 1;
                                                break;
                                        }
                                }
                                if (already_listed != 0) continue;
                                physics->composition[i0][elements_in++]
                                    .element = iel;
                        }
                }
        }

        if (error_->code == PUMAS_RETURN_END_OF_FILE)
                error_->code = PUMAS_RETURN_SUCCESS;
        return error_->code;
}

/**
 * Get the next XML node in the MDF.
 *
 * @param mdf    The MDF handle.
 * @param node   The node handle.
 * @param error_ The error data.
 * @return On success `PUMAS_RETURN_SUCCES` otherwise the corresponding
 * error code.
 *
 * Search the file for the next valid XML node. On success, at return the *node*
 * attributes are updated with links to the XML buffer. The parsing enforces
 * consistency of the `<pumas>` node and children. Openings and closures
 * above are not checked, neither the tag names.
 */
enum pumas_return mdf_get_node(struct mdf_buffer * mdf, struct mdf_node * node,
    struct error_context * error_)
{
        /* Initialise the node. */
        node->at1.name = node->at2.Z = node->at3.A = node->at4.I = NULL;
        node->key = -1;
        node->head = node->tail = 0;

        for (;;) {
                /* Load data if empty buffer. */
                if (mdf->left < 1) {
                        mdf->left = fread(
                            mdf->data, sizeof(char), mdf->size - 1, mdf->fid);
                        if (mdf->left <= 0)
                                return ERROR_REGISTER_EOF(mdf->mdf_path);
                        mdf->pos = mdf->data;
                        *(mdf->data + mdf->left) = '\0';
                }

                /* Locate a node start tag. */
                char * start = mdf->pos;
                while (*start != '<') {
                        if (*start == '\0')
                                break;
                        else if (*start == '\n') {
                                if (start[1] != '\r') mdf->line++;
                        } else if (*start == '\r') {
                                if (start[1] != '\n') mdf->line++;
                        }
                        start++;
                }
                if (*start == '\0') {
                        mdf->left = 0;
                        continue;
                }
                mdf->left -= start - mdf->pos;
                mdf->pos = start;

                /* Relocate the node to the start of the buffer. */
                memmove(mdf->data, mdf->pos, mdf->left);
                mdf->left += fread(mdf->data + mdf->left, sizeof(char),
                    mdf->size - 1 - mdf->left, mdf->fid);
                if (mdf->left <= 1) return ERROR_REGISTER_EOF(mdf->mdf_path);
                mdf->pos = mdf->data;
                *(mdf->data + mdf->left) = '\0';

                /* Parse the node name. */
                char * key = mdf->pos + 1;
                char * tail = key;
                if (strncmp(key, "!--", 3) == 0) {
                        /* We have a comment block. */
                        tail = key + 3;
                } else {
                        /* We have a regular node. */
                        while ((*tail != ' ') && (*tail != '>')) {
                                if (*tail == '\0')
                                        return ERROR_REGISTER_TOO_LONG(
                                            mdf->mdf_path, mdf->line);
                                tail++;
                        }
                }
                char tailler = *tail;
                *tail = '\0';
                tail++;
                mdf->left -= tail - mdf->pos;
                mdf->pos = tail;

                /* Check if this is a comment. */
                if (strcmp(key, "!--") == 0) {
                        enum pumas_return rc =
                            mdf_skip_pattern(mdf, "-->", error_);
                        if (rc != PUMAS_RETURN_SUCCESS) return rc;
                        continue;
                }

                /* Check if we have a head or tail key. */
                if (*key == '/') {
                        node->tail = 1;
                        key++;
                } else {
                        node->head = 1;
                }

                /* Check the node key. */
                node->key = MDF_KEY_OTHER;
                if (mdf->depth == MDF_DEPTH_EXTERN) {
                        if (strcmp(key, "pumas") == 0)
                                node->key = MDF_KEY_PUMAS;
                        else
                                node->key = MDF_KEY_OTHER;
                } else if (mdf->depth == MDF_DEPTH_ROOT) {
                        if (strcmp(key, "pumas") == 0)
                                node->key = MDF_KEY_PUMAS;
                        else if (strcmp(key, "element") == 0)
                                node->key = MDF_KEY_ELEMENT;
                        else if (strcmp(key, "material") == 0)
                                node->key = MDF_KEY_MATERIAL;
                        else if (strcmp(key, "composite") == 0)
                                node->key = MDF_KEY_COMPOSITE;
                        else
                                return ERROR_REGISTER_UNEXPECTED_TAG(
                                    key, mdf->mdf_path, mdf->line);
                } else if (mdf->depth == MDF_DEPTH_ELEMENT) {
                        if (strcmp(key, "element") == 0)
                                node->key = MDF_KEY_ELEMENT;
                        else
                                return ERROR_REGISTER_UNEXPECTED_TAG(
                                    key, mdf->mdf_path, mdf->line);
                } else if (mdf->depth == MDF_DEPTH_MATERIAL) {
                        if (strcmp(key, "material") == 0)
                                node->key = MDF_KEY_MATERIAL;
                        else if (strcmp(key, "component") == 0)
                                node->key = MDF_KEY_ATOMIC_COMPONENT;
                        else
                                return ERROR_REGISTER_UNEXPECTED_TAG(
                                    key, mdf->mdf_path, mdf->line);
                } else if (mdf->depth == MDF_DEPTH_COMPOSITE) {
                        if (strcmp(key, "composite") == 0)
                                node->key = MDF_KEY_COMPOSITE;
                        else if (strcmp(key, "component") == 0)
                                node->key = MDF_KEY_COMPOSITE_COMPONENT;
                        else
                                return ERROR_REGISTER_UNEXPECTED_TAG(
                                    key, mdf->mdf_path, mdf->line);
                } else
                        return ERROR_REGISTER_UNEXPECTED_TAG(
                            key, mdf->mdf_path, mdf->line);

                /* Check if we have an empty node. */
                if (tailler == '>')
                        goto consistency_check;
                else if ((node->head == 0) && (node->tail == 1)) {
                        /* We have a badly finished pure tailler. */
                        return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                            "unmatched XML closer [@%s:%d]", mdf->mdf_path,
                            mdf->line);
                }

                /* Parse any attributes. */
                for (;;) {
                        char * sep = mdf->pos;
                        while ((*sep != '=') && (*sep != '>')) {
                                if (*sep == '\0')
                                        return ERROR_REGISTER_TOO_LONG(
                                            mdf->mdf_path, mdf->line);
                                sep++;
                        }
                        if (*sep == '>') {
                                /* The end tag was reached. */
                                if ((*(sep - 1) == '/') || (*(sep - 1) == '?'))
                                        node->tail = 1;
                                sep++;
                                mdf->left -= sep - mdf->pos;
                                mdf->pos = sep;
                                goto consistency_check;
                        }

                        /* Parse the attribute and its value. */
                        char * attr = mdf->pos;
                        while (*attr == ' ') attr++;
                        char * end = attr + 1;
                        while ((*end != ' ') && (*end != '=')) end++;
                        *end = '\0';

                        char * value = strchr(sep + 1, '"');
                        if (value == NULL)
                                return ERROR_REGISTER_INVALID_XML_VALUE(
                                    sep + 1, mdf->mdf_path, mdf->line);
                        value++;
                        end = strchr(value, '"');
                        if (end == NULL)
                                return ERROR_REGISTER_INVALID_XML_VALUE(
                                    value, mdf->mdf_path, mdf->line);
                        *end = '\0';

                        /* Map the value to the attribute. */
                        if (node->key == MDF_KEY_PUMAS) {
                                return PUMAS_RETURN_FORMAT_ERROR;
                        } else if (node->key == MDF_KEY_ELEMENT) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else if (strcmp(attr, "Z") == 0)
                                        node->at2.Z = value;
                                else if (strcmp(attr, "A") == 0)
                                        node->at3.A = value;
                                else if (strcmp(attr, "I") == 0)
                                        node->at4.I = value;
                                else
                                        return ERROR_REGISTER_INVALID_XML_ATTRIBUTE(
                                            attr, "<element>", mdf->mdf_path,
                                            mdf->line);
                        } else if (node->key == MDF_KEY_MATERIAL) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else if (strcmp(attr, "file") == 0)
                                        node->at2.file = value;
                                else if (strcmp(attr, "density") == 0)
                                        node->at3.density = value;
                                else
                                        return ERROR_REGISTER_INVALID_XML_ATTRIBUTE(
                                            attr, "<material>", mdf->mdf_path,
                                            mdf->line);
                        } else if (node->key == MDF_KEY_COMPOSITE) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else
                                        return ERROR_REGISTER_INVALID_XML_ATTRIBUTE(
                                            attr, "<composite>", mdf->mdf_path,
                                            mdf->line);
                        } else if (node->key == MDF_KEY_ATOMIC_COMPONENT) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else if (strcmp(attr, "fraction") == 0)
                                        node->at2.fraction = value;
                                else
                                        return ERROR_REGISTER_INVALID_XML_ATTRIBUTE(
                                            attr, "<component>", mdf->mdf_path,
                                            mdf->line);
                        } else if (node->key == MDF_KEY_COMPOSITE_COMPONENT) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else if (strcmp(attr, "fraction") == 0)
                                        node->at2.fraction = value;
                                else
                                        return ERROR_REGISTER_INVALID_XML_ATTRIBUTE(
                                            attr, "<component>", mdf->mdf_path,
                                            mdf->line);
                        }

                        /* Update the buffer cursor. */
                        end++;
                        mdf->left -= end - mdf->pos;
                        mdf->pos = end;
                }
        }

consistency_check:
        /* Update the depth. */
        if (mdf->depth == MDF_DEPTH_EXTERN) {
                if (node->key == MDF_KEY_PUMAS) {
                        if (node->tail == 1)
                                return ERROR_VREGISTER(
                                    PUMAS_RETURN_FORMAT_ERROR,
                                    "invalid <pumas> XML element [@%s:%d]",
                                    mdf->mdf_path, mdf->line);
                        mdf->depth = MDF_DEPTH_ROOT;
                }
        } else if (mdf->depth == MDF_DEPTH_ROOT) {
                if (node->key == MDF_KEY_PUMAS) {
                        if (node->head == 1)
                                return ERROR_VREGISTER(
                                    PUMAS_RETURN_FORMAT_ERROR,
                                    "nested <pumas> XML element [@%s:%d]",
                                    mdf->mdf_path, mdf->line);
                        mdf->depth = MDF_DEPTH_EXTERN;
                } else if (node->key == MDF_KEY_ELEMENT) {
                        if (node->tail == 0) mdf->depth = MDF_DEPTH_ELEMENT;
                } else if (node->key == MDF_KEY_MATERIAL) {
                        if (node->tail == 1)
                                return ERROR_VREGISTER(
                                    PUMAS_RETURN_FORMAT_ERROR,
                                    "invalid <material> XML element [@%s:%d]",
                                    mdf->mdf_path, mdf->line);
                        mdf->depth = MDF_DEPTH_MATERIAL;
                } else if (node->key == MDF_KEY_COMPOSITE) {
                        if (node->tail == 1)
                                return ERROR_VREGISTER(
                                    PUMAS_RETURN_FORMAT_ERROR,
                                    "invalid <composite> XML element [@%s:%d]",
                                    mdf->mdf_path, mdf->line);
                        mdf->depth = MDF_DEPTH_COMPOSITE;
                }
        } else if (mdf->depth == MDF_DEPTH_ELEMENT) {
                if (node->key == MDF_KEY_ELEMENT) {
                        if (node->head == 1)
                                return ERROR_VREGISTER(
                                    PUMAS_RETURN_FORMAT_ERROR,
                                    "nested <element> XML element [@%s:%d]",
                                    mdf->mdf_path, mdf->line);
                        mdf->depth = MDF_DEPTH_ROOT;
                }
        } else if (mdf->depth == MDF_DEPTH_MATERIAL) {
                if (node->key == MDF_KEY_MATERIAL) {
                        if (node->head == 1)
                                return ERROR_VREGISTER(
                                    PUMAS_RETURN_FORMAT_ERROR,
                                    "nested <material> XML element [@%s:%d]",
                                    mdf->mdf_path, mdf->line);
                        mdf->depth = MDF_DEPTH_ROOT;
                }
        } else if (mdf->depth == MDF_DEPTH_COMPOSITE) {
                if (node->key == MDF_KEY_COMPOSITE) {
                        if (node->head == 1)
                                return ERROR_VREGISTER(
                                    PUMAS_RETURN_FORMAT_ERROR,
                                    "nested <composite> XML element [@%s:%d]",
                                    mdf->mdf_path, mdf->line);
                        mdf->depth = MDF_DEPTH_ROOT;
                }
        } else
                return PUMAS_RETURN_FORMAT_ERROR;

        /* Check that the node has all its attributes defined. */
        if (node->head == 0) {
                if ((node->key == MDF_KEY_MATERIAL) && (mdf->elements_in == 0))
                        return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                            "missing <element>(s) for XML <material> [@%s:%d]",
                            mdf->mdf_path, mdf->line);
                else if ((node->key == MDF_KEY_COMPOSITE) &&
                    (mdf->materials_in == 0))
                        return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                            "missing <material>(s) for XML <composite> "
                            "[@%s:%d]",
                            mdf->mdf_path, mdf->line);
                return PUMAS_RETURN_SUCCESS;
        }

        if (node->key == MDF_KEY_COMPOSITE) {
                if (node->at1.name == NULL)
                        return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                            "missing attribute(s) for XML <composite> [@%s:%d]",
                            mdf->mdf_path, mdf->line);
        } else if ((node->key == MDF_KEY_ATOMIC_COMPONENT) ||
            (node->key == MDF_KEY_COMPOSITE_COMPONENT)) {
                if ((node->at1.name == NULL) || (node->at2.fraction == NULL))
                        return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                            "missing attribute(s) for XML <component> [@%s:%d]",
                            mdf->mdf_path, mdf->line);
        } else if (node->key == MDF_KEY_MATERIAL) {
                if ((node->at1.name == NULL) || (node->at2.fraction == NULL) ||
                    (node->at3.density == NULL))
                        return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                            "missing attribute(s) for XML <component> [@%s:%d]",
                            mdf->mdf_path, mdf->line);
        } else if (node->key == MDF_KEY_ELEMENT) {
                if ((node->at1.name == NULL) || (node->at2.Z == NULL) ||
                    (node->at3.A == NULL) || (node->at4.I == NULL))
                        return ERROR_VREGISTER(PUMAS_RETURN_FORMAT_ERROR,
                            "missing attribute(s) for XML <element> [@%s:%d]",
                            mdf->mdf_path, mdf->line);
        }

        /* Update analysis data. */
        if (node->key == MDF_KEY_MATERIAL) {
                mdf->elements_in = 0;
        } else if (node->key == MDF_KEY_COMPOSITE) {
                mdf->materials_in = 0;
        } else if (node->key == MDF_KEY_ATOMIC_COMPONENT) {
                mdf->elements_in++;
        } else if (node->key == MDF_KEY_COMPOSITE_COMPONENT) {
                mdf->materials_in++;
        }

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Skip chars in MDF up to and including a given pattern.
 *
 * @param mdf     The MDF handle.
 * @param pattern The pattern to skip.
 * @param error_  The error data.
 * @return On success `PUMAS_RETURN_SUCCES` otherwise the corresponding
 * error code.
 *
 * Search the file for *pattern* and skip it. This routine is used in
 * order to skip XML comments.
 */
enum pumas_return mdf_skip_pattern(struct mdf_buffer * mdf,
    const char * pattern, struct error_context * error_)
{
        const int pattern_size = strlen(pattern);
        for (;;) {
                if (mdf->left < pattern_size) {
                        memmove(mdf->data, mdf->pos, mdf->left);
                        mdf->left += fread(mdf->data + mdf->left, sizeof(char),
                            mdf->size - 1 - mdf->left, mdf->fid);
                        if (mdf->left < pattern_size) {
                                return ERROR_REGISTER_EOF(mdf->mdf_path);
                        }
                        mdf->pos = mdf->data;
                        *(mdf->data + mdf->left) = '\0';
                }

                char * p = mdf->pos;
                while (*p != '\0') {
                        if (*p == '\n') {
                                if (p[1] != '\r') mdf->line++;
                        } else if (*p == '\r') {
                                if (p[1] != '\n') mdf->line++;
                        }
                        if (strncmp(p, pattern, pattern_size) == 0) break;
                        p++;
                }

                if (*p != '\0') {
                        p += pattern_size;
                        mdf->left -= p - mdf->pos;
                        mdf->pos = p;
                        return PUMAS_RETURN_SUCCESS;
                }
                mdf->pos += mdf->left - pattern_size + 1;
                mdf->left = pattern_size - 1;
        }
}

/**
 * Format the path for a dE/dX data file.
 *
 * @param directory The base directory where the dE/dX file is located.
 * @param filename  The formated path for the filename.
 * @param offset    The offset to the file name part of the path.
 * @param size_name The total memory size available for the path string.
 * @param error_    The error data.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
 *
 * The full path to the filename string is formated by prepending the
 * relative or absolute directory path. If directory is a relative path,
 * the directory of the MDF is used as root.
 */
enum pumas_return mdf_format_path(const char * directory, const char * mdf_path,
    char ** filename, int * offset_dir, int * size_name,
    struct error_context * error_)
{
        /* Format the path directory name. */
        const char sep =
#if (defined _WIN32) || (defined WIN32) || (defined __CYGWIN__)
            '\\';
#else
            '/';
#endif
        const int initial_size = 128;
        if (directory[0] != '@') {
                /* We have an absolute path name. */
                const int dir_size = strlen(directory);
                *offset_dir = dir_size + 1;
                *size_name = (*offset_dir) + initial_size;
                *filename = allocate(*size_name);
                if (*filename == NULL) {
                        return ERROR_REGISTER_MEMORY();
                }
                strcpy(*filename, directory);
                (*filename)[(*offset_dir) - 1] = sep;
        } else {
                /* We have a relative path name. */
                const int dir_size = strlen(++directory);
                int n1 = strlen(mdf_path) - 1;
                while ((n1 >= 0) && (mdf_path[n1] != sep)) n1--;
                *offset_dir = n1 + dir_size + 2;
                *size_name = (*offset_dir) + initial_size;
                *filename = allocate(*size_name);
                if (*filename == NULL) {
                        return ERROR_REGISTER_MEMORY();
                }
                if (n1 > 0) strncpy(*filename, mdf_path, n1 + 1);
                strcpy((*filename) + n1 + 1, directory);
                (*filename)[(*offset_dir) - 1] = sep;
        }

        return PUMAS_RETURN_SUCCESS;
}

/*
 * Low level routines: pre-computation of various properties.
 */
/**
 * Precompute the integrals for a deterministic CEL and TT parameters, for
 * composite materials.
 *
 * @param Physics  Handle for physics tables.
 * @param material The composite material index.
 * @param error_   The error data.
 * @return `PUMAS_ERROR` if any computation failled, `PUMAS_RETURN_SUCCESS`
 * otherwise.
 */
enum pumas_return compute_composite(
    struct pumas_physics * physics, int material, struct error_context * error_)
{
        compute_composite_weights(physics, material);
        compute_composite_tables(physics, material);
        compute_regularise_del(physics, material);
        int ikin;
        for (ikin = 0; ikin < physics->n_energies; ikin++) {
                const enum pumas_return rc =
                    compute_coulomb_parameters(physics, material, ikin, error_);
                if (rc != PUMAS_RETURN_SUCCESS) return rc;
        }
        compute_cel_integrals(physics, material);
        compute_csda_magnetic_transport(physics, material);
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Compute various cumulative integrals for a deterministic CEL.
 *
 * @param Physics  Handle for physics tables.
 * @param material The index of the material to tabulate.
 */
void compute_cel_integrals(struct pumas_physics * physics, int material)
{
        compute_cel_grammage_integral(physics, 0, material);
        compute_cel_grammage_integral(physics, 1, material);
        compute_time_integrals(physics, material);
        compute_kinetic_integral(
            physics, table_get_NI_el(physics, PUMAS_MODE_CSDA, material, 0));
        compute_kinetic_integral(physics,
            table_get_NI_el(physics, PUMAS_MODE_HYBRID, material, 0));
        compute_kinetic_integral(
            physics, table_get_NI_in(physics, material, 0));
}

/**
 * Computation of mixture atomic weights for composite materials.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index of the composite to compute.
 */
void compute_composite_weights(struct pumas_physics * physics, int material)
{
        const int icomp =
            material - physics->n_materials + physics->n_composites;
        int i;
        for (i = 0; i < physics->elements_in[material]; i++) {
                struct material_component * component =
                    physics->composition[material] + i;
                component->fraction = 0.;
        }

        for (i = 0; i < physics->composite[icomp]->n_components; i++) {
                struct composite_component * component =
                    physics->composite[icomp]->component + i;
                const int imat = component->material;
                int j;
                for (j = 0; j < physics->elements_in[imat]; j++) {
                        const struct material_component * cij =
                            physics->composition[imat] + j;
                        int k;
                        for (k = 0; k < physics->elements_in[material]; k++) {
                                struct material_component * c =
                                    physics->composition[material] + k;
                                if (c->element == cij->element) {
                                        c->fraction +=
                                            component->fraction * cij->fraction;
                                        break;
                                }
                        }
                }
        }

        /* Compute the relative electron density. */
        compute_ZoA(physics, material);
}

/**
 * Compute the density of a composite material.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index of the composite material.
 * @param error    The error data.
 * @return `PUMAS_ERROR` if any density is not strictly positive,
 * `PUMAS_RETURN_SUCCESS` otherwise.
 */
enum pumas_return compute_composite_density(
    struct pumas_physics * physics, int material, struct error_context * error_)
{
        const int icomp =
            material - physics->n_materials + physics->n_composites;
        int i;
        double rho_inv = 0., nrm = 0.;
        for (i = 0; i < physics->composite[icomp]->n_components; i++) {
                const struct composite_component * component =
                    physics->composite[icomp]->component + i;
                const double component_density =
                    physics->material_density[component->material];
                if (component_density <= 0.)
                        return ERROR_REGISTER_NEGATIVE_DENSITY(
                            physics->material_name[component->material]);
                rho_inv += component->fraction / component_density;
                nrm += component->fraction;
        }

        if (rho_inv <= 0.)
                return ERROR_REGISTER_NEGATIVE_DENSITY(
                    physics->material_name[material]);
        physics->material_density[material] = nrm / rho_inv;
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Computation of tabulated properties for composite materials.
 *
 * @param Physics  Handle for physics tables.
 * @param material The matrial index of the composite to compute.
 */
void compute_composite_tables(struct pumas_physics * physics, int material)
{
        const int icomp =
            material - physics->n_materials + physics->n_composites;

        /* Initialise to zero. */
        int i, row, k0 = 0;
        for (i = 0; i < material; i++) k0 += physics->elements_in[i];
        row = 0;
        *table_get_T(physics, PUMAS_MODE_CSDA, material, row) = 0.;
        *table_get_T(physics, PUMAS_MODE_HYBRID, material, row) = 0.;
        *table_get_NI_in(physics, material, row) = 0.;
        for (row = 0; row < physics->n_energies; row++) {
                *table_get_dE(physics, PUMAS_MODE_CSDA, material, row) = 0.;
                *table_get_dE(physics, PUMAS_MODE_HYBRID, material, row) = 0.;
                *table_get_CS(physics, material, row) = 0.;
                int k, ip;
                for (k = 0; k < physics->elements_in[material]; k++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++)
                                *table_get_CSf(physics, ip, k0 + k, row) = 0.;
        }
        *table_get_Kt(physics, material) = 0.;
        *table_get_a_max(physics, material) = 0.;
        *table_get_b_max(physics, PUMAS_MODE_CSDA, material) = 0.;
        *table_get_b_max(physics, PUMAS_MODE_HYBRID, material) = 0.;

        /* End point statistics */
        for (i = 0; i < physics->composite[icomp]->n_components; i++) {
                struct composite_component * component =
                    physics->composite[icomp]->component + i;
                const int imat = component->material;
                const double kt = *table_get_Kt(physics, imat);

                /* Total cross section and energy loss. */
                for (row = 0; row < physics->n_energies; row++) {
                        *table_get_dE(physics, 0, material, row) +=
                            *table_get_dE(physics, 0, imat, row) *
                            component->fraction;
                        *table_get_dE(physics, 1, material, row) +=
                            *table_get_dE(physics, 1, imat, row) *
                            component->fraction;
                        const double k = *table_get_K(physics, row);
                        const double cs =
                            (k < kt) ? 0. : *table_get_CS(physics, imat, row) *
                                component->fraction;
                        *table_get_CS(physics, material, row) += cs;
                }

                /* Fractional contribution to the cross section. */
                int j, j0 = 0;
                for (j = 0; j < imat; j++) j0 += physics->elements_in[j];

                for (row = 0; row < physics->n_energies; row++) {
                        const double kinetic = *table_get_K(physics, row);
                        if (kinetic < kt) continue;

                        const double cs_tot = *table_get_CS(physics, imat, row);
                        double csf_last = 0.;
                        for (j = 0; j < physics->elements_in[imat]; j++) {
                                /* Locate the element. */
                                const struct material_component * cij =
                                    physics->composition[imat] + j;
                                int k;
                                for (k = 0; k < physics->elements_in[material];
                                     k++) {
                                        struct material_component * c =
                                            k + physics->composition[material];
                                        if (c->element == cij->element) break;
                                }

                                /* Update the fractional cross-section. */
                                int ip;
                                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                        const double csf = *table_get_CSf(
                                            physics, ip, j0 + j, row);
                                        *table_get_CSf(physics, ip, k0 + k,
                                            row) += cs_tot * (csf - csf_last) *
                                            component->fraction;
                                        csf_last = csf;
                                }
                        }
                }

                /* Maximum tabulated energy loss parameters. */
                *table_get_a_max(physics, material) +=
                    *table_get_a_max(physics, imat) * component->fraction;
                *table_get_b_max(physics, 0, material) +=
                    *table_get_b_max(physics, 0, imat) * component->fraction;
                *table_get_b_max(physics, 1, material) +=
                    *table_get_b_max(physics, 1, imat) * component->fraction;
        }

        /* Normalise the fractional contributions to the cross section. */
        for (row = 0; row < physics->n_energies; row++) {
                int k, ip;
                const double cs = *table_get_CS(physics, material, row);
                if (cs <= 0.) {
                        for (k = 0; k < physics->elements_in[material]; k++)
                                for (ip = 0; ip < N_DEL_PROCESSES; ip++)
                                        *table_get_CSf(
                                            physics, ip, k0 + k, row) = 0.;
                        continue;
                }

                double sum = 0.;
                for (k = 0; k < physics->elements_in[material]; k++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                sum += *table_get_CSf(physics, ip, k0 + k, row);
                                const double f = sum / cs;
                                *table_get_CSf(physics, ip, k0 + k, row) =
                                    (f > 1.) ? 1. : f;
                        }

                /* Protect against rounding errors. */
                k = physics->elements_in[material] - 1;
                *table_get_CSf(physics, N_DEL_PROCESSES - 1, k0 + k, row) = 1.;
        }

        /* Weighted integrands */
        for (row = 1; row < physics->n_energies; row++) {
                *table_get_NI_in(physics, material, row) =
                    *table_get_CS(physics, material, row) /
                    *table_get_dE(physics, PUMAS_MODE_HYBRID, material, row);
        }
}

/**
 * Compute the cumulative integral of a table column.
 *
 * @param Physics   Handle for physics tables.
 * @param table     The table column to process.
 *
 * Compute a cumulative integral of the column by integrating over the kinetic
 * energy with a linear interpolation and a trapezoidal rule.
 */
void compute_kinetic_integral(struct pumas_physics * physics, double * table)
{
        double value = 0., dv;
        int i;
        for (i = 1; i < physics->n_energies; i++) {
                const double x0 = *table_get_K(physics, i - 1);
                const double x1 = *table_get_K(physics, i);
                dv = 0.5 * (x1 - x0) * (table[i - 1] + table[i]);
                table[i - 1] = value;
                value += dv;
        }
        table[physics->n_energies - 1] = value;
}

/**
 * Compute cumulative integrals for the proper time.
 *
 * @param Physics   Handle for physics tables.
 * @param table     The table column to process.
 *
 * Compute the proper time cumulative integrals over path length with a
 * trapezoidal rule.
 */
void compute_time_integrals(struct pumas_physics * physics, int material)
{
        static double I0 = 0.;
        if (I0 == 0.) {
                /* Compute the integral of 1/momemtum for the lowest energy bin
                 * using trapezes. */
                const int n = 101;
                int i;
                const double dK = (*table_get_K(physics, 1)) / (n - 1);
                double Ki = dK;
                I0 = 0.5 / sqrt(Ki * (Ki + 2. * physics->mass));
                for (i = 2; i < n - 1; i++) {
                        Ki += dK;
                        const double pi = sqrt(Ki * (Ki + 2. * physics->mass));
                        I0 += 1. / pi;
                }
                Ki += dK;
                I0 += 0.5 / sqrt(Ki * (Ki + 2. * physics->mass));
                I0 /= n - 1;
        }

        /* Compute the cumulative path integrals . */
        double * const K = table_get_K(physics, 0);
        double * const T0 =
            table_get_T(physics, PUMAS_MODE_CSDA, material, 0);
        double * const T1 =
            table_get_T(physics, PUMAS_MODE_HYBRID, material, 0);
        double * const X0 =
            table_get_X(physics, PUMAS_MODE_CSDA, material, 0);
        double * const X1 =
            table_get_X(physics, PUMAS_MODE_HYBRID, material, 0);

        T0[0] = T1[0] = 0.;
        T0[1] = I0 * X0[1] * physics->mass;
        T1[1] = I0 * X1[1] * physics->mass;
        int i;
        for (i = 2; i < physics->n_energies; i++) {
                const double p0 =
                    sqrt(K[i - 1] * (K[i - 1] + 2. * physics->mass));
                const double p1 = sqrt(K[i] * (K[i] + 2. * physics->mass));
                const double psi = 1. / p0 + 1. / p1;
                const double dy0 = 0.5 * (X0[i] - X0[i - 1]) * psi;
                const double dy1 = 0.5 * (X1[i] - X1[i - 1]) * psi;
                T0[i] = T0[i - 1] + dy0 * physics->mass;
                T1[i] = T1[i - 1] + dy1 * physics->mass;
        }
}

/**
 * Compute the CSDA grammage range from the energy loss.
 *
 * @param Physics   Handle for physics tables.
 *  @param scheme   The index of the simulation scheme.
 *  @param material The index of the material to process.
 *
 * Compute the cumulative CSDA grammage integral from the energy loss
 * using a trapezoidal rule.
 */
void compute_cel_grammage_integral(
    struct pumas_physics * physics, int scheme, int material)
{
        const double * const kinetic = table_get_K(physics, 0);
        const double * const dEdX = table_get_dE(physics, scheme, material, 0);
        double * const table = table_get_X(physics, scheme, material, 0);

        /* Compute the cumulative integral. */
        int i;
        table[0] = 0.;
        double y0 = 1. / dEdX[1];
        for (i = 1; i < physics->n_energies; i++) {
                const double y1 = 1. / dEdX[i];
                table[i] = table[i - 1] +
                    0.5 * (kinetic[i] - kinetic[i - 1]) * (y0 + y1);
                y0 = y1;
        }
}

/**
 * Compute the cumulative integrals of the magnetic transport.
 *
 * @param Physics   Handle for physics tables.
 * @param imed      The index of the material to tabulate.
 *
 * Compute the cumulative integrals for the momenta of the magnetic deflection
 * using a trapezoidal rule.
 */
void compute_csda_magnetic_transport(
    struct pumas_physics * physics, int material)
{
        double x[N_LARMOR_ORDERS];
        double dx[N_LARMOR_ORDERS];
        memset(x, 0x0, sizeof(x));

        /* The magnetic phase shift is proportional to the proper time integral.
         * We refer to this table. */
        const double factor = LARMOR_FACTOR / physics->mass;
        double * const T = table_get_T(physics, PUMAS_MODE_CSDA, material, 0);

        /* Compute the deflection starting from max energy down to 0 */
        double * const X0 =
            table_get_X(physics, PUMAS_MODE_CSDA, material, 0);
        const int imax = physics->n_energies - 1;
        int i, j;
        for (i = physics->n_energies - 2; i >= 1; i--) {
                double dX0 = 0.5 * (X0[i + 1] - X0[i]);
                double p1 = (T[imax] - T[i]) * factor;
                double p2 = (T[imax] - T[i + 1]) * factor;

                double f1 = 1., f2 = 1.;
                for (j = 0; j < N_LARMOR_ORDERS; j++) {
                        *table_get_Li(physics, material, j, i + 1) = x[j];
                        dx[j] = dX0 * (f1 + f2);
                        x[j] += dx[j];
                        f1 *= p1;
                        f2 *= p2;
                }
        }

        /* Extrapolate the end points */
        for (j = 0; j < N_LARMOR_ORDERS; j++) {
                *table_get_Li(physics, material, j, 1) = x[j];
                double hx = (X0[1] - X0[0]) / (X0[2] - X0[1]);
                *table_get_Li(physics, material, j, 0) = x[j] + hx * dx[j];
        }
}

/**
 * Compute and tabulate the multiple scattering parameters.
 *
 * @param Physics  Handle for physics tables.
 * @param material The target material.
 * @param row      The kinetic energy index to compute for.
 * @param error    The error data.
 *
 * At output the *row* of `physics::table_Mu0` and `physics::table_Ms1` are
 * updated. This routine manages a temporary workspace memory buffer. Calling
 * it with a negative *material* index causes the workspace memory to be freed.
 */
enum pumas_return compute_coulomb_parameters(struct pumas_physics * physics,
    int material, int row, struct error_context * error_)
{
        /* Handle the memory for the temporary workspace. */
        static struct coulomb_workspace * workspace = NULL;
        if (material < 0) {
                deallocate(workspace);
                workspace = NULL;
                return PUMAS_RETURN_SUCCESS;
        } else if (workspace == NULL) {
                const int work_size = sizeof(struct coulomb_workspace) +
                    physics->max_components * sizeof(struct coulomb_data);
                workspace = allocate(work_size);
                if (workspace == NULL) return ERROR_REGISTER_MEMORY();
        }

        /* Check the kinetic energy. */
        const double kinetic = *table_get_K(physics, row);
        if (kinetic <= 0.) {
                *table_get_Mu0(physics, material, row) = 0.;
                *table_get_Ms1(physics, material, row) = 0.;
                return PUMAS_RETURN_SUCCESS;
        }

        /* Compute the mean free paths, the screening factors and various
         * averages used for the hard scattering.
         */
        double invlb_m = 0., invlb1_m = 0.;
        double s_m_l = 0., s_m_h = 0.;
        int i;
        struct coulomb_data * data;
        for (i = 0, data = workspace->data; i < physics->elements_in[material];
             i++, data++) {
                double G[2];
                const struct material_component * const component =
                    &physics->composition[material][i];
                const struct atomic_element * const element =
                    physics->element[component->element];
                double kinetic0;
                coulomb_frame_parameters(
                    physics, kinetic, element->A, &kinetic0, data->fCM);
                data->fspin = coulomb_spin_factor(physics, kinetic0);
                coulomb_screening_parameters(physics, NULL, kinetic0,
                    component->element, data->screening);
                coulomb_pole_decomposition(data->screening, data->a, data->b);
                coulomb_transport_coefficients(
                    1., data->fspin, data->screening, data->a, data->b, G);
                const double invlb = component->fraction /
                    coulomb_wentzel_path(physics, kinetic0, element->Z,
                                         element->A, data->screening[0]);

                invlb_m += invlb * G[0];
                s_m_h += data->screening[0] * invlb;
                s_m_l += invlb / data->screening[0];
                data->invlambda = invlb;
                const double d = 1. / (data->fCM[0] * (1. + data->fCM[1]));
                invlb1_m += invlb * G[1] * d * d;
        }

        /* Set the hard scattering mean free path. */
        const double lb_m = 1. / invlb_m;
        double lb_h = EHS_OVER_MSC / invlb1_m;
        if (lb_h > EHS_PATH_MAX) lb_h = EHS_PATH_MAX;

        /* Compute the hard scattering cutoff angle, in the CM. */
        const double max_mu0 = 0.5 * (1. - cos(MAX_SOFT_ANGLE * M_PI / 180.));
        double mu0 = 0.;
        if (lb_m < lb_h) {
                /* Initialise the root finder with an asymptotic starting
                 * value. */
                double s_m;
                if (lb_h > 2. * lb_m) /* Asymptotic value when lb_h >> lb_m */
                        s_m = s_m_h * lb_m;
                else
                        /* Asymptotic value when lb_h ~= lb_m */
                        s_m = 1. / (s_m_l * lb_m);
                mu0 = s_m * (lb_h - lb_m) / (s_m * lb_h + lb_m);

                /* Configure for the root solver. */
                workspace->cs_h = 1. / lb_h;
                workspace->material = material;
                workspace->ihard = -1;

                /* Solve for the cut-off angle. We try an initial bracketing in
                 * [0.25*mu0; 4.*mu0], with mu0 the asymptotic estimate. */
                double mubest;
                double mu_max = 4. * mu0;
                if (mu_max > 1.) mu_max = 1.;
                double mu_min = 0.25 * mu0;
                double fmax =
                    compute_cutoff_objective(physics, mu_max, workspace);
                double fmin;
                if (fmax > 0.) {
                        /* This shouldn't occur, but let's be safe and handle
                         * this case. */
                        mu_min = mu_max;
                        fmin = fmax;
                        mu_max = 1.;
                        fmax = -workspace->cs_h;
                } else {
                        fmin = compute_cutoff_objective(
                            physics, mu_min, workspace);
                        if (fmin < 0.) {
                                /* This might occur at high energies when the
                                 * nuclear screening becomes significant. */
                                mu_max = mu_min;
                                fmax = fmin;
                                mu_min = 0.;
                                fmin = compute_cutoff_objective(
                                    physics, mu_min, workspace);
                        }
                }
                if (mu_min < max_mu0) {
                        if (mu_max > max_mu0) mu_max = max_mu0;
                        if (math_find_root(compute_cutoff_objective, physics,
                                mu_min, mu_max, &fmin, &fmax, 1E-06 * mu0,
                                1E-06, 100, workspace, &mubest) == 0)
                                mu0 = mubest;
                }
                if (mu0 > max_mu0) mu0 = max_mu0;
                lb_h = compute_cutoff_objective(physics, mu0, workspace) +
                    workspace->cs_h;
                if (lb_h <= 1. / EHS_PATH_MAX)
                        lb_h = EHS_PATH_MAX;
                else
                        lb_h = 1. / lb_h;
        } else
                lb_h = lb_m;
        *table_get_Mu0(physics, material, row) = mu0;
        *table_get_Lb(physics, material, row) =
            lb_h * kinetic * (kinetic + 2. * physics->mass);
        *table_get_NI_el(physics, PUMAS_MODE_CSDA, material, row) = 1. /
            (*table_get_dE(physics, PUMAS_MODE_CSDA, material, row) * lb_h);
        *table_get_NI_el(physics, PUMAS_MODE_HYBRID, material, row) = 1. /
            (*table_get_dE(physics, PUMAS_MODE_HYBRID, material, row) * lb_h);

        /* Compute the 1st moment of the soft scattering. */
        const int n0 = physics->n_materials - physics->n_composites;
        if (material < n0) {
                /* We have a base material. */
                static double * ms1_table = NULL;
                if (material == 0) {
                        /* Precompute the per element soft scattering terms. */
                        enum pumas_return rc;
                        if ((rc = compute_coulomb_soft(physics, row, &ms1_table,
                                 error_)) != PUMAS_RETURN_SUCCESS)
                                return rc;
                }

                double invlb1 = 0.;
                struct material_component * component =
                    physics->composition[material];
                for (i = 0, data = workspace->data;
                     i < physics->elements_in[material];
                     i++, data++, component++) {
                        /* Screened nucleus contribution to the transport. */
                        double G[2];
                        coulomb_transport_coefficients(mu0, data->fspin,
                            data->screening, data->a, data->b, G);
                        double d = 1. / (data->fCM[0] * (1. + data->fCM[1]));
                        d *= d;
                        invlb1 += data->invlambda * d * G[1];

                        /* Other precomputed soft scattering terms. */
                        const int iel = component->element;
                        invlb1 += *table_get_ms1(physics, iel, row, ms1_table) *
                            component->fraction;
                }
                *table_get_Ms1(physics, material, row) = invlb1;
        } else {
                /* We have a composite material. Let's loop on the base
                 * material components.
                 */
                const struct composite_material * composite =
                    physics->composite[material - n0];
                double invlb1 = 0.;
                int icomp;
                for (icomp = 0; icomp < composite->n_components; icomp++) {
                        const struct composite_component * c =
                            composite->component + icomp;
                        const int imat = c->material;
                        struct material_component * component =
                            physics->composition[imat];

                        /* Compute the variation of the Coulomb scattering. */
                        const double mu0b = *table_get_Mu0(physics, imat, row);
                        double delta_invlb1 = 0.;
                        for (i = 0, data = workspace->data;
                             i < physics->elements_in[imat];
                             i++, data++, component++) {
                                /* Compute the scattering parameters for the
                                 * base material.
                                 */
                                double G[2];
                                const struct atomic_element * const element =
                                    physics->element[component->element];
                                double kinetic0;
                                coulomb_frame_parameters(physics, kinetic,
                                    element->A, &kinetic0, data->fCM);
                                data->fspin =
                                    coulomb_spin_factor(physics, kinetic0);
                                coulomb_screening_parameters(physics, NULL,
                                    kinetic0, component->element,
                                    data->screening);
                                coulomb_pole_decomposition(
                                    data->screening, data->a, data->b);
                                coulomb_transport_coefficients(1., data->fspin,
                                    data->screening, data->a, data->b, G);
                                const double invlb = component->fraction /
                                    coulomb_wentzel_path(physics, kinetic0,
                                                         element->Z, element->A,
                                                         data->screening[0]);

                                /* Variation of the screened nucleus
                                 * contribution to the transport.
                                 */
                                double deltaG;
                                coulomb_transport_coefficients(mu0, data->fspin,
                                    data->screening, data->a, data->b, G);
                                deltaG = G[1];
                                coulomb_transport_coefficients(mu0b,
                                    data->fspin, data->screening, data->a,
                                    data->b, G);
                                deltaG -= G[1];
                                double d =
                                    1. / (data->fCM[0] * (1. + data->fCM[1]));
                                d *= d;
                                delta_invlb1 += invlb * d * deltaG;
                        }

                        /* Base material soft scattering. */
                        invlb1 += (*table_get_Ms1(physics, imat, row) +
                                      delta_invlb1) *
                            c->fraction;
                }
                *table_get_Ms1(physics, material, row) = invlb1;
        }

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Tabulate the multiple scattering per element.
 *
 * @param Physics   Handle for physics tables.
 * @param row       The row index for the kinetic value.
 * @param data      The tabulated data.
 * @param error_    The error data.
 * @return On success `PUMAS_RETRUN_SUCCESS` otherwise an error code.
 *
 * Because this step is time consuming, the multiple scattering 1st path
 * lengths are precomputed per element at initialisation using a temporary
 * buffer. These temporary tables are used for computing per material 1st path
 * length.
 *
 * **Note** This routine handles a static dynamically allocated table. If the
 * *row* index is negative the table is freed.
 */
enum pumas_return compute_coulomb_soft(struct pumas_physics * physics, int row,
    double ** data, struct error_context * error_)
{
        static double * ms1_table = NULL;
        if (data != NULL) *data = NULL;

        if (row < 0) {
                deallocate(ms1_table);
                ms1_table = NULL;
                return PUMAS_RETURN_SUCCESS;
        }

        /* Allocate the temporary table. */
        if (ms1_table == NULL) {
                ms1_table = allocate(
                    physics->n_elements * physics->n_energies * sizeof(double));
                if (ms1_table == NULL) return ERROR_REGISTER_MEMORY();
        }

        /* Loop over atomic elements. */
        const double kinetic = *table_get_K(physics, row);
        int iel;
        for (iel = 0; iel < physics->n_elements; iel++) {
                struct atomic_element * element = physics->element[iel];
                double invlb1 = 0.;

                /* Electron shells contribution to the transverse transport. */
                invlb1 =
                    transverse_transport_ionisation(physics, element, kinetic);

                /* Photonuclear contribution to the transverse transport. */
                invlb1 += transverse_transport_photonuclear(
                    physics, element, kinetic);

                *table_get_ms1(physics, iel, row, ms1_table) = invlb1;
        }

        *data = ms1_table;
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Objective function for solving the soft scattering cut-off angle.
 *
 * @param Physics    Handle for physics tables.
 * @param mu         The proposed cutoff angular parameter.
 * @param parameters A pointer to the temporary workspace.
 * @return The difference between the current and expected restricted cross
 * section.
 *
 * This is a wrapper for the root solver. It provides the objective function
 * to solve for.
 */
double compute_cutoff_objective(
    const struct pumas_physics * physics, double mu, void * parameters)
{
        /* Unpack the workspace. */
        struct coulomb_workspace * workspace = parameters;

        /* Compute the restricted cross section. */
        double cs_tot = 0.;
        int i, n = physics->elements_in[workspace->material];
        struct coulomb_data * data;
        for (i = 0, data = workspace->data; i < n; i++, data++) {
                cs_tot +=
                    data->invlambda * coulomb_restricted_cs(mu, data->fspin,
                                          data->screening, data->a, data->b);
        }

        /* Return the difference with the expectation. */
        return cs_tot - workspace->cs_h;
}

/**
 * Tabulate the deterministic CEL and DEL cross-sections per element.
 *
 * @param Physics   Handle for physics tables.
 * @param row       The row index for the kinetic value.
 * @return          The tabulated data.
 *
 * Because this step is time consuming, the CEL integration is precomputed per
 * element at initialisation using a temporary buffer. These temporary tables
 * are used for computing per material CEL corrections and DEL cross-sections.
 *
 * **Note** This routine handles a static dynamically allocated table. If the
 * *row* index is negative the table is freed.
 */
double * compute_cel_and_del(struct pumas_physics * physics, int row)
{
        static double * cel_table = NULL;

        if (row < 0) {
                deallocate(cel_table);
                return (cel_table = NULL);
        }

        /* Allocate the temporary table. */
        if (cel_table == NULL) {
                cel_table = allocate(N_DEL_PROCESSES * physics->n_elements *
                    physics->n_energies * sizeof(double));
                if (cel_table == NULL) return NULL;
        }

        /* Loop over atomic elements. */
        const double kinetic = *table_get_K(physics, row);
        int iel;
        for (iel = 0; iel < physics->n_elements; iel++) {
                const struct atomic_element * element = physics->element[iel];

                /* Loop over processes. */
                int ip;
                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                        *table_get_CSn(physics, ip, iel, row) =
                            compute_dcs_integral(physics, 0, element, kinetic,
                                dcs_get(ip), physics->cutoff, 180);
                        *table_get_cel(physics, ip, iel, row, cel_table) =
                            compute_dcs_integral(physics, 1, element, kinetic,
                                dcs_get(ip), physics->cutoff, 180);
                }
        }

        return cel_table;
}

/**
 * Regularise the cross-section table and related quantities with a
 * *do nothing* process.
 *
 * @param material The index of the propagation material.
 */
void compute_regularise_del(struct pumas_physics * physics, int material)
{
        /* Find the kinetic threshold above which the total cross-section is
         * not null. */
        int it;
        double cs0;
        for (it = 1; it < physics->n_energies; it++)
                if ((cs0 = *table_get_CS(physics, material, it)) != 0) break;
        *table_get_Kt(physics, material) = *table_get_K(physics, it);

        /*
         * Regularise the cross section and the total number of interaction
         * lengths.
         */
        int row;
        for (row = 1; row < it; row++) {
                *table_get_CS(physics, material, row) = cs0;
                const double dEdX =
                    *table_get_dE(physics, PUMAS_MODE_HYBRID, material, row);
                *table_get_NI_in(physics, material, row) = cs0 / dEdX;
        }

        if (material > 0)
                return; /* The computations below need to be done only
once. */

        /* Find the kinetic threshold per process and atomic element. */
        int iel, ip;
        for (iel = 0; iel < physics->n_elements; iel++) {
                const struct atomic_element * element = physics->element[iel];
                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                        for (row = 0; row < it; row++)
                                *table_get_Xt(physics, ip, iel, row) = 1.;
                        dcs_function_t * dcs_func = dcs_get(ip);
                        for (row = it; row < physics->n_energies; row++) {
                                const double k = *table_get_K(physics, row);
                                double x = physics->cutoff;
                                while ((x < 1.) && (dcs_func(physics, element,
                                                        k, k * x) <= 0.))
                                        x *= 2;
                                if (x >= 1.)
                                        x = 1.;
                                else if (x > physics->cutoff) {
                                        const double eps = 1E-02 *
                                            physics->cutoff;
                                        double x0 = 0.5 * x;
                                        double dcs = 0.;
                                        for (;;) {
                                                if (dcs == 0.)
                                                        x0 += 0.5 * (x - x0);
                                                else {
                                                        const double dx =
                                                            x - x0;
                                                        x = x0;
                                                        x0 -= 0.5 * dx;
                                                }
                                                if ((x - x0) <= eps) break;
                                                dcs = dcs_func(physics, element,
                                                    k, k * x0);
                                        }
                                }
                                *table_get_Xt(physics, ip, iel, row) = x;
                        }
                }
        }
}

/**
 * Compute integrals of DCSs.
 *
 * @param Physics Handle for physics tables.
 * @param forward The MC flow direction.
 * @param mode    Flag to select the integration mode.
 * @param element The target atomic element.
 * @param kinetic The initial or final kinetic energy.
 * @param dcs     Handle to the dcs function.
 * @param xlow    The lower bound of the fractional energy transfer.
 * @param nint    The requested number of point for the integral.
 * @return The integrated dcs, in m^2/kg or the energy loss in GeV m^2/kg.
 *
 * The parameter *mode* controls the integration mode as following. If *mode*
 * is `0` the restricted cross section is computed, for energy losses larger
 * than *xlow*. Else, the restricted energy loss is computed, for losses
 * larger than *xlow*
 */
double compute_dcs_integral(struct pumas_physics * physics, int mode,
    const struct atomic_element * element, double kinetic, dcs_function_t * dcs,
    double xlow, int nint)
{

        /* Let us use the analytical form for ionisation when radiative
         * corrections can be neglected.
         */
        const double m1 = physics->mass - ELECTRON_MASS;
        if ((dcs == &dcs_ionisation) &&
            (kinetic <= 0.5 * m1 * m1 / ELECTRON_MASS))
                return dcs_ionisation_integrate(
                    physics, mode, element, kinetic, xlow);

        /* We integrate over the recoil energy using a logarithmic sampling. */
        double dcsint = 0.;
        double x0 = log(kinetic * xlow), x1 = log(kinetic);
        math_gauss_quad(nint, &x0, &x1); /* Initialisation. */

        double xi, wi;
        while (math_gauss_quad(0, &xi, &wi) == 0) { /* Iterations. */
                const double qi = exp(xi);
                double y = dcs(physics, element, kinetic, qi) * qi;
                if (mode != 0) y *= qi;
                dcsint += y * wi;
        }
        dcsint /= kinetic + physics->mass;

        return dcsint;
}

/**
 * Compute or update the relative electron density of a material.
 *
 * @param Physics  Handle for physics tables.
 * @param material The material index.
 */
void compute_ZoA(struct pumas_physics * physics, int material)
{
        int i;
        double ZoA = 0.;
        for (i = 0; i < physics->elements_in[material]; i++) {

                const struct material_component * const c =
                    physics->composition[material] + i;
                const struct atomic_element * const e =
                    physics->element[c->element];
                ZoA += e->Z / e->A * c->fraction;
        }
        physics->material_ZoA[material] = ZoA;
}

/**
 * Compute the coefficients of the polynomial approximation to an
 * inelastic DCS.
 *
 * @param Physics  Handle for physics tables.
 * @param dcs_func The DCS function.
 * @param element  The target atomic element.
 * @param error_   The error data.
 *
 * The DCS is also tabulated at specific values for the ziggurat algorithm.
 */
enum pumas_return compute_dcs_model(struct pumas_physics * physics,
    dcs_function_t * dcs_func, struct atomic_element * element,
    struct error_context * error_)
{
        /* Manage temporary storage. */
        static int m = 0;
        static double * tmp;
        const int n = DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1;

        if (dcs_func == NULL) {
                deallocate(tmp);
                tmp = NULL;
                m = 0;
                dcs_model_fit(0, 0, NULL, NULL, NULL, NULL);
                return PUMAS_RETURN_SUCCESS;
        } else if (m <= 0) {
                m = (int)(100. * log10(DCS_MODEL_MAX_FRACTION /
                    physics->cutoff)) + 1 + DCS_SAMPLING_N;
                if (m < 0) return PUMAS_RETURN_INTERNAL_ERROR;
                tmp = allocate((3 * m + n + DCS_SAMPLING_N) * sizeof(double));
                if (tmp == NULL) return ERROR_REGISTER_MEMORY();
        }
        double * x = tmp;
        double * y = x + m;
        double * w = y + m;
        double * c_ = w + m;

        /*  Loop over the kinetic energy values. */
        const int index = dcs_get_index(dcs_func);
        const int nkeff = physics->n_energies - physics->dcs_model_offset;
        float * const coeff = table_get_dcs_coeff(physics, element, index, 0);
        const double x0 = log(physics->cutoff);
        const double dx = log(DCS_MODEL_MAX_FRACTION / physics->cutoff) /
            (m - 1 - DCS_SAMPLING_N);
        int i;
        float * c;
        for (i = 0, c = coeff; i < nkeff; i++, c += n + DCS_SAMPLING_N) {
                /* Prepare the fit values using a log sampling. */
                const double K =
                    *table_get_K(physics, i + physics->dcs_model_offset);
                int j, first = 1, j0 = -1, j1 = -1;
                for (j = 0; j < m - DCS_SAMPLING_N; j++) {
                        const double nu = exp(x0 + j * dx);
                        x[j] = nu;
                        const double dcs =
                            dcs_func(physics, element, K, nu * K) * K /
                            (K + physics->mass);
                        if (dcs > 0.) {
                                y[j] = dcs;
                                if (first) {
                                        first = 0;
                                        j0 = j;
                                } else
                                        j1 = j;
                                w[j] = 1.;
                        } else {
                                y[j] = w[j] = 0.;
                        }
                }
                if ((j0 < 0) || (j1 < 0))
                        return ERROR_REGISTER(PUMAS_RETURN_INTERNAL_ERROR,
                            "failed to fit the DCS");

                /* Add the tabulated values in linear scale. */
                const double dnu = (1. - x[j0]) / DCS_SAMPLING_N;
                double nu;
                for (j = m - DCS_SAMPLING_N, nu = x[j0]; j - m < 0;
                     j++, nu += dnu)
                /* Patch a gcc warning here. With -O2 j < m != j-m < 0 ... */
                {
                        x[j] = nu;
                        const double tmp =
                            dcs_func(physics, element, K, nu * K) * K /
                            (K + physics->mass);
                        c[j + n + DCS_SAMPLING_N - m] = (float)tmp;
                        if (nu <= x[j1]) {
                                y[j] = tmp;
                                w[j] = 1.;
                        } else {
                                y[j] = w[j] = 0.;
                        }
                }

                w[j0] = w[j1] = 1E+06; /*  Constrain the end points. */
                dcs_model_fit(m - DCS_SAMPLING_N, n, x, y, w, c_);
                for (j = 0; j < n; j++) c[j] = (float)c_[j];
        }

        return PUMAS_RETURN_SUCCESS;
}

/*
 * Low level routines: Models for the differential cross-sections and helper
 * functions.
 *
 * Note that the macroscopic DCS per fractional energy loss and per mass unit
 * are actualy used, i.e. dSigma/dnu. The unit is an inverse column density,
 * expressed in m^2/kg. The corresponding interaction length is
 * `Lint = 1/(Sigma*rho)` with rho the medium density.
 */
/**
 * Get a differential cross-section given a process index.
 *
 * @param process The process index.
 * @return The differential cross section function.
 *
 * `index = 0` corresponds to Bremsstrahlung, `index = 1` to Pair Production,
 * `index = 2` to Photonuclear interactions and `index = 3` to Ionisation.
 */
dcs_function_t * dcs_get(int process)
{
        dcs_function_t * const dcs_func[] = { dcs_bremsstrahlung,
                dcs_pair_production, dcs_photonuclear, dcs_ionisation };
        return dcs_func[process];
}

/**
 * Get the table index for a given DCS function.
 *
 * @param dcs_function The DCS function.
 * @return The corresponding table index.
 *
 * This is the inverse mapping than dcs_get.
 */
int dcs_get_index(dcs_function_t * dcs_func)
{
        if (dcs_func == dcs_bremsstrahlung)
                return 0;
        else if (dcs_func == dcs_pair_production)
                return 1;
        else if (dcs_func == dcs_photonuclear)
                return 2;
        return 3;
}

/**
 * Wrapper for the Bremsstrahlung differential cross section.
 *
 * @param Physics Handle for physics tables.
 * @param element The target atomic element.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The differential cross section in m^2/kg.
 */
double dcs_bremsstrahlung(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q)
{
        return physics->dcs_bremsstrahlung(element->Z, element->A,
            physics->mass, K, q) * 1E+03 * AVOGADRO_NUMBER *
            (physics->mass + K) / element->A;
}

/**
 * Wrapper for the pair production differential cross section.
 *
 * @param Physics Handle for physics tables.
 * @param element The target atomic element.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The differential cross section in m^2/kg.
 */
double dcs_pair_production(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q)
{
        return physics->dcs_pair_production(element->Z, element->A,
            physics->mass, K, q) * 1E+03 * AVOGADRO_NUMBER *
            (physics->mass + K) / element->A;
}

/**
 * Wrapper for the photonuclear differential cross section.
 *
 * @param Physics Handle for physics tables.
 * @param element The target atomic element.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The differential cross section in m^2/kg.
 */
double dcs_photonuclear(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q)
{
        if (dcs_photonuclear_check(K, q)) { /* Check the kinematic range. */
                return 0.;
        } else {
                return physics->dcs_photonuclear(element->Z, element->A,
                    physics->mass, K, q) * 1E+03 * AVOGADRO_NUMBER *
                    (physics->mass + K) / element->A;
        }
}

/**
 * Utility function for checking the consistency of the Photonuclear model.
 *
 * @param K The kinetic energy.
 * @param q The kinetic energy lost to the photon.
 * @return `0` if the model is valid.
 *
 * Check for the limit of the PDG model. Below this kinetic transfer a
 * tabulation is used, which we don't do. Therefore we set the cross-section to
 * zero below 1 GeV where it deviates from the model. In addition, for x < 2E-03
 * the model is unstable and can lead to osccilations and negative cross-section
 * values.
 */
int dcs_photonuclear_check(double K, double q)
{
        return (q < 1.) || (q < 2E-03 * K);
}

/**
 * The ionisation differential cross section.
 *
 * @param Physics Handle for physics tables.
 * @param element The target atomic element.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The differential cross section in m^2/kg
 *
 * The differential cross section for ionisation is computed following
 * Salvat et al., NIMB316 (2013) 144-159, considering only close interactions
 * for DELs. In addition a radiative correction is applied according to
 * Sokalski et al., Phys.Rev.D64 (2001) 074015 (MUM).
 */
double dcs_ionisation(const struct pumas_physics * physics,
    const struct atomic_element * element, double K, double q)
{
        const double Z = element->Z;
        const double A = element->A;
        const double P2 = K * (K + 2. * physics->mass);
        const double E = K + physics->mass;
        const double Wmax = 2. * ELECTRON_MASS * P2 /
            (physics->mass * physics->mass +
                                ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
        if ((Wmax < physics->cutoff * K) || (q > Wmax)) return 0.;
        const double Wmin = 0.62 * element->I;
        if (q <= Wmin) return 0.;

        /* Close interactions for Q >> atomic binding energies. */
        const double a0 = 0.5 / P2;
        const double a1 = -1. / Wmax;
        const double a2 = E * E / P2;
        const double cs =
            1.535336E-05 * E * Z / A * (a0 + 1. / q * (a1 + a2 / q));

        /* Radiative correction. */
        double Delta = 0.;
        const double m1 = physics->mass - ELECTRON_MASS;
        if (K >= 0.5 * m1 * m1 / ELECTRON_MASS) {
                const double L1 = log(1. + 2. * q / ELECTRON_MASS);
                Delta = 1.16141E-03 * L1 *
                    (log(4. * E * (E - q) / (physics->mass * physics->mass)) -
                            L1);
        }

        return cs * (1. + Delta);
}

/**
 * The analytical form for the partial integral of the ionisation DCS.
 *
 * @param Physics Handle for physics tables.
 * @param mode    Flag to select the integration mode.
 * @param element The target atomic element.
 * @param K       The projectile initial kinetic energy.
 * @param xlow    The lower bound of the fractional energy transfer.
 * @return The integrated dcs, in m^2/kg or the energy loss in GeV m^2/kg.
 *
 * The parameter *mode* controls the integration mode as following. If *mode*
 * is `0` the restricted cross section is computed, for energy losses larger
 * than *xlow*. Else, the restricted energy loss is computed, for losses
 * larger than *xlow*
 */
double dcs_ionisation_integrate(const struct pumas_physics * physics, int mode,
    const struct atomic_element * element, double K, double xlow)
{
        const double Z = element->Z;
        const double A = element->A;
        const double P2 = K * (K + 2. * physics->mass);
        const double E = K + physics->mass;
        const double Wmax = 2. * ELECTRON_MASS * P2 /
            (physics->mass * physics->mass +
                                ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
        if (Wmax < physics->cutoff * K) return 0.;
        double Wmin = 0.62 * element->I;
        const double qlow = K * xlow;
        if (qlow >= Wmin) Wmin = qlow;

        /* Check the bounds. */
        if (Wmax <= Wmin) return 0.;

        /* Close interactions for Q >> atomic binding energies. */
        const double a0 = 0.5 / P2;
        const double a1 = -1. / Wmax;
        const double a2 = E * E / P2;

        double I;
        if (mode == 0) {
                I = a0 * (Wmax - Wmin) + a1 * log(Wmax / Wmin) +
                    a2 * (1. / Wmin - 1. / Wmax);
        } else {
                I = 0.5 * a0 * (Wmax * Wmax - Wmin * Wmin) +
                    a1 * (Wmax - Wmin) + a2 * log(Wmax / Wmin);
        }

        return 1.535336E-05 * Z / A * I;
}

double dcs_ionisation_randomise(const struct pumas_physics * physics,
    struct pumas_context * context, const struct atomic_element * element,
    double K, double xlow)
{
        const double P2 = K * (K + 2. * physics->mass);
        const double E = K + physics->mass;
        const double Wmax = 2. * ELECTRON_MASS * P2 /
            (physics->mass * physics->mass +
                                ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
        if (Wmax < physics->cutoff * K) return K;
        double Wmin = 0.62 * element->I;
        const double qlow = K * xlow;
        if (qlow >= Wmin) Wmin = qlow;

        /* Close interactions for Q >> atomic binding energies. */
        const double a0 = 0.5 / P2;
        const double a1 = -1. / Wmax;
        const double a2 = E * E / P2;

        const double p0 = a0 * (Wmax - Wmin);
        const double p2 = a2 * (1. / Wmin - 1. / Wmax);
        double q;
        for (;;) {
                /* Inverse sampling using an enveloppe. */
                const double z = context->random(context) * (p0 + p2);
                if (z <= p0) {
                        q = (Wmax - Wmin) * context->random(context);
                } else {
                        q = Wmin /
                            (1. -
                                context->random(context) * (1. - Wmin / Wmax));
                }

                /* Rejection sampling. */
                const double r0 = a0 + a2 / (q * q);
                const double r1 = r0 + a1 / q;
                if (context->random(context) * r0 < r1) break;
        }

        return K - q;
}

/**
 * Encapsulation for the evaluation of DCS.
 *
 * @param Physics  Handle for physics tables.
 * @param context  The simulation context.
 * @param dcs_func The DCS function to evaluate.
 * @param element  The target atomic element.
 * @param K        The projectile initial kinetic energy.
 * @param q        The transfered energy.
 * @return The DCS value or `0`.
 *
 * This routine encapsulate the evaluation of DCS during the MC. It takes care
 * of checking whether an approximate model can be used or not. In addition it
 * applies a Jacobian weight factor for changing from nu = q / E to x = q / K.
 */
double dcs_evaluate(const struct pumas_physics * physics,
    struct pumas_context * context, dcs_function_t * dcs_func,
    const struct atomic_element * element, double K, double q)
{
        /* Compute the Jacobian factor. */
        const double wj = K / (K + physics->mass);

        /* Check if the process has a valid tabulated model. */
        if (dcs_func == dcs_ionisation)
                return dcs_func(physics, element, K, q) * wj;
        else if ((dcs_func == dcs_photonuclear) && dcs_photonuclear_check(K, q))
                return 0.;

        /* Check if the exact computation should be used. */
        const double min_k = physics ?
            *table_get_K(physics, physics->dcs_model_offset) :
            DCS_MODEL_MIN_KINETIC;
        if (!physics || (K <= min_k) || (q < physics->cutoff * K) ||
            (q > DCS_MODEL_MAX_FRACTION * K))
                return dcs_func(physics, element, K, q) * wj;

        /* Get the coefficients of the polynomial approximation. */
        const int index = dcs_get_index(dcs_func);
        const int offset = physics->dcs_model_offset;
        const int n = DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1;
        const float * const coeff =
            table_get_dcs_coeff(physics, element, index, 0);
        const int imax = physics->n_energies - 1;
        const float *c0, *c1;
        int i0;
        float h1;
        const double Kmax = *table_get_K(physics, imax);
        if ((K >= Kmax) || ((i0 = table_index(physics, context,
                                 table_get_K(physics, 0), K)) >= imax)) {
                /*  Rescale and use the last tabulated value. */
                h1 = 0.;
                c0 = c1 = coeff + (imax - offset) * (n + DCS_SAMPLING_N);
        } else {
                const float K0 = (float)(*table_get_K(physics, i0));
                const float K1 = (float)(*table_get_K(physics, i0 + 1));
                h1 = (float)((K - K0) / (K1 - K0));
                c0 = coeff + (i0 - offset) * (n + DCS_SAMPLING_N);
                c1 = c0 + n + DCS_SAMPLING_N;
        }

        /* Compute the polynomial approximations. */
        const float nu = (float)(q / K);
        float dcs0 = 0., dcs1 = 0.;
        float u = 1.;
        float lx = logf(nu);
        int i;
        for (i = 0; i < DCS_MODEL_ORDER_P + 1; i++) {
                dcs0 += c0[i] * u;
                dcs1 += c1[i] * u;
                u *= lx;
        }
        lx = logf(1.f - nu);
        u = lx;
        for (i = DCS_MODEL_ORDER_P + 1; i < n; i++) {
                dcs0 += c0[i] * u;
                dcs1 += c1[i] * u;
                u *= lx;
        }
        dcs0 = expf(dcs0);
        dcs1 = (h1 > 0.f) ? expf(dcs1) : 0.f;

        /* Return the interpolated result. */
        return (1. - h1) * dcs0 + h1 * dcs1;
}

/**
 * Fit the DCS model.
 *
 * @param m The number of constraints.
 * @param n The number of parameters.
 * @param x The fractional energy transfer.
 * @param y The DCS values.
 * @param w The fit weights.
 * @param c The model coefficients.
 *
 * This routines fits the log of the DCS with a polynomial model in log(nu)
 * and log(1-nu).
 */
void dcs_model_fit(int m, int n, const double * x, const double * y,
    const double * w, double * c)
{
        /* Memory management. */
        static int size = 0;
        static double * memory = NULL;
        const int size_ = m * n;
        if (size != size_) {
                memory =
                    reallocate(memory, (n * (m + n + 2) + m) * sizeof(double));
                size = size_;
        }
        if ((m <= 0) || (n <= 0)) return;

        double * A = memory;
        double * B = A + m * n;
        double * W = B + m;
        double * V = W + n;
        double * work = V + n * n;

        /* Prepare the linear system. */
        int i, j, ij;
        for (i = 0, ij = 0; i < m; i++) {
                double xi = 1.;
                double lxi = log(x[i]);
                for (j = 0; j < DCS_MODEL_ORDER_P + 1; j++, ij++) {
                        A[ij] = w[i] * xi;
                        xi *= lxi;
                }
                lxi = log(1. - x[i]);
                xi = lxi;
                for (j = 0; j < DCS_MODEL_ORDER_Q; j++, ij++) {
                        A[ij] = w[i] * xi;
                        xi *= lxi;
                }
                B[i] = (w[i] > 0.) ? log(y[i]) * w[i] : 0.;
        }

        /* Solve using SVD. */
        if (math_svd(m, n, A, W, V, work) != 0) return;
        math_svdsol(m, n, B, A, W, V, work, c);
}

/*
 * Low level routines: sampling of the polar angle in a DEL.
 *
 * All the sampling routines below are implementations of the methods
 * described in the Geant4 Physics Reference Manual (PRM).
 */
/**
 * Get a polar angle distrubution given a process index.
 *
 * @param process The process index.
 * @return The polar angle function.
 *
 * `index = 0` corresponds to Bremsstrahlung, `index = 1` to Pair Production,
 * `index = 2` to Photonuclear interactions and `index = 3` to Ionisation.
 */
polar_function_t * polar_get(int process)
{
        polar_function_t * const polar_func[] = { polar_bremsstrahlung,
                polar_pair_production, polar_photonuclear, polar_ionisation };
        return polar_func[process];
}

/**
 * Sample the polar angle in a Bremsstrahlung event.
 *
 * @param Physics Handle for physics tables.
 * @param context The simulation context.
 * @param ki      The initial kinetic energy.
 * @param kf      The final kinetic energy.
 * @return The cosine of the polar angle.
 */
double polar_bremsstrahlung(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf)
{
        const double e = ki + physics->mass;
        const double q = ki - kf;
        double rmax2 = e / q - 1.;
        rmax2 = (rmax2 > 1.) ? 1. : rmax2 * rmax2;
        double r = context->random(context) * rmax2 / (1. + rmax2);
        r = sqrt(r / (1. - r));

        return cos(r * physics->mass * q / (e * (kf + physics->mass)));
}

/**
 * Sample the polar angle in a Pair Production event.
 *
 * @param Physics Handle for physics tables.
 * @param context The simulation context.
 * @param ki      The initial kinetic energy.
 * @param kf      The final kinetic energy.
 * @return The cosine of the polar angle.
 *
 * The polar angle is sampled assuming a virtual Bremsstrahlung event.
 */
double polar_pair_production(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf)
{
        return polar_bremsstrahlung(physics, context, ki, kf);
}

/**
 * Sample the polar angle in a photonuclear event.
 *
 * @param Physics Handle for physics tables.
 * @param context The simulation context.
 * @param ki      The initial kinetic energy.
 * @param kf      The final kinetic energy.
 * @return The cosine of the polar angle.
 */
double polar_photonuclear(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf)
{
        const double q = ki - kf;
        const double e = ki + physics->mass;
        const double y = q / e;
        const double tmin = physics->mass * physics->mass * y * y / (1. - y);
        const double tmax = 1.87914 * q;
        const double m02 = 0.4;
        const double q2 = q * q;
        const double t1 = (q2 < m02) ? q2 : m02;

        const double p = context->random(context);
        const double r = tmax * (tmin + t1) / (tmin * (tmax + t1));
        const double tp = tmax * t1 / ((tmax + t1) * pow(r, p) - tmax);
        const double ct = 1. -
            (tp - tmin) / (2. * (e * (kf + physics->mass) -
                                    physics->mass * physics->mass) -
                              tmin);

        if (ct < 0.)
                return 0.;
        else if (ct > 1.)
                return 1.;
        else
                return ct;
}

/**
 * Sample the polar angle in an ionisation event.
 *
 * @param Physics Handle for physics tables.
 * @param context The simulation context.
 * @param ki      The initial kinetic energy.
 * @param kf      The final kinetic energy.
 * @return The cosine of the polar angle.
 *
 * The polar angle is set from energy-momentum conservation assuming that the
 * electron is initially at rest. See for example appendix A of Fernandez-Varea
 * et al., NIMB 229 (2005) 185-218.
 */
double polar_ionisation(const struct pumas_physics * physics,
    struct pumas_context * context, double ki, double kf)
{
        if (kf > ki) {
                const double tmp = ki;
                ki = kf, kf = tmp;
        }
        const double p02 = ki * (ki + 2. * physics->mass);
        const double p12 = kf * (kf + 2. * physics->mass);
        const double ke = ki - kf;
        const double pe2 = ke * (ke + 2. * ELECTRON_MASS);
        return 0.5 * (p02 + p12 - pe2) / sqrt(p02 * p12);
}

/*
 * Low level routines: various algorithms for handling some specific
 * mathematic problems.
 */
/**
 * Root solver for a function of a scalar variable.
 *
 * @param Physics Handle for physics tables.
 * @param f       The objective function to resolve.
 * @param xa      The lower bound of the search interval.
 * @param xb      The upper bound of the search interval.
 * @param fa_p    The initial value at *a* if already computed.
 * @param fb_p    The initial value at *b* if already computed.
 * @param xtol    The absolute tolerance on the root value.
 * @param rtol    The absolute tolerance on the root value.
 * @param params  A handle for passing additional parameters to the
 *                objective function.
 * @param x0      An estimate of the root over `[xa; xb]`.
 * @return On success `0` is returned. Otherwise a negative number.
 *
 * The root is searched for over `[xa; xb]` using Ridder's method
 * (https://en.wikipedia.org/wiki/Ridders%27_method). If initial values of the
 * objective have already been computed they can be passed over as *fa_p* or
 * *fb_p*. Otherwise these values must be `NULL`. The relative tolerance is
 * relative to min(xa, xb).
 */
int math_find_root(
    double (*f)(const struct pumas_physics * physics, double x, void * params),
    const struct pumas_physics * physics, double xa, double xb,
    const double * fa_p, const double * fb_p, double xtol, double rtol,
    int max_iter, void * params, double * x0)
{
        /*  Check the initial values. */
        double fa, fb;
        if (fa_p == NULL)
                fa = (*f)(physics, xa, params);
        else
                fa = *fa_p;
        if (fb_p == NULL)
                fb = (*f)(physics, xb, params);
        else
                fb = *fb_p;
        if (fa * fb > 0) {
                *x0 = 0.;
                return -1;
        }
        if (fa == 0) {
                *x0 = xa;
                return 0;
        }
        if (fb == 0) {
                *x0 = xb;
                return 0;
        }

        /* Set the tolerance for the root finding*/
        const double tol =
            xtol + rtol * ((fabs(xa) < fabs(xb)) ? fabs(xa) : fabs(xb));

        /* Do the bracketing using Ridder's update rule. */
        double xn = 0.;
        int i = 0;
        for (i = 0; i < max_iter; i++) {
                double dm = 0.5 * (xb - xa);
                const double xm = xa + dm;
                const double fm = (*f)(physics, xm, params);
                double sgn = (fb > fa) ? 1. : -1.;
                double dn = sgn * dm * fm / sqrt(fm * fm - fa * fb);
                sgn = (dn > 0.) ? 1. : -1.;
                dn = fabs(dn);
                dm = fabs(dm) - 0.5 * tol;
                if (dn < dm) dm = dn;
                xn = xm - sgn * dm;
                const double fn = (*f)(physics, xn, params);
                if (fn * fm < 0.0) {
                        xa = xn;
                        fa = fn;
                        xb = xm;
                        fb = fm;
                } else if (fn * fa < 0.0) {
                        xb = xn;
                        fb = fn;
                } else {
                        xa = xn;
                        fa = fn;
                }
                if (fn == 0.0 || fabs(xb - xa) < tol) {
                        /*  A valid bracketing was found*/
                        *x0 = xn;
                        return 0;
                }
        }

        /* The maximum number of iterations was reached*/
        *x0 = xn;
        return -2;
}

/**
 * Iterator function for integration with a Gaussian quadrature.
 *
 * @param n  The minimum number of integration points.
 * @param p1 The integration range lower bound or the sampling point.
 * @param p2 The integration range upper bound or the sampling weight.
 * @return The number of points used for the integration or a return code.
 *
 * If *n* is strictly positive the iterator is initialised for a new
 * integration with at least *n* points. The integration range is `[p1, p2]`
 * and at return the number of integration points is returned.
 *
 * Otherwise, the next integration point is provided in *p1* and *p2* is
 * filled with the corresponding weight. If all integration steps have been
 * done `1` is returned, otherwise `0`.
 */
int math_gauss_quad(int n, double * p1, double * p2)
{
/*
 * Coefficients for the Gaussian quadrature from:
 * https://pomax.github.io/bezierinfo/legendre-gauss.html.
 */
#define N_GQ 6
        const double xGQ[N_GQ] = { 0.03376524, 0.16939531, 0.38069041,
                0.61930959, 0.83060469, 0.96623476 };
        const double wGQ[N_GQ] = { 0.08566225, 0.18038079, 0.23395697,
                0.23395697, 0.18038079, 0.08566225 };

        /* Initialisation step. */
        static int i, j, n_itv;
        static double h;
        static double x0;
        if (n > 0) {
                n_itv = (n + N_GQ - 1) / N_GQ;
                i = j = 0;
                h = (*p2 - (x0 = *p1)) / n_itv;
                return N_GQ * n_itv;
        }

        /* Iteration step. */
        if (i == n_itv) return 1;

        *p1 = x0 + xGQ[j] * h;
        *p2 = wGQ[j] * h;

        if (++j == N_GQ) {
                i++;
                j = 0;
                x0 += h;
        }
        return 0;

#undef N_GQ
}

/*
 * Low level routines: Pseudo inverse from Singular Value Decomposition (SVD).
 *
 * The routines below are adapted from SLALIB/C (Copyright P.T.Wallace) which
 * is distributed under the LGPL license.
 */
/**
 * Singular value decomposition of a matrix as A = U x W x V^t.
 *
 * @param m     The number of rows in A.
 * @param n     The number of columns in A.
 * @param a     The A matrix to decompose or the U matrix at output.
 * @param w     The W matrix.
 * @param v     The V matrix, **not** its tranpose.
 * @param work  A length *n* temporary array.
 * @return On success, `0` is returned. Otherwise an error occured.
 *
 * This routine expresses a given matrix *A* as the product of three matrices
 * *U*, *W*, *V*, as `A = U x W x V^t`, with U, V orthogonal and W diagonal.
 * The algorithm is an adaptation of the routine SVD in the EISPACK library
 * (Garbow et al 1977, Eispack guide extension, Springer Verlag), which is a
 * Fortran 66 implementation of the Algol routine SVD of Wilkinson & Reinsch
 * 1971 (Handbook for Automatic Computation, vol 2, Ed Bauer et al, Springer
 * Verlag). For the non-specialist, probably the clearest general account of
 * the use of SVD in least squares problems is given in Numerical Recipes
 * (Press et al 1986, Cambridge University Press).
 *
 * **Warning** the *A* matrix is ovewritten at output.
 */
int math_svd(int m, int n, double * a, double * w, double * v, double * work)
{
#define dsign(A, B) ((B) < 0.0 ? -(A) : (A))
#define gmax(A, B) ((A) > (B) ? (A) : (B))
/* Maximum number of iterations in QR phase. */
#define ITMAX 30

        int i, k, l = 0, j, k1, its, l1 = 0, i1, cancel;
        double g, scale, an, s, x, f, h, cn, c, y, z;
        double *ai, *aj, *ak;
        double *vi, *vj, *vk;

        /* Check that the matrix is the right size. */
        if (m < n) return -1;

        /* Householder reduction to bidiagonal form. */
        int jstat = 0;
        g = 0.0;
        scale = 0.0;
        an = 0.0;
        for (i = 0, ai = a; i < n; i++, ai += n) {
                l = i + 1;
                work[i] = scale * g;
                g = 0.0;
                s = 0.0;
                scale = 0.0;
                if (i < m) {
                        for (k = i, ak = ai; k < m; k++, ak += n) {
                                scale += fabs(ak[i]);
                        }
                        if (scale != 0.0) {
                                for (k = i, ak = ai; k < m; k++, ak += n) {
                                        x = ak[i] / scale;
                                        ak[i] = x;
                                        s += x * x;
                                }
                                f = ai[i];
                                g = -dsign(sqrt(s), f);
                                h = f * g - s;
                                ai[i] = f - g;
                                if (i != n - 1) {
                                        for (j = l; j < n; j++) {
                                                s = 0.0;
                                                for (k = i, ak = ai; k < m;
                                                     k++, ak += n) {
                                                        s += ak[i] * ak[j];
                                                }
                                                f = s / h;
                                                for (k = i, ak = ai; k < m;
                                                     k++, ak += n) {
                                                        ak[j] += f * ak[i];
                                                }
                                        }
                                }
                                for (k = i, ak = ai; k < m; k++, ak += n) {
                                        ak[i] *= scale;
                                }
                        }
                }
                w[i] = scale * g;
                g = 0.0;
                s = 0.0;
                scale = 0.0;
                if (i < m && i != n - 1) {
                        for (k = l; k < n; k++) {
                                scale += fabs(ai[k]);
                        }
                        if (scale != 0.0) {
                                for (k = l; k < n; k++) {
                                        x = ai[k] / scale;
                                        ai[k] = x;
                                        s += x * x;
                                }
                                f = ai[l];
                                g = -dsign(sqrt(s), f);
                                h = f * g - s;
                                ai[l] = f - g;
                                for (k = l; k < n; k++) {
                                        work[k] = ai[k] / h;
                                }
                                if (i != m - 1) {
                                        for (j = l, aj = a + l * n; j < m;
                                             j++, aj += n) {
                                                s = 0.0;
                                                for (k = l; k < n; k++) {
                                                        s += aj[k] * ai[k];
                                                }
                                                for (k = l; k < n; k++) {
                                                        aj[k] += s * work[k];
                                                }
                                        }
                                }
                                for (k = l; k < n; k++) {
                                        ai[k] *= scale;
                                }
                        }
                }

                /* Overestimate of largest column norm for convergence test. */
                cn = fabs(w[i]) + fabs(work[i]);
                an = gmax(an, cn);
        }

        /* Accumulation of right-hand transformations. */
        for (i = n - 1, ai = a + (n - 1) * n, vi = v + (n - 1) * n; i >= 0;
             i--, ai -= n, vi -= n) {
                if (i != n - 1) {
                        if (g != 0.0) {
                                for (j = l, vj = v + l * n; j < n;
                                     j++, vj += n) {
                                        vj[i] = (ai[j] / ai[l]) / g;
                                }
                                for (j = l; j < n; j++) {
                                        s = 0.0;
                                        for (k = l, vk = v + l * n; k < n;
                                             k++, vk += n) {
                                                s += ai[k] * vk[j];
                                        }
                                        for (k = l, vk = v + l * n; k < n;
                                             k++, vk += n) {
                                                vk[j] += s * vk[i];
                                        }
                                }
                        }
                        for (j = l, vj = v + l * n; j < n; j++, vj += n) {
                                vi[j] = 0.0;
                                vj[i] = 0.0;
                        }
                }
                vi[i] = 1.0;
                g = work[i];
                l = i;
        }

        /* Accumulation of left-hand transformations. */
        for (i = n - 1, ai = a + i * n; i >= 0; i--, ai -= n) {
                l = i + 1;
                g = w[i];
                if (i != n - 1) {
                        for (j = l; j < n; j++) {
                                ai[j] = 0.0;
                        }
                }
                if (g != 0.0) {
                        if (i != n - 1) {
                                for (j = l; j < n; j++) {
                                        s = 0.0;
                                        for (k = l, ak = a + l * n; k < m;
                                             k++, ak += n) {
                                                s += ak[i] * ak[j];
                                        }
                                        f = (s / ai[i]) / g;
                                        for (k = i, ak = a + i * n; k < m;
                                             k++, ak += n) {
                                                ak[j] += f * ak[i];
                                        }
                                }
                        }
                        for (j = i, aj = ai; j < m; j++, aj += n) {
                                aj[i] /= g;
                        }
                } else {
                        for (j = i, aj = ai; j < m; j++, aj += n) {
                                aj[i] = 0.0;
                        }
                }
                ai[i] += 1.0;
        }

        /* Diagonalization of the bidiagonal form. */
        for (k = n - 1; k >= 0; k--) {
                k1 = k - 1;

                /* Iterate until converged. */
                for (its = 1; its <= ITMAX; its++) {

                        /* Test for splitting into submatrices. */
                        cancel = 1;
                        for (l = k; l >= 0; l--) {
                                l1 = l - 1;
                                if (an + fabs(work[l]) == an) {
                                        cancel = 0;
                                        break;
                                }
                                /* (Following never attempted for l=0 because
                                 * work[0] is zero).
                                 */
                                if (an + fabs(w[l1]) == an) {
                                        break;
                                }
                        }

                        /* Cancellation of work[l] if l>0. */
                        if (cancel) {
                                c = 0.0;
                                s = 1.0;
                                for (i = l; i <= k; i++) {
                                        f = s * work[i];
                                        if (an + fabs(f) == an) {
                                                break;
                                        }
                                        g = w[i];
                                        h = math_rms(f, g);
                                        w[i] = h;
                                        c = g / h;
                                        s = -f / h;
                                        for (j = 0, aj = a; j < m;
                                             j++, aj += n) {
                                                y = aj[l1];
                                                z = aj[i];
                                                aj[l1] = y * c + z * s;
                                                aj[i] = -y * s + z * c;
                                        }
                                }
                        }

                        /* Converged? */
                        z = w[k];
                        if (l == k) {

                                /* Yes: ensure singular values non-negative. */
                                if (z < 0.0) {
                                        w[k] = -z;
                                        for (j = 0, vj = v; j < n;
                                             j++, vj += n) {
                                                vj[k] = -vj[k];
                                        }
                                }

                                /* Stop iterating. */
                                break;

                        } else {

                                /* Not converged yet: set status if iteration
                                 * limit reached.
                                 */
                                if (its >= ITMAX) {
                                        jstat = k + 1;
                                }

                                /* Shift from bottom 2 x 2 minor. */
                                x = w[l];
                                y = w[k1];
                                g = work[k1];
                                h = work[k];
                                f = ((y - z) * (y + z) + (g - h) * (g + h)) /
                                    (2.0 * h * y);
                                g = (fabs(f) <= 1e15) ? math_rms(f, 1.0) :
                                                        fabs(f);
                                f = ((x - z) * (x + z) +
                                        h * (y / (f + dsign(g, f)) - h)) /
                                    x;

                                /* Next QR transformation. */
                                c = 1.0;
                                s = 1.0;
                                for (i1 = l; i1 <= k1; i1++) {
                                        i = i1 + 1;
                                        g = work[i];
                                        y = w[i];
                                        h = s * g;
                                        g = c * g;
                                        z = math_rms(f, h);
                                        work[i1] = z;
                                        if (z != 0.0) {
                                                c = f / z;
                                                s = h / z;
                                        } else {
                                                c = 1.0;
                                                s = 0.0;
                                        }
                                        f = x * c + g * s;
                                        g = -x * s + g * c;
                                        h = y * s;
                                        y = y * c;
                                        for (j = 0, vj = v; j < n;
                                             j++, vj += n) {
                                                x = vj[i1];
                                                z = vj[i];
                                                vj[i1] = x * c + z * s;
                                                vj[i] = -x * s + z * c;
                                        }
                                        z = math_rms(f, h);
                                        w[i1] = z;
                                        if (z != 0.0) {
                                                c = f / z;
                                                s = h / z;
                                        }
                                        f = c * g + s * y;
                                        x = -s * g + c * y;
                                        for (j = 0, aj = a; j < m;
                                             j++, aj += n) {
                                                y = aj[i1];
                                                z = aj[i];
                                                aj[i1] = y * c + z * s;
                                                aj[i] = -y * s + z * c;
                                        }
                                }
                                work[l] = 0.0;
                                work[k] = f;
                                w[k] = x;
                        }
                }
        }

        /* Return the status flag. */
        return jstat;

#undef dsign
#undef gmax
#undef ITMAX
}

/**
 * Compute sqrt(a*a+b*b) with protection against under/overflow.
 */
double math_rms(double a, double b)
{
        double wa = fabs(a);
        double wb = fabs(b);
        if (wa > wb) {
                const double tmp = wa;
                wa = wb;
                wb = tmp;
        }

        if (wb == 0.) return 0.;
        const double w = wa / wb;
        return wb * sqrt(1.0 + w * w);
}

/**
 * Solve a linear problem `A X = B` using a pseudo-inverse from SVD.
 *
 * @param m     The number of rows in A.
 * @param n     The number of columns in A.
 * @param a     The B array.
 * @param a     The U matrix from the SVD decomposition of A.
 * @param w     The W matrix from the SVD decomposition of A.
 * @param v     The V matrix from the SVD decomposition of A.
 * @param work  A length *n* temporary array.
 * @param X     The computed pseudo inverse.
 */
void math_svdsol(int m, int n, double * b, double * u, double * w, double * v,
    double * work, double * x)
{
        /* Calculate [diag(1/Wj)] . U^T . B (or zero for zero Wj). */
        int j, i, jj;
        const double * ui;
        for (j = 0; j < n; j++) {
                double s = 0.0;
                if (w[j] != 0.0) {
                        for (i = 0, ui = u; i < m; i++, ui += n) {
                                s += ui[j] * b[i];
                        }
                        s /= w[j];
                }
                work[j] = s;
        }

        /* Multiply by matrix V to get result. */
        const double * vj;
        for (j = 0, vj = v; j < n; j++, vj += n) {
                double s = 0.0;
                for (jj = 0; jj < n; jj++) {
                        s += vj[jj] * work[jj];
                }
                x[j] = s;
        }
}

/** Dilogarithm for real valued arguments
 *
 * Ref: CERNLIB RDILOG function (C332)
 * http://cds.cern.ch/record/2050865
 */
static double math_dilog(double x)
{
        const double C[20] = {
             0.42996693560813697,  0.40975987533077105, -0.01858843665014592,
             0.00145751084062268, -0.00014304184442340,  0.00001588415541880,
            -0.00000190784959387,  0.00000024195180854, -0.00000003193341274,
             0.00000000434545063, -0.00000000060578480,  0.00000000008612098,
            -0.00000000001244332,  0.00000000000182256, -0.00000000000027007,
             0.00000000000004042, -0.00000000000000610,  0.00000000000000093,
            -0.00000000000000014,  0.00000000000000002};

        const double PI3 = M_PI * M_PI / 3.;
        const double PI6 = M_PI * M_PI / 6.;
        const double PI12 = M_PI * M_PI / 12.;

        if (x == 10) {
                return PI6;
        } else if (x == -1.) {
                return -PI12;
        } else {
                double T = -x;
                double A, S, Y;

                if (T <= -2) {
                        Y = -1. / (1. + T);
                        S = 1.;
                        A = -PI3 + 0.5 * (pow(log(-T), 2.) -
                            pow(log(1. + 1. / T), 2.));
                } else if (T < -1.) {
                        Y = -1. - T;
                        S = -1.;
                        A = log(-T);
                        A = -PI6 + A * (A + log(1. + 1. / T));
                } else if (T <= -0.5) {
                        Y = -(1. + T) / T;
                        S = 1.;
                        A = log(-T);
                        A = -PI6 + A * (-0.5 * A + log(1. + T));
                } else if (T < 0.) {
                        Y = -T / (1. + T);
                        S = -1.;
                        A = 0.5 * pow(log(1 + T), 2.);
                } else if (T <= 1.) {
                        Y = T;
                        S = 1.;
                        A = 0.;
                } else {
                        Y = 1. / T;
                        S = -1.;
                        A = PI6 + 0.5 * pow(log(T), 2);
                }

                const double H = Y + Y - 1;
                const double ALFA = H + H;
                double B0, B1 = 0., B2 = 0.;
                int i;
                for (i = 19; i >= 0; i--) {
                        B0 = C[i] + ALFA * B1 - B2;
                        B2 = B1;
                        B1 = B0;
                }

                return -(S * (B0 - H * B2) + A);
        }
}

/**
 * Routines for computing energy loss tabulations.
 */

/** Tabulation mode initialisation, without loading the energy loss tables. */
enum pumas_return pumas_physics_create_tabulation(
    struct pumas_physics ** physics, enum pumas_particle particle,
    const char * mdf_path, const struct pumas_physics_settings * settings)
{
        return _initialise(physics, particle, mdf_path, NULL, 1, settings);
}

/* The density effect for the electronic energy loss. */
static void electronic_density_effect(const struct pumas_physics * physics,
    const struct pumas_physics_material * m, double kinetic, double * delta)
{
        const double c = 2. * log(10.);
        const double r = kinetic / physics->mass;
        const double x = 0.5 * log10(r * (r + 2.));
        const struct pumas_physics_density_effect * const d =
            &m->density_effect;
        if (x < d->x0)
                *delta = d->delta0 > 0. ?
                    d->delta0 * pow(10., 2. * (x - d->x0)) :
                    0.;
        else if (x < d->x1)
                *delta = c * x - d->Cbar + d->a * pow(d->x1 - x, d->k);
        else
                *delta = c * x - d->Cbar;
}

/* The average energy loss from atomic electrons. */
static double electronic_energy_loss(const struct pumas_physics * physics,
    const struct pumas_physics_material * m, double kinetic, double * delta)
{
        /* Kinematic factors. */
        const double E = kinetic + physics->mass;
        const double P2 = kinetic * (kinetic + 2. * physics->mass);
        const double beta2 = P2 / (E * E);

        /* Electronic Bremsstrahlung correction. */
        const double r = ELECTRON_MASS / physics->mass;
        const double Qmax =
            2. * r * P2 / (physics->mass * (1. + r * r) + 2. * r * E);
        const double lQ = log(2. * Qmax / ELECTRON_MASS);
        const double Delta =
            5.8070487E-04 * (log(2. * E / physics->mass) - lQ / 3.) * lQ * lQ;

        /* Density effect. */
        electronic_density_effect(physics, m, kinetic, delta);

        /* Bethe-Bloch equation. */
        return 0.307075E-04 * physics->material_ZoA[m->index] *
            (0.5 / beta2 * (log(2. * ELECTRON_MASS * P2 * Qmax /
                                (physics->mass * physics->mass * m->I * m->I)) -
                               *delta) -
                   1. + 0.125 * Qmax * Qmax / P2 + Delta);
}

/* Container for atomic element tabulation data */
struct tabulation_element {
        /* The API proxy */
        struct pumas_physics_element api;
        /* Placeholder for tabulation data */
        double data[];
};

static void tabulate_element(struct pumas_physics * physics,
    struct tabulation_element * data, int n_energies, double * kinetic)
{
        const struct atomic_element * element =
            physics->element[data->api.index];

        /* Loop over the kinetic energy values. */
        double * v;
        int ik, ip;
        for (ik = 0, v = data->data; ik < n_energies; ik++) {
                const double k = kinetic[ik];
                double x = 1E-06 / k;
                if (x < 1E-05) x = 1E-05;
                const int n =
                    (int)(-1E+02 * log10(x)); /* 100 pts per decade. */
                for (ip = 0; ip < N_DEL_PROCESSES - 1; ip++, v++)
                        *v = compute_dcs_integral(
                            physics, 1, element, k, dcs_get(ip), x, n);
        }
}

/*
 * Create a new energy loss table for an element and add it to the stack of
 * temporary data.
 */
static struct pumas_physics_element * tabulation_element_create(
    struct pumas_physics_tabulation_data * data, int element)
{
        /* Allocate memory for the new element. */
        struct pumas_physics_element * e =
            allocate(sizeof(struct tabulation_element) +
                3 * data->n_energies * sizeof(double));
        if (e == NULL) return NULL;
        e->index = element;
        e->fraction = 0.;

        /* Add the element's data on top of the stack. */
        if (data->elements != NULL) data->elements->next = e;
        e->prev = data->elements;
        e->next = NULL;
        data->elements = e;

        return e;
}

/*
 * Get the energy loss table for an element from the temporary data and
 * put it on top of the stack.
 */
static struct pumas_physics_element * tabulation_element_get(
    struct pumas_physics_tabulation_data * data, int element)
{
        struct pumas_physics_element * e;
        for (e = data->elements; e != NULL; e = e->prev) {
                if (e->index == element) {
                        struct pumas_physics_element * next = e->next;
                        if (next != NULL) {
                                /* Put the element on top of the stack. */
                                struct pumas_physics_element * prev = e->prev;
                                if (prev != NULL) prev->next = next;
                                next->prev = prev;
                                data->elements->next = e;
                                e->prev = data->elements;
                                e->next = NULL;
                                data->elements = e;
                        }
                        return e;
                }
        }

        /* The element wasn't found, return `NULL`. */
        return NULL;
}

/**
 * Tabulate the energy loss for the given material and kinetic energies.
 */
enum pumas_return pumas_physics_tabulate(
    struct pumas_physics * physics, struct pumas_physics_tabulation_data * data)
{
        ERROR_INITIALISE(pumas_physics_create);

        /* Check the material index */
        struct pumas_physics_material * m = &data->material;
        if ((m->index < 0) ||
            (m->index >= physics->n_materials - physics->n_composites)) {
                return ERROR_FORMAT(PUMAS_RETURN_INDEX_ERROR,
                    "invalid material index [%d]", m->index);
        }

        /* Set the energy grid */
        if ((data->n_energies <= 0) || (data->energy == NULL)) {
                static double energy_[201];
                data->energy = energy_;

                if (physics->particle == PUMAS_PARTICLE_MUON) {
                        /* For muons the PDG energy grid is used by default. */
                        energy_[0] = 1.000E-03;
                        energy_[1] = 1.200E-03;
                        energy_[2] = 1.400E-03;
                        energy_[3] = 1.700E-03;
                        energy_[4] = 2.000E-03;
                        energy_[5] = 2.500E-03;
                        energy_[6] = 3.000E-03;
                        energy_[7] = 3.500E-03;
                        energy_[8] = 4.000E-03;
                        energy_[9] = 4.500E-03;
                        energy_[10] = 5.000E-03;
                        energy_[11] = 5.500E-03;
                        energy_[12] = 6.000E-03;
                        energy_[13] = 7.000E-03;
                        energy_[14] = 8.000E-03;
                        energy_[15] = 9.000E-03;

                        const int n_per_decade = 16;
                        const int n_decades = 8;
                        data->n_energies = (n_decades + 1) * n_per_decade + 1;

                        int i;
                        for (i = 1; i <= n_decades; i++) {
                                int j;
                                for (j = 0; j < n_per_decade; j++) {
                                        energy_[i * n_per_decade + j] = 10 *
                                            energy_[(i - 1) * n_per_decade + j];
                                }
                        }
                        energy_[(n_decades + 1) * n_per_decade] = 10 *
                            energy_[n_decades * n_per_decade];
                } else {
                        /* For taus a logarithmic energy grid is used by
                         * default.
                         */
                        data->n_energies = 201;
                        double emin = 1E+02, emax = 1E+12;
                        const double dlnk = log(emax / emin) /
                            (data->n_energies - 1);
                        int i;
                        for (i = 0; i < data->n_energies; i++) {
                                energy_[i] = emin * exp(dlnk * i);
                        }
                }
        }

        /* Compute the mean excitation energy, if not provided. */
        if (m->I <= 0.) {
                double lnI = 0., Z = 0.;
                struct material_component * component;
                int iel;
                for (iel = 0, component = physics->composition[m->index];
                     iel < physics->elements_in[m->index]; iel++, component++) {
                        struct atomic_element * e =
                            physics->element[component->element];
                        const double nZ = component->fraction * e->Z / e->A;
                        lnI += nZ * log(e->I);
                        Z += nZ;
                }
                m->I = exp(lnI / Z);
        }

        /* Compute the density effect coefficients, if not provided. */
        if (m->density_effect.a <= 0.) {
                struct pumas_physics_density_effect * const d =
                    &m->density_effect;

                /* Use the Sternheimer and Peierls recipee. */
                d->k = 3.;
                const double density = physics->material_density[m->index];
                const double hwp = 28.816E-09 *
                    sqrt(density * 1E-03 * physics->material_ZoA[m->index]);
                d->Cbar = 2. * log(m->I / hwp);
                if (m->state == PUMAS_PHYSICS_STATE_GAS) {
                        if (d->Cbar < 10.) {
                                d->x0 = 1.6, d->x1 = 4.;
                        } else if (d->Cbar < 10.5) {
                                d->x0 = 1.7, d->x1 = 4.;
                        } else if (d->Cbar < 11.) {
                                d->x0 = 1.8, d->x1 = 4.;
                        } else if (d->Cbar < 11.5) {
                                d->x0 = 1.9, d->x1 = 4.;
                        } else if (d->Cbar < 12.25) {
                                d->x0 = 2.0, d->x1 = 4.;
                        } else if (d->Cbar < 13.804) {
                                d->x0 = 2.0, d->x1 = 5.;
                        } else {
                                d->x0 = 0.326 * d->Cbar - 1.5;
                                d->x1 = 5.0;
                        }
                } else {
                        if (m->I < 100.E-09) {
                                d->x1 = 2.;
                                if (d->Cbar < 3.681)
                                        d->x0 = 0.2;
                                else
                                        d->x0 = 0.326 * d->Cbar - 1.;
                        } else {
                                d->x1 = 3.;
                                if (d->Cbar < 5.215)
                                        d->x0 = 0.2;
                                else
                                        d->x0 = 0.326 * d->Cbar - 1.5;
                        }
                }
                const double dx = d->x1 - d->x0;
                d->a = (d->Cbar - 2. * log(10.) * d->x0) / (dx * dx * dx);
        }

        /*
         * Tabulate the radiative energy losses for the constitutive atomic
         * elements, if not already done. Note that the element are also sorted
         * on top of the *data* stack.
         */
        struct material_component * component;
        int iel;
        for (iel = 0, component = physics->composition[m->index];
             iel < physics->elements_in[m->index]; iel++, component++) {
                /*
                 * Get the requested element's data and put them on top of
                 * the stack for further usage.
                 */
                struct pumas_physics_element * e =
                    tabulation_element_get(data, component->element);

                if (e == NULL) {
                        /* Create and tabulate the new element. */
                        e = tabulation_element_create(data, component->element);
                        if (e == NULL) return PUMAS_RETURN_MEMORY_ERROR;
                        tabulate_element(physics,
                            (struct tabulation_element *)e, data->n_energies,
                            data->energy);
                }

                /* Set the fraction in the current material. */
                e->fraction = component->fraction;
        }

        /* Check and open the output file. */
        int n_d = (data->outdir == NULL) ? 0 : strlen(data->outdir) + 1;
        int n_f = strlen(physics->dedx_filename[m->index]) + 1;
        char * path = reallocate(data->path, n_d + n_f);
        if (path == NULL) return PUMAS_RETURN_MEMORY_ERROR;
        data->path = path;
        if (n_d) {
                memcpy(path, data->outdir, n_d);
                path[n_d - 1] = '/';
        }
        memcpy(path + n_d, physics->dedx_filename[m->index], n_f);

        FILE * stream;
        if (data->overwrite == 0) {
                /* Check if the file already exists. */
                stream = fopen(path, "r");
                if (stream != NULL) {
                        fclose(stream);
                        return PUMAS_RETURN_IO_ERROR;
                }
        }
        stream = fopen(path, "w+");
        if (stream == NULL) return PUMAS_RETURN_PATH_ERROR;

        /* Print the header. */
        const char * type =
            (physics->particle == PUMAS_PARTICLE_MUON) ? "Muon" : "Tau";
        fprintf(stream, " Incident particle is a %s with M = %.5lf MeV\n", type,
            (double)(physics->mass * 1E+03));
        fprintf(stream, " Index = %d: %s\n", m->index,
            physics->material_name[m->index]);
        fprintf(stream, "      Absorber with <Z/A> = %.5lf\n",
            physics->material_ZoA[m->index]);
        fprintf(stream, " Sternheimer coef:  a     k=m_s   x_0    x_1    "
                        "I[eV]   Cbar  delta0\n");
        const struct pumas_physics_density_effect * const d =
            &m->density_effect;
        fprintf(stream, "                %7.4lf %7.4lf %7.4lf %7.4lf %6.1lf "
                        "%7.4lf %.2lf\n",
            d->a, d->k, d->x0, d->x1, m->I * 1E+09, d->Cbar, d->delta0);
        fprintf(stream, "\n *** Table generated with PUMAS v%.2f ***\n\n",
            PUMAS_VERSION);
        fprintf(stream,
            "      T         p     Ionization  brems     pair     "
            "photonuc  Radloss    dE/dx   CSDA Range  delta   beta\n");
        fprintf(stream, "    [MeV]    [MeV/c]  -----------------------"
                        "[MeV cm^2/g]------------------------  [g/cm^2]\n");

        /* Loop on the kinetic energy values and print the table. */
        double X = 0., dedx_last = 0.;
        int i;
        for (i = 0; i < data->n_energies; i++) {
                /* Compute the electronic energy loss. */
                double delta;
                double elec = electronic_energy_loss(
                    physics, m, data->energy[i], &delta);

                double brad[N_DEL_PROCESSES - 1];
                memset(brad, 0x0, sizeof(brad));
                struct tabulation_element * e;
                int iel;
                for (iel = 0, e = (struct tabulation_element *)data->elements;
                     iel < physics->elements_in[m->index];
                     iel++, e = (struct tabulation_element *)e->api.prev) {
                        int j;
                        for (j = 0; j < N_DEL_PROCESSES - 1; j++)
                                brad[j] += e->api.fraction *
                                    e->data[(N_DEL_PROCESSES - 1) * i + j];
                }

                /* Update the CDSA range. */
                const double radloss = brad[0] + brad[1] + brad[2];
                const double dedx = radloss + elec;
                if (i == 0)
                        X = 0.5 * data->energy[i] / dedx;
                else
                        X += 0.5 * (data->energy[i] - data->energy[i - 1]) *
                            (1. / dedx + 1. / dedx_last);
                dedx_last = dedx;

                const double p = sqrt(
                    data->energy[i] * (data->energy[i] + 2. * physics->mass));
                const double beta = p / (data->energy[i] + physics->mass);
                const double MeV = 1E+03;
                const double cmgs = 1E+04;
                fprintf(stream, "  %.3lE %.3lE %.3lE %.3lE %.3lE %.3lE "
                                "%.3lE %.3lE %.3lE %7.4lf %7.5lf\n",
                    data->energy[i] * MeV, p * MeV, elec * cmgs,
                    brad[0] * cmgs, brad[1] * cmgs, brad[2] * cmgs,
                    radloss * cmgs, dedx * cmgs, X * MeV / cmgs, delta, beta);
        }

        /* Close and return. */
        fclose(stream);
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Clear the temporary memory used for the tabulation of materials.
 */
void pumas_physics_tabulation_clear(const struct pumas_physics * physics,
    struct pumas_physics_tabulation_data * data)
{
        deallocate(data->path);
        data->path = NULL;

        struct pumas_physics_element * e;
        for (e = data->elements; e != NULL;) {
                struct pumas_physics_element * prev = e->prev;
                deallocate(e);
                e = prev;
        }
        data->elements = NULL;
}


/**
 * Get the physics differential Cross-Section (DCS) for a given process.
 */
enum pumas_return pumas_physics_dcs(
    const struct pumas_physics * physics, enum pumas_process process,
    const char ** model, pumas_dcs_t ** dcs)
{
        ERROR_INITIALISE(pumas_physics_dcs);

        if (process == PUMAS_PROCESS_BREMSSTRAHLUNG) {
                if (dcs != NULL) *dcs = physics->dcs_bremsstrahlung;
                if (model != NULL) *model = physics->model_bremsstrahlung;
        } else if (process == PUMAS_PROCESS_PAIR_PRODUCTION) {
                if (dcs != NULL) *dcs = physics->dcs_pair_production;
                if (model != NULL) *model = physics->model_pair_production;
        } else if (process == PUMAS_PROCESS_PHOTONUCLEAR) {
                if (dcs != NULL) *dcs = physics->dcs_photonuclear;
                if (model != NULL) *model = physics->model_photonuclear;
        } else {
                return ERROR_FORMAT(PUMAS_RETURN_INDEX_ERROR,
                    "bad process (expected a value in [0, 2], got %u)",
                    process);
        }

        return PUMAS_RETURN_SUCCESS;
}

/* Radiation logarithm calculated with Hartree-Fock model
 * Ref: Kelner, Kokoulin & Petrukhin (199), Physics of Atomic Nuclei, 62(11),
 *      1894-1898. doi:101134/1855464
 *
 * Values have been taken from Koehne et al.
 * (https://doi.org/10.1016/j.cpc.2013.04.001) since the original paper does not
 * seem to be available online.
 */
static double radiation_logarithm(double Z)
{
        const int i = (int)Z;
        if ((i >= 1) && (i <= 92)) {
                static double Lz[92] = {
                    202.4, 151.9, 159.9, 172.3, 177.9,
                    178.3, 176.6, 173.4, 170.0, 165.8,
                    165.8, 167.1, 169.1, 170.8, 172.2,
                    173.4, 174.3, 174.8, 175.1, 175.6,
                    176.2, 176.8,   0.0,   0.0,   0.0,
                    175.8,   0.0,   0.0, 173.1,   0.0,
                      0.0, 173.0,   0.0,   0.0, 173.5,
                      0.0,   0.0,   0.0,   0.0,   0.0,
                      0.0, 175.9,   0.0,   0.0,   0.0,
                      0.0,   0.0,   0.0,   0.0, 177.4,
                      0.0,   0.0, 178.6,   0.0,   0.0,
                      0.0,   0.0,   0.0,   0.0,   0.0,
                      0.0,   0.0,   0.0,   0.0,   0.0,
                      0.0,   0.0,   0.0,   0.0,   0.0,
                      0.0,   0.0,   0.0, 177.6,   0.0,
                      0.0,   0.0,   0.0,   0.0,   0.0,
                      0.0, 178.0,   0.0,   0.0,   0.0,
                      0.0,   0.0,   0.0,   0.0,   0.0,
                      0.0, 179.8
                };
                const double l = Lz[i - 1];
                if (l > 0) return l;
        }
        return 182.7;
}

/**
 * The Bremsstrahlung differential cross section according to
 * Kelner, Kokoulin & Petrukhin (KKP).
 *
 * @param Z       The charge number of the target atom.
 * @param A       The mass number of the target atom.
 * @param mu      The projectile rest mass, in GeV
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The corresponding value of the atomic DCS, in m^2 / GeV.
 *
 * The KKP Bremsstrahlung differential cross-section was initially implemented
 * following Groom et al. (see e.g.
 * http://pdg.lbl.gov/2020/AtomicNuclearProperties/adndt.pdf). Then, it has been
 * refined following Koehne et al. (https://doi.org/10.1016/j.cpc.2013.04.001)
 * by taking into account the nucleus excitation term (See e.g. Kelner et al.
 * https://cds.cern.ch/record/288828/files/MEPHI-024-95.pdf) and more accurate
 * radiation logarithm computations.
 */
static double dcs_bremsstrahlung_KKP(
    double Z, double A, double mu, double K, double q)
{
        if ((Z <= 0) || (A <= 0) || (mu <= 0) || (K <= 0) || (q <= 0) ||
            (q >= K + mu))
                return 0.;

        const double me = ELECTRON_MASS;
        const double sqrte = 1.648721271;
        const double phie_factor = mu / (me * me * sqrte);
        const double rem = 5.63588E-13 * me / mu;

        const double BZ_n = radiation_logarithm(Z) * pow(Z, -1. / 3.);
        const double BZ_e = (Z == 1.) ? 446. : 1429. * pow(Z, -2. / 3.);
        const double D_n = 1.54 * pow(A, 0.27);
        const double E = K + mu;
        const double dcs_factor = 7.297182E-07 * rem * rem * Z / E;

        const double delta_factor = 0.5 * mu * mu / E;
        const double qe_max = E / (1. + 0.5 * mu * mu / (me * E));

        const double nu = q / E;
        const double delta = delta_factor * nu / (1. - nu);
        const double muD_factor = mu + delta * (D_n * sqrte - 2.);
        double Phi_n, Phi_x, Phi_e;
        Phi_n = log(BZ_n * muD_factor /
            (D_n * (me + delta * sqrte * BZ_n)));
        if (Phi_n < 0.) Phi_n = 0.;
        if (Z >= 2) {
                Phi_x = log(mu * D_n / muD_factor);
                if (Phi_x < 0.) Phi_x = 0.;
        } else {
                Phi_x = 0.;
        }
        if (q < qe_max) {
                Phi_e = log(BZ_e * mu /
                    ((1. + delta * phie_factor) * (me + delta * sqrte * BZ_e)));
                if (Phi_e < 0.) Phi_e = 0.;
        } else
                Phi_e = 0.;

        const double dcs = dcs_factor *
            (Z * Phi_n + Phi_x + Phi_e) * (4. / 3. * (1. / nu - 1.) + nu);
        return (dcs < 0.) ? 0. : dcs;
}

/* Bremsstahlung DCS according to Andreev, Bezrukov & Bugaev.
 *
 * @param Z       The charge number of the target atom.
 * @param A       The mass number of the target atom.
 * @param mu      The projectile rest mass, in GeV
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The corresponding value of the atomic DCS, in m^2 / GeV.
 *
 * Ref: https://arxiv.org/abs/hep-ph/0010322 (MUM)
 *
 * PROPOSAL implementation converted to C
 * Ref: https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/private/PROPOSAL/crossection/parametrization/Bremsstrahlung.cxx
 */
static double dcs_bremsstrahlung_ABB(
    double Z, double A, double m, double K, double q)
{
#define ALPHA 0.0072973525664
#define SQRTE 1.648721270700128
#define ME    0.5109989461
#define RE    2.8179403227E-13
#define MMU   105.6583745

        /* Check inputs */
        if ((Z <= 0) || (A <= 0) || (m <= 0) || (K <= 0) || (q <= 0) ||
            (q >= K + m))
                return 0.;

        /* Convert from GeV to MeV */
        const double energy = (K + m) * 1E+03;
        const double v = q * 1E+03 / energy;
        m *= 1E+03;

        /* Least momentum transferred to the nucleus (eq. 2.2) */
        const double Z3 = pow(Z, -1. / 3);
        const double a1 = 184.15 * Z3 / (SQRTE * ME);    /* eq 2.18 */
        const double a2 = 1194 * Z3 * Z3 / (SQRTE * ME); /* eq.2.19 */

        /* Calculating the contribution of elastic nuclear and atomic form
         * factors (eq. 2.30)
         */
        const double qc   = 1.9 * MMU * Z3;
        double aux        = 2 * m / qc;
        const double zeta = sqrt(1 + aux * aux);

        const double delta = m * m * v / (2 * energy * (1 - v));
        const double x1    = a1 * delta;
        const double x2    = a2 * delta;

        double aux1, aux2, d1, d2, psi1, psi2;

        if (Z == 1) {
                d1 = 0;
                d2 = 0;
        } else {
                aux1 = log(m / qc);
                aux2 = 0.5 * zeta * log((zeta + 1) / (zeta - 1));
                d1   = aux1 + aux2;
                d2   = aux1 + 0.5 * ((3 - zeta * zeta) * aux2 + aux * aux);
        }

        /* eq. 2.20 and 2.21 */
        aux  = m * a1;
        aux1 = log(aux * aux / (1 + x1 * x1));
        aux  = m * a2;
        aux2 = log(aux * aux / (1 + x2 * x2));
        psi1 = 0.5 * ((1 + aux1) + (1 + aux2) / Z);
        psi2 = 0.5 * ((2. / 3 + aux1) + (2. / 3 + aux2) / Z);

        aux1 = x1 * atan(1 / x1);
        aux2 = x2 * atan(1 / x2);
        psi1 -= aux1 + aux2 / Z;
        aux = x1 * x1;
        psi2 += 2 * aux * (1 - aux1 + 0.75 * log(aux / (1 + aux)));
        aux = x2 * x2;
        psi2 += 2 * aux * (1 - aux2 + 0.75 * log(aux / (1 + aux))) / Z;

        psi1 -= d1;
        psi2 -= d2;
        const double result = (2 - 2 * v + v * v) * psi1 -
            (2. / 3) * (1 - v) * psi2;

        if (result < 0) return 0;

        aux = 2 * (ME / m) * RE * Z;
        return aux * aux * (ALPHA / q) * result * 1E-04;
}

/* Bremsstahlung DCS according to Sandrock, Soedingrekso & Rhode.
 *
 * @param Z       The charge number of the target atom.
 * @param A       The mass number of the target atom.
 * @param mu      The projectile rest mass, in GeV
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The corresponding value of the atomic DCS, in m^2 / GeV.
 *
 * Ref: https://arxiv.org/abs/1910.07050
 *
 * PROPOSAL implementation converted to C
 * Ref: https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/private/PROPOSAL/crossection/parametrization/Bremsstrahlung.cxx
 */
static double dcs_bremsstrahlung_SSR(
    double Z, double A, double m, double K, double q)
{
        /* Check inputs */
        if ((Z <= 0) || (A <= 0) || (m <= 0) || (K <= 0) || (q <= 0) ||
            (q >= K + m))
                return 0.;

        const double a[3] = {
            -0.00349, 148.84, -987.531};
        const double b[4] = {
            0.1642, 132.573, -585.361, 1407.77};
        const double c[6] = {
            -2.8922, -19.0156, 57.698, -63.418, 14.1166, 1.84206};
        const double d[6] = {
            2134.19, 581.823, -2708.85, 4767.05, 1.52918, 0.361933};

        /* Convert from GeV to MeV */
        const double energy = (K + m) * 1E+03;
        const double v = q * 1E+03 / energy;
        m *= 1E+03;

        const double Z13 = pow(Z, -1. / 3);
        const double rad_log = radiation_logarithm(Z);
        const double rad_log_inel = (Z == 1.) ? 446. : 1429.;
        const double Dn = 1.54 * pow(A, 0.27);

        const double mu_qc = m / (MMU * exp(1.) / Dn);
        const double rho = sqrt(1. + 4. * mu_qc * mu_qc);

        const double log_rho = log((rho + 1.) / (rho - 1.));
        const double delta1 = log(mu_qc) + 0.5 * rho * log_rho;
        const double delta2 = log(mu_qc) +
            0.25 * (3. * rho - rho * rho * rho) * log_rho + 2. * mu_qc * mu_qc;

        /* Least momentum transferred to the nucleus (eq. 7) */
        const double delta = m * m * v / (2. * energy * (1. - v));

        double phi1 = log(rad_log * Z13 * (m / ME) /
            (1. + rad_log * Z13 * exp(0.5) * delta / ME));
        double phi2 = log(rad_log * Z13 * exp(-1. / 6.) * (m / ME) /
            (1. + rad_log * Z13 * exp(1. / 3.) * delta / ME));
        phi1 -= delta1 * (1. - 1. / Z);
        phi2 -= delta2 * (1. - 1. / Z);

        /* s_atomic */
        const double square_momentum  = (energy - m) * (energy + m);
        const double particle_momentum = (square_momentum > 0) ?
            sqrt(square_momentum) : 0.;
        const double maxV = ME * (energy - m) /
            (energy * (energy - particle_momentum + ME));

        double s_atomic = 0.0;

        if (v < maxV) {
                const double s_atomic_1 = log(m / delta /
                    (m * delta / (ME * ME) + SQRTE));
                const double s_atomic_2 = log(1. + ME /
                    (delta * rad_log_inel * Z13 * Z13 * SQRTE));
                s_atomic = (4. / 3. * (1. - v) + v * v) *
                    (s_atomic_1 - s_atomic_2);
        }

        /* s_rad */
        double s_rad;

        if (v < .0 || v > 1.) {
                s_rad = 0.;
        } else if (v < 0.02) {
                s_rad = a[0] + a[1] * v + a[2] * v * v;
        } else if (v >= 0.02 && v < 0.1) {
                s_rad = b[0] + b[1] * v + b[2] * v * v + b[3] * v * v * v;
        } else if (v >= 0.01 && v < 0.9) {
                s_rad = c[0] + c[1] * v + c[2] * v * v;

                const double tmp = log(1. - v);
                s_rad += c[3] * v * log(v) + c[4] * tmp + c[5] * tmp * tmp;
        } else {
                s_rad = d[0] + d[1] * v + d[2] * v * v;

                const double tmp = log(1. - v);
                s_rad += d[3] * v * log(v) + d[4] * tmp + d[5] * tmp * tmp;
        }

        const double result = ((2. - 2. * v + v * v) * phi1 - 2. / 3. *
            (1. - v) * phi2) + 1. / Z * s_atomic + 0.25 * ALPHA * phi1 * s_rad;

        if (result <= 0.) return 0.;

        const double aux = 2 * (ME / m) * RE * Z;
        return aux * aux * (ALPHA / q) * result * 1E-04;

#undef ALPHA
#undef SQRTE
#undef ME
#undef RE
#undef MMU
}

/**
 * The e+e- pair production differential cross section according to Kelner,
 * Kokoulin & Petrukhin.
 *
 * @param Z       The charge number of the target atom.
 * @param A       The mass number of the target atom.
 * @param mu      The projectile rest mass, in GeV
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The corresponding value of the atomic DCS, in m^2 / GeV.
 *
 * Geant4 like implementation, see e.g. section 11.3.1 of the Geant4 Physics
 * Reference Manual.
 */
static double dcs_pair_production_KKP(
    double Z, double A_, double mass, double K, double q)
{
        if ((Z <= 0) || (A_ <= 0) || (mass <= 0) || (K <= 0) || (q <= 0))
                return 0.;

/*
 * Coefficients for the Gaussian quadrature from:
 * https://pomax.github.io/bezierinfo/legendre-gauss.html.
 */
#define N_GQ 8
        const double xGQ[N_GQ] = { 0.01985507, 0.10166676, 0.2372338,
                0.40828268, 0.59171732, 0.7627662, 0.89833324, 0.98014493 };
        const double wGQ[N_GQ] = { 0.05061427, 0.11119052, 0.15685332,
                0.18134189, 0.18134189, 0.15685332, 0.11119052, 0.05061427 };

        /*  Check the bounds of the energy transfer. */
        if (q <= 4. * ELECTRON_MASS) return 0.;
        const double sqrte = 1.6487212707;
        const double Z13 = pow(Z, 1. / 3.);
        if (q >= K + mass * (1. - 0.75 * sqrte * Z13)) return 0.;

        /*  Precompute some constant factors for the integration. */
        const double nu = q / (K + mass);
        const double r = mass / ELECTRON_MASS;
        const double beta = 0.5 * nu * nu / (1. - nu);
        const double xi_factor = 0.5 * r * r * beta;
        const double A = (Z == 1.) ? 202.4 : 183.;
        const double AZ13 = A / Z13;
        const double cL = 2. * sqrte * ELECTRON_MASS * AZ13;
        const double cLe = 2.25 * Z13 * Z13 / (r * r);

        /*  Compute the bound for the integral. */
        const double gamma = 1. + K / mass;
        const double x0 = 4. * ELECTRON_MASS / q;
        const double x1 = 6. / (gamma * (gamma - q / mass));
        const double argmin =
            (x0 + 2. * (1. - x0) * x1) / (1. + (1. - x1) * sqrt(1. - x0));
        if ((argmin >= 1.) || (argmin <= 0.)) return 0.;
        const double tmin = log(argmin);

        /*  Compute the integral over t = ln(1-rho). */
        double I = 0.;
        int i;
        for (i = 0; i < N_GQ; i++) {
                const double eps = exp(xGQ[i] * tmin);
                const double rho = 1. - eps;
                const double rho2 = rho * rho;
                const double rho21 = eps * (2. - eps);
                const double xi = xi_factor * rho21;
                const double xi_i = 1. / xi;

                /* Compute the e-term. */
                double Be;
                if (xi >= 1E+03)
                        Be =
                            0.5 * xi_i * ((3 - rho2) + 2. * beta * (1. + rho2));
                else
                        Be = ((2. + rho2) * (1. + beta) + xi * (3. + rho2)) *
                                log(1. + xi_i) +
                            (rho21 - beta) / (1. + xi) - 3. - rho2;
                const double Ye = (5. - rho2 + 4. * beta * (1. + rho2)) /
                    (2. * (1. + 3. * beta) * log(3. + xi_i) - rho2 -
                                      2. * beta * (2. - rho2));
                const double xe = (1. + xi) * (1. + Ye);
                const double cLi = cL / rho21;
                const double Le = log(AZ13 * sqrt(xe) * q / (q + cLi * xe)) -
                    0.5 * log(1. + cLe * xe);
                double Phi_e = Be * Le;
                if (Phi_e < 0.) Phi_e = 0.;

                /* Compute the mu-term. */
                double Bmu;
                if (xi <= 1E-03)
                        Bmu = 0.5 * xi * (5. - rho2 + beta * (3. + rho2));
                else
                        Bmu = ((1. + rho2) * (1. + 1.5 * beta) -
                                  xi_i * (1. + 2. * beta) * rho21) *
                                log(1. + xi) +
                            xi * (rho21 - beta) / (1. + xi) +
                            (1. + 2. * beta) * rho21;
                const double Ymu = (4. + rho2 + 3. * beta * (1. + rho2)) /
                    ((1. + rho2) * (1.5 + 2. * beta) * log(3. + xi) + 1. -
                                       1.5 * rho2);
                const double xmu = (1. + xi) * (1. + Ymu);
                const double Lmu =
                    log(r * AZ13 * q / (1.5 * Z13 * (q + cLi * xmu)));
                double Phi_mu = Bmu * Lmu;
                if (Phi_mu < 0.) Phi_mu = 0.;

                /* Update the t-integral. */
                I -= (Phi_e + Phi_mu / (r * r)) * (1. - rho) * wGQ[i] * tmin;
        }

        /* Atomic electrons form factor. */
        double zeta;
        if (gamma <= 35.)
                zeta = 0.;
        else {
                double gamma1, gamma2;
                if (Z == 1.) {
                        gamma1 = 4.4E-05;
                        gamma2 = 4.8E-05;
                } else {
                        gamma1 = 1.95E-05;
                        gamma2 = 5.30E-05;
                }
                zeta = 0.073 * log(gamma / (1. + gamma1 * gamma * Z13 * Z13)) -
                    0.26;
                if (zeta <= 0.)
                        zeta = 0.;
                else {
                        zeta /=
                            0.058 * log(gamma / (1. + gamma2 * gamma * Z13)) -
                            0.14;
                }
        }

        /* Gather the results and return the macroscopic DCS. */
        const double E = K + mass;
        const double dcs = 1.794664E-34 * Z * (Z + zeta) * (E - q) * I /
            (q * E);
        return (dcs < 0.) ? 0. : dcs;

#undef N_GQ
}

/**
 * The e+e- pair production doubly differential cross section according to
 * Sandrock, Soedingrekso & Rhode.
 *
 * @param Z       The charge number of the target atom.
 * @param A       The mass number of the target atom.
 * @param mu      The projectile rest mass, in GeV
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @param rho     The e+e- asymmetry.
 * @return The corresponding value of the atomic DCS, in m^2 / GeV.
 *
 * Ref: https://arxiv.org/abs/1910.07050
 *
 * PROPOSAL implementation converted to C
 * Ref: https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/private/PROPOSAL/crossection/parametrization/EPairProduction.cxx
 */
static inline double dcs_pair_production_d2_SSR(
    double Z, double A, double m, double K, double q, double rho)
{
#define ALPHA 0.0072973525664
#define ME    0.5109989461
#define RE    2.8179403227E-13

        const double energy = (K + m) * 1E+03;
        const double v = q / (K + m);
        m *= 1E+03;
        const double rad_log = radiation_logarithm(Z);

        const double const_prefactor = 4. / (3. * M_PI) * Z *
            pow(ALPHA * RE, 2.);
        const double Z13 = pow(Z, -1. / 3.);
        const double d_n = 1.54 * pow(A, 0.27);

        rho = 1 - rho;
        const double rho2 = rho * rho;

        /* Zeta */
        double g1, g2;
        if (Z == 1.) {
                g1 = 4.4E-05;
                g2 = 4.8E-05;
        } else {
                g1 = 1.95E-05;
                g2 = 5.3E-05;
        }

        const double zeta1 = (0.073 * log(energy / m /
            (1. + g1 * pow(Z, 2. / 3.) * energy / m)) - 0.26);
        const double zeta2 = (0.058 * log(energy / m /
            (1 + g2 / Z13 * energy / m)) - 0.14);

        double zeta;
        if ((zeta1 > 0.) && (zeta2 > 0.)) {
                zeta = zeta1 / zeta2;
        } else {
                zeta = 0.;
        }

        const double beta = v * v / (2. * (1. - v));
        const double xi = pow(m * v / (2. * ME), 2.) * (1. - rho2) / (1. - v);

        /* Diagram e */
        const double Be = ((2. + rho2) * (1. + beta) +
            xi * (3. + rho2)) * log(1. + 1. / xi) +
            (1. - rho2 - beta) / (1. + xi) - (3. + rho2);

        const double Ce2 = ((1. - rho2) * (1. + beta) +
            xi * (3. - rho2)) * log(1. + 1. / xi) +
            2. * (1. - beta - rho2) / (1. + xi) - (3. - rho2);
        const double Ce1 = Be - Ce2;

        const double De = ((2. + rho2) * (1. + beta) +
            xi * (3. + rho2)) * math_dilog(1. / (1. + xi)) -
            (2. + rho2) * xi * log(1. + 1. / xi) -
            (xi + rho2 + beta) / (1. + xi);

        double Le1, Le2;
        if (De / Be > 0.) {
                const double Xe = exp(-De / Be);
                Le1 = log(rad_log * Z13 * sqrt(1. + xi) /
                    (Xe + 2. * ME * exp(0.5) * rad_log * Z13 * (1. + xi) /
                    (energy * v * (1. - rho2)))) - De / Be -
                    0.5 * log(Xe + pow(ME / m * d_n, 2.) * (1. + xi));

                Le2 = log(rad_log * Z13 * exp(-1. / 6.) * sqrt(1 + xi) /
                    (Xe + 2. * ME * exp(1. / 3.) * rad_log * Z13 * (1. + xi) /
                    (energy * v * (1. - rho2)))) - De / Be -
                    0.5 * log(Xe + pow(ME / m * d_n, 2.) *
                    exp(-1. / 3.) * (1. + xi));
        } else {
                const double Xe_inv = exp(De / Be);
                Le1 = log(rad_log * Z13 * sqrt(1. + xi) /
                    (1. + Xe_inv * 2. * ME * exp(0.5) * rad_log * Z13 *
                    (1. + xi) / (energy * v * (1. - rho2)))) - 0.5 * De / Be -
                    0.5 * log(1. + Xe_inv * pow(ME / m * d_n, 2.) * (1. + xi));

                Le2 = log(rad_log * Z13 * exp(-1. / 6.) * sqrt(1 + xi) /
                    (1. + Xe_inv * 2. * ME * exp(1. / 3.) * rad_log * Z13 *
                    (1. + xi) / (energy * v * (1. - rho2)))) - 0.5 * De / Be -
                    0.5 * log(1. + Xe_inv * pow(ME / m * d_n, 2.) *
                    exp(-1. / 3.) * (1. + xi));
        }

        double diagram_e = const_prefactor * (Z + zeta) * (1. - v) / v *
            (Ce1 * Le1 + Ce2 * Le2);
        if (diagram_e < 0.) diagram_e = 0.;

        /* Diagram mu */
        const double Bm = ((1. + rho2) * (1. + (3. * beta) / 2) - 1. / xi *
            (1. + 2. * beta) * (1. - rho2)) * log(1. + xi) +
            xi * (1. - rho2 - beta) / (1. + xi) +
            (1. + 2. * beta) * (1. - rho2);

        const double Cm2 = ((1. - beta) * (1. - rho2) -
            xi * (1. + rho2)) * log(1. + xi) / xi -
            2. * (1. - beta - rho2) / (1. + xi) +
            1. - beta - (1. + beta) * rho2;
        const double Cm1 = Bm - Cm2;

        const double Dm = ((1. + rho2) * (1. + (3. * beta) / 2.) -
            1. / xi * (1. + 2. * beta) * (1. - rho2)) *
            math_dilog(xi / (1. + xi)) + (1. + (3. * beta) / 2.) *
            (1. - rho2) / xi * log(1. + xi) + (1. - rho2 - beta / 2. *
            (1. + rho2) + (1. - rho2) / (2. * xi) * beta) * xi / (1. + xi);

        double Lm1, Lm2;
        if (Dm / Bm > 0.) {
                const  double Xm = exp(-Dm / Bm);
                Lm1 = log(Xm * m / ME * rad_log * Z13 / d_n /
                    (Xm + 2. * ME * exp(0.5) * rad_log * Z13 * (1. + xi) /
                    (energy * v * (1. - rho2))));
                Lm2 = log(Xm * m / ME * rad_log * Z13 / d_n /
                    (Xm + 2. * ME * exp(1. / 3.) * rad_log * Z13 * (1. + xi) /
                    (energy * v * (1. - rho2))));
        } else {
                const double Xmv = exp(Dm / Bm);
                Lm1 = log(m / ME * rad_log * Z13 / d_n /
                    (1. + 2. * ME * exp(0.5) * rad_log * Z13 * (1. + xi) /
                    (energy * v * (1. - rho2)) * Xmv));
                Lm2 = log(m / ME * rad_log * Z13 / d_n /
                    (1. + 2. * ME * exp(1. / 3.) * rad_log * Z13 * (1. + xi) /
                    (energy * v * (1. - rho2)) * Xmv));
        }

        double diagram_mu = const_prefactor * (Z + zeta) * (1. - v) / v *
            pow(ME / m, 2.) * (Cm1 * Lm1 + Cm2 * Lm2);
        if (diagram_mu < 0.) diagram_mu = 0.;

        return (diagram_e + diagram_mu) * 1E-01 / energy;

#undef ALPHA
#undef ME
#undef RE
}

/**
 * The e+e- pair production differential cross section according to Sandrock,
 * Soedingrekso & Rhode.
 *
 * @param Z       The charge number of the target atom.
 * @param A       The mass number of the target atom.
 * @param mu      The projectile rest mass, in GeV
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The corresponding value of the atomic DCS, in m^2 / GeV.
 *
 * Mixed implementation. The DDCS of PROPOSAL is used but the numeric
 * integration is done with a Gaussian quadrature a la Geant4.
 */
static double dcs_pair_production_SSR(
    double Z, double A_, double mass, double K, double q)
{
        if ((Z <= 0) || (A_ <= 0) || (mass <= 0) || (K <= 0) || (q <= 0))
                return 0.;
/*
 * Coefficients for the Gaussian quadrature from:
 * https://pomax.github.io/bezierinfo/legendre-gauss.html.
 */
#define N_GQ 8
        const double xGQ[N_GQ] = { 0.01985507, 0.10166676, 0.2372338,
                0.40828268, 0.59171732, 0.7627662, 0.89833324, 0.98014493 };
        const double wGQ[N_GQ] = { 0.05061427, 0.11119052, 0.15685332,
                0.18134189, 0.18134189, 0.15685332, 0.11119052, 0.05061427 };

        /*  Check the bounds of the energy transfer. */
        if (q <= 4. * ELECTRON_MASS) return 0.;
        const double sqrte = 1.6487212707;
        const double Z13 = pow(Z, 1. / 3.);
        if (q >= K + mass * (1. - 0.75 * sqrte * Z13)) return 0.;

        /* Compute the bound for the integral */
        const double gamma = 1. + K / mass;
        const double x0 = 1. - 4. * ELECTRON_MASS / q;
        const double x1 = 1. - 6. / (gamma * (gamma - q / mass));

        double rmax;
        if ((x0 > 0) && (x1 > 0)) {
                rmax = sqrt(x0) * x1;
        } else {
                return 0.;
        }

        const double ri = 1. - rmax;
        if ((ri <= 0.) || (ri >= 1.)) return 0.;
        const double tmin = log(ri);

        /*  Compute the integral over t = ln(rho) */
        double I = 0.;
        int i;
        for (i = 0; i < N_GQ; i++) {
                const double rho = exp(xGQ[i] * tmin);
                I -= dcs_pair_production_d2_SSR(Z, A_, mass, K, q, rho) *
                    rho * wGQ[i] * tmin;
        }

        return (I < 0) ? 0. : I;

#undef N_GQ
}

/** ALLM97 parameterisation of the proton structure function, F2.
 *
 * @param x       The Bjorken x parameter.
 * @param Q2      The negative four momentum squared.
 * @return The corresponding value of the proton structure function, F2.
 *
 * References:
 *      DESY 97-251 [arXiv:hep-ph/9712415].
 */
static inline double dcs_photonuclear_f2p_ALLM97(double x, double Q2)
{
        const double m02 = 0.31985;
        const double mP2 = 49.457;
        const double mR2 = 0.15052;
        const double Q02 = 0.52544;
        const double Lambda2 = 0.06527;

        const double cP1 = 0.28067;
        const double cP2 = 0.22291;
        const double cP3 = 2.1979;
        const double aP1 = -0.0808;
        const double aP2 = -0.44812;
        const double aP3 = 1.1709;
        const double bP1 = 0.36292;
        const double bP2 = 1.8917;
        const double bP3 = 1.8439;

        const double cR1 = 0.80107;
        const double cR2 = 0.97307;
        const double cR3 = 3.4942;
        const double aR1 = 0.58400;
        const double aR2 = 0.37888;
        const double aR3 = 2.6063;
        const double bR1 = 0.01147;
        const double bR2 = 3.7582;
        const double bR3 = 0.49338;

        const double M = 0.5 * (PROTON_MASS + NEUTRON_MASS);
        const double M2 = M * M;
        const double W2 = M2 + Q2 * (1.0 / x - 1.0);
        const double t = log(log((Q2 + Q02) / Lambda2) / log(Q02 / Lambda2));
        const double xP = (Q2 + mP2) / (Q2 + mP2 + W2 - M2);
        const double xR = (Q2 + mR2) / (Q2 + mR2 + W2 - M2);
        const double lnt = log(t);
        const double cP =
            cP1 + (cP1 - cP2) * (1.0 / (1.0 + exp(cP3 * lnt)) - 1.0);
        const double aP =
            aP1 + (aP1 - aP2) * (1.0 / (1.0 + exp(aP3 * lnt)) - 1.0);
        const double bP = bP1 + bP2 * exp(bP3 * lnt);
        const double cR = cR1 + cR2 * exp(cR3 * lnt);
        const double aR = aR1 + aR2 * exp(aR3 * lnt);
        const double bR = bR1 + bR2 * exp(bR3 * lnt);

        const double F2P = cP * exp(aP * log(xP) + bP * log(1 - x));
        const double F2R = cR * exp(aR * log(xR) + bR * log(1 - x));

        return Q2 / (Q2 + m02) * (F2P + F2R);
}

/* The F2 nuclear structure function following Dutta et al.
 *
 * @param A            The atomic charge number.
 * @param A            The atomic weight.
 * @param F2p          The proton structure function, F2.
 * @param shadowing    The shadowing factor.
 * @param x            The Bjorken x parameter.
 * @return The corresponding value of the nuclear structure function, F2a.
 *
 * The F2a structure function for a nucleus of charge number Z and atomic
 * weight A is computed according to Dutta et al.
 *
 * References:
 *      Dutta et al., Phys.Rev. D63 (2001) 094020 [arXiv:hep-ph/0012350].
 */
static inline double dcs_photonuclear_f2a_DRSS(
    double Z, double A, double F2p, double shadowing, double x)
{
        return F2p * shadowing * (Z + (A - Z) *
            (1.0 + x * (-1.85 + x * (2.45 + x * (-2.35 + x)))));
}

/* Shadowing factor following Dutta et al.
 *
 * @param A       The atomic charge number.
 * @param A       The atomic weight.
 * @param x       The Bjorken x parameter.
 * @return The corresponding value of the shadowing factor.
 *
 * The shadowing factor for a nucleus of atomic weight A is computed according
 * to Dutta et al.
 *
 * References:
 *      Dutta et al., Phys.Rev. D63 (2001) 094020 [arXiv:hep-ph/0012350].
 */
static inline double dcs_photonuclear_shadowing_DRSS(
    double Z, double A, double x)
{
        if (Z == 1.) return 1.;

        if (x < 0.0014)
                return exp(-0.1 * log(A));
        else if (x < 0.04)
                return exp((0.069 * log10(x) + 0.097) * log(A));
        else
                return 1.;
}

/** The doubly differential cross sections d^2S/(dq*dQ2) for photonuclear
 * interactions following Dutta et al.
 *
 * @param Z       The target charge number.
 * @param A       The target atomic weight.
 * @param ml      The projectile mass.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @param Q2      The negative four momentum squared.
 * @return The doubly differential cross section in m^2/kg/GeV^3.
 *
 * References:
 *      Dutta et al., Phys.Rev. D63 (2001) 094020 [arXiv:hep-ph/0012350].
 */
static inline double dcs_photonuclear_d2_DRSS(
    double Z, double A, double ml, double K, double q, double Q2)
{
        const double cf = 2.6056342605319227E-35;
        const double M = 0.5 * (PROTON_MASS + NEUTRON_MASS);
        const double E = K + ml;

        const double y = q / E;
        const double x = 0.5 * Q2 / (M * q);
        const double F2p = dcs_photonuclear_f2p_ALLM97(x, Q2);
        const double shadowing = dcs_photonuclear_shadowing_DRSS(Z, A, x);
        const double F2A = dcs_photonuclear_f2a_DRSS(Z, A, F2p, shadowing, x);
        const double R = 0.;

        const double dds = (1 - y +
                               0.5 * (1 - 2 * ml * ml / Q2) *
                                   (y * y + Q2 / (E * E)) / (1 + R)) /
                (Q2 * Q2) -
            0.25 / (E * E * Q2);

        return cf * F2A * dds / q;
}

/**
 * The photonuclear differential cross section following Dutta et al.
 *
 * @param Z       The charge number of the target atom.
 * @param A       The mass number of the target atom.
 * @param mu      The projectile rest mass, in GeV
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The corresponding value of the atomic DCS, in m^2 / GeV.
 *
 * The photonuclear differential cross-section is computed following DRSS,
 * with ALLM97 parameterisation of the structure function F2.
 *
 * References:
 *      Dutta et al., Phys.Rev. D63 (2001) 094020 [arXiv:hep-ph/0012350].
 */
static double dcs_photonuclear_DRSS(
    double Z, double A, double ml, double K, double q)
{
        if ((Z <= 0) || (A <= 0) || (ml <= 0) || (K <= 0) || (q <= 0))
                return 0.;
/*
 * Coefficients for the Gaussian quadrature from:
 * https://pomax.github.io/bezierinfo/legendre-gauss.html.
 */
#define N_GQ 9
        const double xGQ[N_GQ] = { 0.0000000000000000, -0.8360311073266358,
                0.8360311073266358, -0.9681602395076261, 0.9681602395076261,
                -0.3242534234038089, 0.3242534234038089, -0.6133714327005904,
                0.6133714327005904 };
        const double wGQ[N_GQ] = { 0.3302393550012598, 0.1806481606948574,
                0.1806481606948574, 0.0812743883615744, 0.0812743883615744,
                0.3123470770400029, 0.3123470770400029, 0.2606106964029354,
                0.2606106964029354 };

        const double M = 0.931494;
        const double mpi = 0.134977;
        const double E = K + ml;

        double ds = 0.;
        if ((q >= (E - ml)) || (q <= (mpi * (1.0 + 0.5 * mpi / M)))) return ds;

        const double y = q / E;
        const double Q2min = ml * ml * y * y / (1 - y);
        const double Q2max = 2.0 * M * (q - mpi) - mpi * mpi;
        if ((Q2max < Q2min) | (Q2min < 0)) return ds;

        /* Set the binning. */
        const double pQ2min = log(Q2min);
        const double pQ2max = log(Q2max);
        const double dpQ2 = pQ2max - pQ2min;
        const double pQ2c = 0.5 * (pQ2max + pQ2min);

        /*
         * Integrate the doubly differential cross-section over Q2 using
         * a Gaussian quadrature. Note that 9 points are enough to get a
         * better than 0.1 % accuracy.
         */
        int i;
        for (i = 0; i < N_GQ; i++) {
                const double Q2 = exp(pQ2c + 0.5 * dpQ2 * xGQ[i]);
                ds += dcs_photonuclear_d2_DRSS(Z, A, ml, K, q, Q2) *
                    Q2 * wGQ[i];
        }

        if (ds < 0.) ds = 0.;
        return 0.5 * ds * dpQ2;

#undef N_GQ
}

/** Data structure for caracterising a DCS model */
struct dcs_entry {
        enum pumas_process process;
        const char * model;
        pumas_dcs_t * dcs;
};

/** Stack (library) of available DCS models
 *
 * Note that the first entry for a given process is taken as the default
 * model for the corresponding process.
 */
#define DCS_STACK_SIZE 64
static struct dcs_entry dcs_stack[DCS_STACK_SIZE] = {
    {PUMAS_PROCESS_BREMSSTRAHLUNG,  "KKP",  &dcs_bremsstrahlung_KKP},
    {PUMAS_PROCESS_BREMSSTRAHLUNG,  "ABB",  &dcs_bremsstrahlung_ABB},
    {PUMAS_PROCESS_BREMSSTRAHLUNG,  "SSR",  &dcs_bremsstrahlung_SSR},
    {PUMAS_PROCESS_PAIR_PRODUCTION, "KKP",  &dcs_pair_production_KKP},
    {PUMAS_PROCESS_PAIR_PRODUCTION, "SSR",  &dcs_pair_production_SSR},
    {PUMAS_PROCESS_PHOTONUCLEAR,    "DRSS", &dcs_photonuclear_DRSS}
};

/** Mapping between enum and names for processes */
static const char * process_name[3] = {
        "bremsstrahlung", "pair production", "photonuclear"};

/** Routine for checking if a model's name exists */
static enum pumas_return dcs_check_model(enum pumas_process process,
     const char * model, struct error_context * error_)
{
        int i;
        struct dcs_entry * entry;
        for (i = 0, entry = dcs_stack;
            (i < DCS_STACK_SIZE) && (entry->model != NULL); i++, entry++) {
                if ((entry->process == process) &&
                    (strcmp(entry->model, model) == 0)) {
                        return PUMAS_RETURN_SUCCESS;
                }
        }

        return ERROR_VREGISTER(PUMAS_RETURN_MODEL_ERROR,
            "cannot find %s model for %s process", model,
            process_name[process]);
}

/** Routine for checking a process index */
static enum pumas_return dcs_check_process(
    enum pumas_process process, struct error_context * error_)
{
        if ((process < 0) || (process >= 3)) {
                return ERROR_VREGISTER(PUMAS_RETURN_INDEX_ERROR,
                    "bad process (expected an index in [0, 3], got %u)",
                    process);
        } else {
                return PUMAS_RETURN_SUCCESS;
        }
}

/* API function for registering a DCS model */
enum pumas_return pumas_dcs_register(
    enum pumas_process process, const char * model, pumas_dcs_t * dcs)
{
        ERROR_INITIALISE(pumas_dcs_register);

        /* Check the process index */
        if (dcs_check_process(process, error_) != PUMAS_RETURN_SUCCESS) {
                return ERROR_RAISE();
        }

        /* Check that a DCS function was actually provided */
        if (dcs == NULL) {
                return ERROR_MESSAGE(PUMAS_RETURN_VALUE_ERROR,
                    "bad dcs (expected a function, got nil)");
        }

        /* Check if the model is already registered */
        if (model == NULL) {
                return ERROR_MESSAGE(PUMAS_RETURN_VALUE_ERROR,
                    "bad model (expected a string, got nil)");
        }

        int i;
        struct dcs_entry * entry;
        for (i = 0, entry = dcs_stack;
            (i < DCS_STACK_SIZE) && (entry->model != NULL); i++, entry++) {
                if ((entry->process == process) &&
                    (strcmp(entry->model, model) == 0)) {
                        return ERROR_FORMAT(PUMAS_RETURN_MODEL_ERROR,
                            "model %s already registered for %s process",
                            model, process_name[process]);
                }
        }
        if (i == DCS_STACK_SIZE) {
                return ERROR_MESSAGE(PUMAS_RETURN_MEMORY_ERROR,
                    "max stack size reached");
        }

        /* Append the new DCS */
        entry->process = process;
        entry->model = model;
        entry->dcs = dcs;

        return PUMAS_RETURN_SUCCESS;
}

/* API function for getting a DCS model */
enum pumas_return pumas_dcs_get(
    enum pumas_process process, const char * model, pumas_dcs_t ** dcs)
{
        ERROR_INITIALISE(pumas_dcs_get);

        /* Check the process index */
        if (dcs_check_process(process, error_) != PUMAS_RETURN_SUCCESS) {
                return ERROR_RAISE();
        }

        /* Set the default model if none provided */
        if (model == NULL) {
                if (process == PUMAS_PROCESS_BREMSSTRAHLUNG)
                        model = DEFAULT_BREMSSTRAHLUNG;
                else if (process == PUMAS_PROCESS_PAIR_PRODUCTION)
                        model = DEFAULT_PAIR_PRODUCTION;
                else
                        model = DEFAULT_PHOTONUCLEAR;
        }

        /* Look for the model */
        int i;
        struct dcs_entry * entry;
        for (i = 0, entry = dcs_stack;
            (i < DCS_STACK_SIZE) && (entry->model != NULL); i++, entry++) {
                if ((entry->process == process) && 
                    (strcmp(entry->model, model) == 0)) {
                        *dcs = entry->dcs;
                        return PUMAS_RETURN_SUCCESS;
                }
        }
        *dcs = NULL;

        return ERROR_FORMAT(PUMAS_RETURN_MODEL_ERROR,
            "model %s not found for %s process",
            model, process_name[process]);
}

/* API function for getting the name of the default DCS model */
const char * pumas_dcs_default(enum pumas_process process)
{
        if (process == PUMAS_PROCESS_BREMSSTRAHLUNG)
                return DEFAULT_BREMSSTRAHLUNG;
        else if (process == PUMAS_PROCESS_PAIR_PRODUCTION)
                return DEFAULT_PAIR_PRODUCTION;
        else if (process == PUMAS_PROCESS_PHOTONUCLEAR)
                return DEFAULT_PHOTONUCLEAR;
        else
                return NULL;
}
