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

/* The PuMAS API. */
#include "pumas.h"

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
#define PUMAS_VERSION 0
#define PUMAS_SUBVERSION 11

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
 * Relative switch between Continuous Energy Loss (CEL) and DELs.
 */
#define X_FRACTION 5E-02
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
 * Relative step length limit due to Multiple Scattering (MS).
 */
#define RATIO_MSC 1E-02
/**
 * Relative step length limit due to CEL.
 */
#define RATIO_ENERGY_LOSS 1E-02
/**
 * Relative step length limit due to magnetic bending.
 */
#define RATIO_MAGNETIC 1E-02
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
#ifndef M_PI
/**
 * Define pi, if unknown.
 */
#define M_PI 3.14159265358979323846
#endif

/* Helper macro for returning an encapsulated error code. */
#define PUMAS_RETURN(rc, caller)                                               \
        return pumas_return(rc, (pumas_function_t *)caller)

/* Function prototypes for the DCS implementations. */
/**
 * Handle for a DCS computation function.
 */
struct atomic_element;
typedef double(dcs_function_t)(
    const struct atomic_element * element, double K, double q);
/**
 * Handle for a polar angle sampling function.
 */
typedef double(polar_function_t)(
    struct pumas_context * context, double ki, double kf);

/* A collection of low level flags. */
/**
 * Flags for the stepping.
 */
enum stepping_event {
        /** No event occured sofar or none is foreseen. */
        EVENT_NONE = 0,
        /** A kinetic limit occured. */
        EVENT_KINETIC = 1,
        /** A distance or grammage limit is reached or foreseen. */
        EVENT_RANGE = 2,
        /** A grammage limit is reached or foreseen. */
        EVENT_GRAMMAGE = 4,
        /** A proper time limit is reached or foreseen. */
        EVENT_TIME = 8,
        /** A change of medium occured. */
        EVENT_MEDIUM = 16,
        /** A DEL occured or is foreseen. */
        EVENT_DEL = 32,
        /** An EHS occured or is foreseen. */
        EVENT_EHS = 64,
        /** A null weight occured. */
        EVENT_WEIGHT = 128
};
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
 * Tags for operations relative to the the parsing of materials in MDFs.
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
};
/**
 * The local data managed by a simulation context.
 */
struct simulation_context {
        /** The public API settings exposed to the end user. */
        struct pumas_context api;
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
        enum stepping_event step_event;
        /** The expected next event during the stepping. */
        enum stepping_event step_foreseen;
        /** The kinetic limit converted to grammage. */
        double step_X_limit;
        /** The scaterring 1st transport path length of the previous step. */
        double step_invlb1;
        /** The larmor radius of the previous step. */
        double step_rLarmor;
        /** The magnetic transverse direction of the previous step. */
        double step_uT[3];
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
        /*! The memory size left, in bytes. */
        int size;
        /*! Pointer to the next memory segment. */
        struct frame_stack * next;
        /*! Pointer to the first frame in the stack. */
        struct pumas_frame * frame;
        /*! Placeholder for frames */
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
        /** The base material density. */
        double density;
        /** The element-wise weight. */
        double weight;
};
/**
 * Handle for a composite material.
 */
struct composite_material {
        /** The composite reference density. */
        double density;
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
        /** The number of kinetic rows in a dE/dX file. */
        int n_kinetics;
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
struct shared_data {
/*
 * Version tag for the shared data format. Increment whenever the
 * structure changes.
 */
#define BINARY_DUMP_TAG 0

        /** The total byte size of the shared data. */
        int size;
        /** The number of kinetic values in the dE/dX tables. */
        int n_kinetics;
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
        /** The relative electronic density of a material. */
        double * material_ZoA;
        /** The properties of an atomic element . */
        struct atomic_element ** element;
        /** The composition of a base material. */
        struct material_component ** composition;
        /** The composition of a composite material. */
        struct composite_material ** composite;
        /** The material names. */
        char ** material_name;
        /**
         * Placeholder for shared data storage with -double- memory alignment.
         */
        double data[];
};
/**
 * Handle for the shared data, allocated on the heap.
 */
static struct shared_data * s_shared = NULL;

/**
 * Shared data for the error handling.
 */
static struct {
        pumas_handler_cb * handler;
        int catch;
        enum pumas_return catch_rc;
        pumas_function_t * catch_function;
        int line;
#define ERROR_FILE_LENGTH 1024
        char file[ERROR_FILE_LENGTH];
} s_error = { NULL, 0, PUMAS_RETURN_SUCCESS, NULL, 0, "\0" };

/* Prototypes of low level static functions. */
/**
 * Encapsulations of the tabulated CEL and DEL properties.
 */
static double cel_grammage(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double kinetic);
static double cel_grammage_as_time(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double time);
static double cel_proper_time(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double kinetic);
static double cel_kinetic_energy(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double grammage);
static double cel_energy_loss(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double kinetic);
static double cel_magnetic_rotation(
    struct pumas_context * context, int material, double kinetic);
static double del_cross_section(
    struct pumas_context * context, int material, double kinetic);
static double del_interaction_length(
    struct pumas_context * context, int material, double kinetic);
static double del_kinetic_from_interaction_length(
    struct pumas_context * context, int material, double nI);
static double ehs_interaction_length(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double kinetic);
static double ehs_kinetic_from_interaction_length(
    struct pumas_context * context, enum pumas_scheme scheme, int material,
    double nI);
/**
 * Routines related to DCS: implementation and handling.
 */
static inline dcs_function_t * dcs_get(int process);
static inline int dcs_get_index(dcs_function_t * dcs_func);
static double dcs_bremsstrahlung(
    const struct atomic_element * element, double K, double q);
static double dcs_pair_production(
    const struct atomic_element * element, double K, double q);
static double dcs_photonuclear(
    const struct atomic_element * element, double K, double q);
static inline double dcs_photonuclear_d2(
    double A, double K, double q, double Q2);
static inline double dcs_photonuclear_f2_allm(double x, double Q2);
static inline double dcs_photonuclear_f2a_drss(double x, double F2p, double A);
static inline double dcs_photonuclear_r_whitlow(double x, double Q2);
static inline int dcs_photonuclear_check(double K, double q);
static double dcs_ionisation(
    const struct atomic_element * element, double K, double q);
static double dcs_ionisation_integrate(
    int mode, const struct atomic_element * element, double K, double xlow);
static double dcs_ionisation_randomise(struct pumas_context * context,
        const struct atomic_element * element, double K, double xlow);
static double dcs_evaluate(struct pumas_context * context,
    dcs_function_t * dcs_func, const struct atomic_element * element, double K,
    double q);
static void dcs_model_fit(int m, int n, const double * x, const double * y,
    const double * w, double * c);
/**
 * Implementations of polar angle distributions and accessor.
 */
static inline polar_function_t * polar_get(int process);
static double polar_bremsstrahlung(
    struct pumas_context * context, double ki, double kf);
static double polar_pair_production(
    struct pumas_context * context, double ki, double kf);
static double polar_photonuclear(
    struct pumas_context * context, double ki, double kf);
static double polar_ionisation(
    struct pumas_context * context, double ki, double kf);
/**
 * Low level routines for the propagation in matter.
 */
static enum pumas_return transport_with_csda(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    struct medium_locals * locals);
static enum pumas_return transport_csda_deflect(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    struct medium_locals * locals, double ki, double distance);
static enum pumas_return csda_magnetic_transport(struct pumas_context * context,
    int material, double density, double magnet, double charge, double kinetic,
    double phase, double * x, double * y, double * z);
static enum pumas_return transport_with_stepping(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    struct medium_locals * locals, double step_max_medium,
    double step_max_locals);
static double transport_set_locals(struct pumas_medium * medium,
    struct pumas_state * state, struct medium_locals * locals);
static void transport_limit(struct pumas_context * context,
    const struct pumas_state * state, int material, double di, double Xi,
    double * distance_max);
static void transport_do_del(
    struct pumas_context * context, struct pumas_state * state, int material);
static void transport_do_ehs(
    struct pumas_context * context, struct pumas_state * state, int material);
/**
 * Low level routines for randomising DELs.
 */
static polar_function_t * del_randomise_forward(
    struct pumas_context * context, struct pumas_state * state, int material);
static polar_function_t * del_randomise_reverse(
    struct pumas_context * context, struct pumas_state * state, int material);
static void del_randomise_power_law(struct pumas_context * context,
    double alpha, double xmin, double xmax, double * p_r, double * p_w);
static void del_randomise_ziggurat(struct pumas_context * context,
    struct pumas_state * state, dcs_function_t * dcs_func,
    const struct atomic_element * element, double xmin, double xmax,
    float * cdf_sampling);
static void del_randomise_target(struct pumas_context * context,
    struct pumas_state * state, int material, struct del_info * info);
/**
 * Helper routine for recording a state.
 */
static void record_state(struct pumas_recorder * recorder,
    struct pumas_medium * medium, struct pumas_state * state);
/**
 * For memory padding.
 */
static int memory_padded_size(int size, int pad_size);
/**
 * For error handling.
 */
static enum pumas_return pumas_return(
    enum pumas_return rc, pumas_function_t * caller);
/**
 * Routines for the Coulomb scattering and Transverse Transport (TT).
 */
static void coulomb_screening_parameters(struct pumas_context * context,
    double kinetic, int element, double * screening);
static double coulomb_wentzel_path(
    double kinetic, double Z, double A, double screening);
static double coulomb_ehs_length(
    struct pumas_context * context, int material, double kinetic);
static double coulomb_spin_factor(double kinetic);
static void coulomb_frame_parameters(
    double kinetic, double Ma, double * kinetic0, double * parameters);
static void coulomb_pole_decomposition(
    double * screening, double * a, double * b);
static double coulomb_restricted_cs(
    double mu0, double fspin, double * screening, double * a, double * b);
static void coulomb_transport_coefficients(double mu, double fspin,
    double * screening, double * a, double * b, double * coefficient);
static double transverse_transport_ionisation(
    const struct atomic_element * element, double kinetic);
static double transverse_transport_photonuclear(
    const struct atomic_element * element, double kinetic);
/**
 * Routines for handling tables: interpolation and utility accessors.
 */
static void table_bracket(
    const double * table, double value, int * p1, int * p2);
static int table_index(
    struct pumas_context * context, const double * table, double value);
static double table_interpolate(struct pumas_context * context,
    const double * table_X, const double * table_Y, double x);
static void table_get_msc(struct pumas_context * context, int material,
    double kinetic, double * mu0, double * invlb1);
static inline double * table_get_K(int row);
static inline double * table_get_X(int scheme, int material, int row);
static inline double * table_get_T(int scheme, int material, int row);
static inline double * table_get_dE(int scheme, int material, int row);
static inline double * table_get_NI_el(int scheme, int material, int row);
static inline double * table_get_NI_in(int material, int row);
static inline double * table_get_CS(int material, int row);
static inline double * table_get_CSf(int process, int component, int row);
static inline double * table_get_CSn(int process, int element, int row);
static inline double * table_get_Xt(int process, int element, int row);
static inline double * table_get_Kt(int material);
static inline double * table_get_cel(
    int process, int element, int row, double * table);
static inline double * table_get_Li(int material, int order, int row);
static inline double * table_get_a_max(int material);
static inline double * table_get_b_max(int scheme, int material);
static inline double * table_get_Mu0(int material, int row);
static inline double * table_get_Lb(int material, int row);
static inline double * table_get_Ms1(int material, int row);
static inline double * table_get_ms1(int element, int row, double * table);
static inline float * table_get_dcs_coeff(
    const struct atomic_element * element, int process, int kinetic);
static inline float * table_get_dcs_value(
    const struct atomic_element * element, int process, int kinetic);
/**
 * Low level routines for the stepping.
 */
static enum pumas_return step_transport(struct pumas_context * context,
    struct pumas_state * state, int straight, struct pumas_medium * medium,
    struct medium_locals * locals, double grammage_max,
    double * step_max_medium, double * step_max_locals,
    struct pumas_medium ** out_medium);
static void step_fluctuate(struct pumas_context * context,
    struct pumas_state * state, int material, double Xtot, double dX,
    double * kf, double * dE);
static double step_fluctuations2(int material, double kinetic);
static double step_randn(struct pumas_context * context);
static double transport_hard_coulomb_objective(double mu, void * parameters);
static void step_rotate_direction(struct pumas_context * context,
    struct pumas_state * state, double cos_theta);
/**
 * I/O utility routines.
 */
static enum pumas_return io_parse_dedx_file(
    FILE * fid, int material, int * line);
static enum pumas_return io_parse_dedx_row(
    char * buffer, int material, int * row);
static enum pumas_return io_read_line(FILE * fid, char ** buffer);
/**
 * Routines for the parsing of MDFs.
 */
static enum pumas_return mdf_parse_settings(
    struct mdf_buffer * mdf, const char * dedx_path);
static int mdf_settings_index(int operation, int value);
static int mdf_settings_name(int size, char prefix, const char * name);
static enum pumas_return mdf_parse_kinetic(
    struct mdf_buffer * mdf, const char * path);
static enum pumas_return mdf_parse_elements(struct mdf_buffer * mdf);
static enum pumas_return mdf_parse_materials(struct mdf_buffer * mdf);
static enum pumas_return mdf_parse_composites(struct mdf_buffer * mdf);
static enum pumas_return mdf_get_node(
    struct mdf_buffer * mdf, struct mdf_node * node);
static enum pumas_return mdf_skip_pattern(
    struct mdf_buffer * mdf, const char * pattern);
static enum pumas_return mdf_format_path(const char * directory,
    const char * mdf_path, char ** filename, int * offset_dir, int * size_name);
/**
 * Routines for the pre-computation of various properties: CEL, DCS, ...
 */
static enum pumas_return compute_composite(int material);
static enum pumas_return compute_composite_density(int material);
static void compute_composite_weights(int material);
static void compute_composite_tables(int material);
static void compute_cel_integrals(int imed);
static void compute_kinetic_integral(double * table);
static void compute_time_integrals(int material);
static void compute_cel_grammage_integral(int scheme, int material);
static void compute_csda_magnetic_transport(int imed);
static enum pumas_return compute_coulomb_parameters(int medium_index, int row);
static enum pumas_return compute_coulomb_soft(int row, double ** data);
static double compute_cutoff_objective(double mu, void * workspace);
static double * compute_cel_and_del(int row);
static void compute_regularise_del(int material);
static double compute_dcs_integral(int mode,
    const struct atomic_element * element, double kinetic, dcs_function_t * dcs,
    double xlow, int nint);
static void compute_ZoA(int material);
static enum pumas_return compute_dcs_model(
    dcs_function_t * dcs_func, struct atomic_element * element);
/**
 * Helper function for mapping an atomic element from its name.
 */
static int element_index(const char * name);
/**
 * Various math utilities, for integration, root finding and SVD.
 */
static int math_find_root(double (*f)(double x, void * params), double xa,
    double xb, const double * fa_p, const double * fb_p, double xtol,
    double rtol, int iter, void * params, double * x0);
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

/*
 * Low level initialisation. If *dry_mode* is not null the energy loss tables
 * are not loaded and processed.
 */
static enum pumas_return _initialise(enum pumas_particle particle,
    const char * mdf_path, const char * dedx_path, struct pumas_error * error,
    int dry_mode)
{
        /* Reset the error additional info. */
        s_error.line = 0;
        s_error.file[0] = '\0';
        if (error != NULL) {
                error->file = NULL;
                error->line = 0;
        }

        /* Check if the library is already initialised. */
        if (s_shared != NULL)
                PUMAS_RETURN(
                    PUMAS_RETURN_INITIALISATION_ERROR, pumas_initialise);
#if (GDB_MODE)
        /* Save the floating points exceptions status and enable them. */
        fe_status = fegetexcept();
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
        enum pumas_return rc = PUMAS_RETURN_SUCCESS;
        FILE * fid_mdf = NULL;
        struct mdf_buffer * mdf = NULL;
        const int pad_size = sizeof(*(s_shared->data));
#define N_DATA_POINTERS 26
        int size_data[N_DATA_POINTERS];

        /* Check the particle type. */
        if ((particle != PUMAS_PARTICLE_MUON) &&
            (particle != PUMAS_PARTICLE_TAU))
                PUMAS_RETURN(PUMAS_RETURN_UNKNOWN_PARTICLE, pumas_initialise);

        /* Check the path to energy loss tables. */
        if (!dry_mode) {
                if (dedx_path == NULL) dedx_path = getenv("PUMAS_DEDX");
                if (dedx_path == NULL) {
                        rc = PUMAS_RETURN_UNDEFINED_DEDX;
                        goto clean_and_exit;
                }
        }

        /* Parse the MDF. */
        const int size_mdf = 2048;
        const char * file_mdf =
            (mdf_path != NULL) ? mdf_path : getenv("PUMAS_MDF");
        if (file_mdf == NULL) {
                rc = PUMAS_RETURN_UNDEFINED_MDF;
                goto clean_and_exit;
        }
        int size_path = strlen(file_mdf) + 1;
        fid_mdf = fopen(file_mdf, "r");
        if (fid_mdf == NULL) {
                strncpy(s_error.file, file_mdf, ERROR_FILE_LENGTH);
                s_error.file[ERROR_FILE_LENGTH - 1] = '\0';
                rc = PUMAS_RETURN_PATH_ERROR;
                goto clean_and_exit;
        }
        mdf = allocate(size_mdf);
        if (mdf == NULL) {
                rc = PUMAS_RETURN_MEMORY_ERROR;
                goto clean_and_exit;
        }
        mdf->dry_mode = dry_mode;
        mdf->mdf_path = file_mdf;
        mdf->fid = fid_mdf;
        mdf->size = size_mdf - sizeof(*mdf);
        if ((rc = mdf_parse_settings(mdf, dedx_path)) != PUMAS_RETURN_SUCCESS)
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
            memory_padded_size(sizeof(double) * settings.n_kinetics, pad_size);
        /* table_X. */
        size_data[imem++] = memory_padded_size(sizeof(double) * N_SCHEMES *
                settings.n_materials * settings.n_kinetics,
            pad_size);
        /* table_T. */
        size_data[imem++] = memory_padded_size(sizeof(double) * N_SCHEMES *
                settings.n_materials * settings.n_kinetics,
            pad_size);
        /* table_dE. */
        size_data[imem++] = memory_padded_size(sizeof(double) * N_SCHEMES *
                settings.n_materials * settings.n_kinetics,
            pad_size);
        /* table_NI_el. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * 2 * settings.n_materials * settings.n_kinetics,
            pad_size);
        /* table_NI_in. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_kinetics,
            pad_size);
        /* table_CS. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_kinetics,
            pad_size);
        /* table_CSf. */
        size_data[imem++] = memory_padded_size(sizeof(double) *
                N_DEL_PROCESSES * settings.n_components * settings.n_kinetics,
            pad_size);
        /* table_CSn. */
        size_data[imem++] = memory_padded_size(sizeof(double) *
                N_DEL_PROCESSES * settings.n_elements * settings.n_kinetics,
            pad_size);
        /* table_Xt */
        size_data[imem++] = memory_padded_size(sizeof(double) *
                N_DEL_PROCESSES * settings.n_elements * settings.n_kinetics,
            pad_size);
        /* table_Kt. */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials, pad_size);
        /* table_Li. */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials *
                    (N_LARMOR_ORDERS + 1) * settings.n_kinetics,
                pad_size);
        /* table_a_max. */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials, pad_size);
        /* table_b_max. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * N_SCHEMES * settings.n_materials, pad_size);
        /* table_Mu0. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_kinetics,
            pad_size);
        /* table_Lb. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_kinetics,
            pad_size);
        /* table_Ms1. */
        size_data[imem++] = memory_padded_size(
            sizeof(double) * settings.n_materials * settings.n_kinetics,
            pad_size);
        /* elements_in. */
        size_data[imem++] =
            memory_padded_size(sizeof(int) * settings.n_materials, pad_size);
        /* material_ZoA */
        size_data[imem++] =
            memory_padded_size(sizeof(double) * settings.n_materials, pad_size);
        /* element. */
        const int n_dcs = (N_DEL_PROCESSES - 1) *
            (DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1 + DCS_SAMPLING_N) *
            (settings.n_kinetics - settings.dcs_model_offset);
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

        /* Allocate the shared memory. */
        int size_total = 0;
        for (imem = 0; imem < N_DATA_POINTERS; imem++)
                size_total += size_data[imem];
        const int size_shared = sizeof(*s_shared) + size_total;
        void * tmp_ptr = reallocate(mdf, size_shared);
        if (tmp_ptr == NULL) {
                rc = PUMAS_RETURN_MEMORY_ERROR;
                goto clean_and_exit;
        }
        s_shared = tmp_ptr;
        memset(s_shared, 0x0, size_shared);
        s_shared->size = size_shared;
        mdf = NULL;

        /* Map the data pointers. */
        double * p = s_shared->data;
        void ** ptr = (void **)(&(s_shared->mdf_path));
        for (imem = 0; imem < N_DATA_POINTERS; imem++) {
                *ptr = p;
                ptr++;
                p += size_data[imem] / pad_size;
        }

        /* Copy the global settings. */
        s_shared->particle = particle;
        if (particle == PUMAS_PARTICLE_MUON) {
                s_shared->ctau = MUON_C_TAU;
                s_shared->mass = MUON_MASS;
        } else {
                s_shared->ctau = TAU_C_TAU;
                s_shared->mass = TAU_MASS;
        }
        s_shared->n_kinetics = settings.n_kinetics;
        s_shared->n_materials = settings.n_materials;
        s_shared->n_composites = settings.n_composites;
        s_shared->n_elements = settings.n_elements;
        s_shared->n_components = settings.n_components;
        s_shared->max_components = settings.max_components;
        s_shared->n_energy_loss_header = settings.n_energy_loss_header;
        s_shared->dcs_model_offset = settings.dcs_model_offset;
        strcpy(s_shared->mdf_path, file_mdf);

        /* Allocate a new MDF buffer. */
        if ((mdf = allocate(sizeof(struct mdf_buffer) + size_mdf)) == NULL) {
                rc = PUMAS_RETURN_MEMORY_ERROR;
                goto clean_and_exit;
        }
        memcpy(mdf, &settings, sizeof(settings));

        /* Set the path to the dE/dX files. */
        if (!dry_mode)
                strcpy(s_shared->dedx_path, dedx_path);
        else
                s_shared->dedx_path = NULL;

        /* Parse the elements. */
        if ((rc = mdf_parse_elements(mdf)) != PUMAS_RETURN_SUCCESS)
                goto clean_and_exit;

        /* Parse the base materials. */
        if ((rc = mdf_parse_materials(mdf)) != PUMAS_RETURN_SUCCESS)
                goto clean_and_exit;

        /* Parse the composite materials. */
        if ((rc = mdf_parse_composites(mdf)) != PUMAS_RETURN_SUCCESS)
                goto clean_and_exit;

        /* All done if in dry mode. */
        if (dry_mode) goto clean_and_exit;

        /* Precompute the CEL integrals and the TT parameters. */
        int imat;
        for (imat = 0; imat < s_shared->n_materials - s_shared->n_composites;
             imat++) {
                int ikin;
                for (ikin = 0; ikin < s_shared->n_kinetics; ikin++) {
                        rc = compute_coulomb_parameters(imat, ikin);
                        if (rc != PUMAS_RETURN_SUCCESS) goto clean_and_exit;
                }
                compute_cel_integrals(imat);
                compute_csda_magnetic_transport(imat);
        }

        /* Precompute the same properties for composite materials. */
        for (imat = s_shared->n_materials - s_shared->n_composites;
             imat < s_shared->n_materials; imat++) {
                if (((rc = compute_composite(imat)) != PUMAS_RETURN_SUCCESS) ||
                    ((rc = compute_composite_density(imat)) !=
                        PUMAS_RETURN_SUCCESS))
                        goto clean_and_exit;
        }

        /* Tabulate and fit the DCS for atomic elements. */
        int iel;
        for (iel = 0; iel < s_shared->n_elements; iel++) {
                int ip;
                for (ip = 0; ip < N_DEL_PROCESSES - 1; ip++) {
                        if ((rc = compute_dcs_model(
                                 dcs_get(ip), s_shared->element[iel])) !=
                            PUMAS_RETURN_SUCCESS)
                                goto clean_and_exit;
                }
        }

clean_and_exit:
        if (error != NULL) {
                error->file = (s_error.file[0] != '\0') ? s_error.file : NULL;
                error->line = s_error.line;
        }
        if (fid_mdf != NULL) fclose(fid_mdf);
        deallocate(mdf);
        io_read_line(NULL, NULL);
        compute_coulomb_parameters(-1, -1);
        compute_coulomb_soft(-1, NULL);
        compute_cel_and_del(-1);
        compute_dcs_model(NULL, NULL);
        if ((rc != PUMAS_RETURN_SUCCESS) && (s_shared != NULL)) {
                deallocate(s_shared);
                s_shared = NULL;
        }

        PUMAS_RETURN(rc, pumas_initialise);
}

/* The standard API initialisation. */
enum pumas_return pumas_initialise(enum pumas_particle particle,
    const char * mdf_path, const char * dedx_path, struct pumas_error * error)
{
        return _initialise(particle, mdf_path, dedx_path, error, 0);
}

enum pumas_return pumas_load(FILE * stream)
{
        /* Reset the error additional info. */
        s_error.line = 0;
        s_error.file[0] = '\0';

        /* Check if the library is already initialised. */
        if (s_shared != NULL)
                PUMAS_RETURN(PUMAS_RETURN_INITIALISATION_ERROR, pumas_load);
#if (GDB_MODE)
        /* Save the floating points exceptions status and enable them. */
        fe_status = fegetexcept();
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
#endif
        /* Check the binary dump tag. */
        enum pumas_return rc = PUMAS_RETURN_IO_ERROR;
        int tag;
        if (fread(&tag, sizeof(tag), 1, stream) != 1) goto error;
        if (tag != BINARY_DUMP_TAG) {
                rc = PUMAS_RETURN_FORMAT_ERROR;
                goto error;
        }

        /* Allocate the container. */
        int size;
        if (fread(&size, sizeof(size), 1, stream) != 1) goto error;
        s_shared = allocate(size);
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_MEMORY_ERROR;
                goto error;
        }

        /* Load the data and remap the addresses. */
        if (fread(s_shared, size, 1, stream) != 1) goto error;

        void ** ptr = (void **)(&(s_shared->mdf_path));
        ptrdiff_t delta = (char *)(s_shared->data) - (char *)(*ptr);
        int i;
        for (i = 0; i < N_DATA_POINTERS; i++, ptr++)
                *ptr = ((char *)(*ptr)) + delta;

        struct atomic_element ** element = s_shared->element;
        for (i = 0; i < s_shared->n_elements; i++) {
                element[i] =
                    (struct atomic_element *)(((char *)element[i]) + delta);
                element[i]->name += delta;
                element[i]->dcs_data =
                    (float *)(((char *)(element[i]->dcs_data)) + delta);
        }

        struct material_component ** composition = s_shared->composition;
        for (i = 0; i < s_shared->n_materials; i++)
                composition[i] =
                    (struct material_component *)(((char *)composition[i]) +
                        delta);

        struct composite_material ** composite = s_shared->composite;
        for (i = 0; i < s_shared->n_composites; i++)
                composite[i] =
                    (struct composite_material *)(((char *)composite[i]) +
                        delta);

        char ** material_name = s_shared->material_name;
        for (i = 0; i < s_shared->n_materials; i++) material_name[i] += delta;

        return PUMAS_RETURN_SUCCESS;

error:
        deallocate(s_shared);
        s_shared = NULL;
        PUMAS_RETURN(rc, pumas_load);

#undef N_DATA_POINTERS
}

enum pumas_return pumas_dump(FILE * stream)
{
        /* Check if the library is initialised. */
        if (s_shared == NULL)
                PUMAS_RETURN(PUMAS_RETURN_INITIALISATION_ERROR, pumas_dump);

        /* Dump the configuration. */
        int tag = BINARY_DUMP_TAG;
        if (fwrite(&tag, sizeof(tag), 1, stream) != 1) goto error;
        if (fwrite(&s_shared->size, sizeof(s_shared->size), 1, stream) != 1)
                goto error;
        if (fwrite(s_shared, s_shared->size, 1, stream) != 1) goto error;

        return PUMAS_RETURN_SUCCESS;

error:
        PUMAS_RETURN(PUMAS_RETURN_IO_ERROR, pumas_dump);

#undef BINARY_DUMP_TAG
}

void pumas_finalise()
{
#if (GDB_MODE)
        /* Restore the floating points exceptions status. */
        feclearexcept(FE_ALL_EXCEPT);
        feenableexcept(fe_status);
#endif

        if (s_shared == NULL) return;

        /* Free the shared data. */
        int i;
        for (i = 0; i < s_shared->n_materials - s_shared->n_composites; i++) {
                deallocate(s_shared->dedx_filename[i]);
                s_shared->dedx_filename[i] = NULL;
        }
        deallocate(s_shared);
        s_shared = NULL;
}

const char * pumas_error_string(enum pumas_return rc)
{
        static const char * msg[PUMAS_N_RETURNS] = { "Operation succeeded",
                "Unexpected end of file", "Invalid decay mode",
                "Negative or null density value", "Missing node(s) in MDF file",
                "Invalid index value", "Library not/already initialised",
                "An internal library error occured",
                "Couldn't read or write file", "Unexpected format in file",
                "Missing propagation medium", "Not enough memory",
                "Missing user limit", "Missing random callback",
                "No such file or directory", "Raise without catching enabled",
                "MDF node is too long", "Missing energy loss path",
                "Missing MDF file", "Unknown atomic element",
                "Unknown material", "Unknown particle type",
                "Invalid value encountered" };

        if ((rc < PUMAS_RETURN_SUCCESS) || (rc >= PUMAS_N_RETURNS))
                return NULL;
        else
                return msg[rc];
}

const char * pumas_error_function(pumas_function_t * caller)
{
#define TOSTRING(function)                                                     \
        if (caller == (pumas_function_t *)function) return #function;

        /* Library functions with an error code. */
        TOSTRING(pumas_initialise)
        TOSTRING(pumas_dump)
        TOSTRING(pumas_load)
        TOSTRING(pumas_transport)
        TOSTRING(pumas_particle)
        TOSTRING(pumas_context_create)
        TOSTRING(pumas_recorder_create)
        TOSTRING(pumas_material_name)
        TOSTRING(pumas_material_index)
        TOSTRING(pumas_composite_update)
        TOSTRING(pumas_composite_properties)
        TOSTRING(pumas_print)
        TOSTRING(pumas_error_raise)
        TOSTRING(pumas_property_grammage)
        TOSTRING(pumas_property_proper_time)
        TOSTRING(pumas_property_magnetic_rotation)
        TOSTRING(pumas_property_kinetic_energy)
        TOSTRING(pumas_property_energy_loss)
        TOSTRING(pumas_property_scattering_length)
        TOSTRING(pumas_property_cross_section)
        TOSTRING(pumas_table_value)
        TOSTRING(pumas_table_index)

        /* Other library functions. */
        TOSTRING(pumas_finalise)
        TOSTRING(pumas_context_destroy)
        TOSTRING(pumas_recorder_clear)
        TOSTRING(pumas_recorder_destroy)
        TOSTRING(pumas_tag)
        TOSTRING(pumas_error_string)
        TOSTRING(pumas_error_function)
        TOSTRING(pumas_error_handler_set)
        TOSTRING(pumas_error_handler_get)
        TOSTRING(pumas_error_catch)
        TOSTRING(pumas_error_print)
        TOSTRING(pumas_material_length)
        TOSTRING(pumas_composite_length)
        TOSTRING(pumas_table_length)
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
                s_error.catch_rc = PUMAS_RETURN_SUCCESS;
                s_error.catch_function = NULL;
        } else {
                s_error.catch = 0;
        }
}

enum pumas_return pumas_error_raise(void)
{
        if (s_error.catch == 0)
                PUMAS_RETURN(PUMAS_RETURN_RAISE_ERROR, pumas_error_raise);
        s_error.catch = 0;
        PUMAS_RETURN(s_error.catch_rc, s_error.catch_function);
}

void pumas_error_print(FILE * stream, enum pumas_return rc,
    pumas_function_t * function, struct pumas_error * error)
{
        fprintf(stream, "{\"code\" : %d, \"message\" : \"%s\"", rc,
            pumas_error_string(rc));
        if (function != NULL) {
                fprintf(stream, ", \"function\" : \"%s\"",
                    pumas_error_function(function));
        }
        if (error != NULL) {
                if (error->file != NULL)
                        fprintf(stream, ", \"file\" : \"%s\"", error->file);
                if (error->line > 0)
                        fprintf(stream, ", \"line\" : %d", error->line);
        }
        fprintf(stream, "}");
}

/* Public library functions: simulation context management. */
enum pumas_return pumas_context_create(
    int extra_memory, struct pumas_context ** context_)
{
        *context_ = NULL;

        /* Check the library initialisation. */
        if (s_shared == NULL)
                PUMAS_RETURN(
                    PUMAS_RETURN_INITIALISATION_ERROR, pumas_context_create);

        /* Allocate the new context. */
        struct simulation_context * context;
        const int pad_size = sizeof(*(context->data));
        const int work_size =
            memory_padded_size(sizeof(struct coulomb_workspace) +
                    s_shared->max_components * sizeof(struct coulomb_data),
                pad_size);
        if (extra_memory <= 0)
                extra_memory = 0;
        else
                extra_memory = memory_padded_size(extra_memory, pad_size);
        context = allocate(sizeof(*context) + work_size + extra_memory);
        if (context == NULL)
                PUMAS_RETURN(PUMAS_RETURN_MEMORY_ERROR, pumas_context_create);

        /* Set the default configuration. */

        *context_ = (struct pumas_context *)context;
        context->extra_memory = extra_memory;
        if (extra_memory > 0)
                (*context_)->user_data = context->data + work_size / pad_size;
        else
                (*context_)->user_data = NULL;

        int imax = s_shared->n_kinetics - 2;
        context->index_K_last[0] = context->index_K_last[1] = imax;
        context->index_X_last[0] = context->index_X_last[1] = imax;

        (*context_)->medium = NULL;
        (*context_)->random = NULL;
        (*context_)->recorder = NULL;

        (*context_)->longitudinal = 0;
        (*context_)->forward = 1;
        (*context_)->scheme = PUMAS_SCHEME_DETAILED;
        (*context_)->decay = (s_shared->particle == PUMAS_PARTICLE_MUON) ?
            PUMAS_DECAY_WEIGHT :
            PUMAS_DECAY_PROCESS;

        (*context_)->kinetic_limit = 0.; /* GeV */
        (*context_)->distance_max = 0.;  /* m */
        (*context_)->grammage_max = 0.;  /* kg/m^2 */
        (*context_)->time_max = 0.;      /* m/c */

        /* Initialise the Gaussian transform of the random stream. */
        context->randn_done = 0;
        context->randn_next = 0.;

        /* Initialise the work space. */
        context->workspace = (struct coulomb_workspace *)context->data;

        return PUMAS_RETURN_SUCCESS;
}

void pumas_context_destroy(struct pumas_context ** context)
{
        /* Check that the context hasn't already been destroyed. */
        if (*context == NULL) return;

        /* Release the memory. */
        deallocate(*context);
        *context = NULL;
}

/* Public library functions: global print routines. */
enum pumas_return pumas_print(
    FILE * stream, const char * tabulation, const char * newline)
{
        const char * tab = (tabulation == NULL) ? "" : tabulation;
        const char * cr = (newline == NULL) ? "" : newline;

        /* Check the library initialisation. */
        if (s_shared == NULL)
                PUMAS_RETURN(PUMAS_RETURN_INITIALISATION_ERROR, pumas_print);

        /* Print the particle info. */
        if (fprintf(stream,
                "{%s%s\"particle\" : {%s%s%s\"mass (GeV/c^2)\""
                " : %.6lf",
                cr, tab, cr, tab, tab, s_shared->mass) < 0)
                goto error;
        if (fprintf(stream, ",%s%s%s\"lifetime (m/c)\" : %.3lf%s%s}", cr, tab,
                tab, s_shared->ctau, cr, tab) < 0)
                goto error;

        /* Print the atomic elements. */
        if (fprintf(stream, ",%s%s\"elements\" : {", cr, tab) < 0) goto error;
        int iel = 0;
        for (; iel < s_shared->n_elements; iel++) {
                const char * head = (iel == 0) ? "" : ",";
                const struct atomic_element * element = s_shared->element[iel];
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

        /* Print the materials. */
        if (fprintf(stream, ",%s%s\"materials\" : {", cr, tab) < 0) goto error;
        int material = 0;
        for (; material < s_shared->n_materials - s_shared->n_composites;
             material++) {
                const char * head = (material == 0) ? "" : ",";
                if (fprintf(stream, "%s%s%s%s\"%s\" : {", head, cr, tab, tab,
                        s_shared->material_name[material]) < 0)
                        goto error;
                int iel = 0;
                for (; iel < s_shared->elements_in[material]; iel++) {
                        const char * head2 = (iel == 0) ? "" : ",";
                        int element =
                            s_shared->composition[material][iel].element;
                        if (fprintf(stream, "%s%s%s%s%s\"%s (%%)\" : %.5lg",
                                head2, cr, tab, tab, tab,
                                s_shared->element[element]->name,
                                100. *
                                    s_shared->composition[material][iel]
                                        .fraction) < 0)
                                goto error;
                }
                if (fprintf(stream, "%s%s%s}", cr, tab, tab) < 0) goto error;
        }
        if (fprintf(stream, "%s%s}", cr, tab) < 0) goto error;
        if (s_shared->n_composites <= 0) goto closure;

        /* Print the composites. */
        if (fprintf(stream, ",%s%s\"composites\" : {", cr, tab) < 0) goto error;
        const int material0 = s_shared->n_materials - s_shared->n_composites;
        material = material0;
        for (; material < s_shared->n_materials; material++) {
                const char * head = (material == material0) ? "" : ",";
                struct composite_material * composite =
                    s_shared->composite[material - material0];
                if (fprintf(stream, "%s%s%s%s\"%s\" : {", head, cr, tab, tab,
                        s_shared->material_name[material]) < 0)
                        goto error;
                if (fprintf(stream, "%s%s%s%s\"density\" : %.5lg", cr, tab, tab,
                        tab, composite->density * 1E-03) < 0)
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
                                s_shared->material_name[c->material]) < 0)
                                goto error;
                        if (fprintf(stream,
                                "%s%s%s%s%s%s\"density (g/cm^3)\" "
                                ": %.5lg",
                                cr, tab, tab, tab, tab, tab,
                                c->density * 1E-03) < 0)
                                goto error;
                        if (fprintf(stream,
                                ",%s%s%s%s%s%s\"fraction (%%)\" "
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
        PUMAS_RETURN(PUMAS_RETURN_IO_ERROR, pumas_print);
}

int pumas_tag() { return 1000 * PUMAS_VERSION + PUMAS_SUBVERSION; }

/* Public library functions: recorder handling. */
enum pumas_return pumas_recorder_create(
    struct pumas_context * context, struct pumas_recorder ** recorder_)
{
        /*  Allocate memory for the new recorder. */
        struct frame_recorder * recorder = NULL;
        recorder = allocate(sizeof(*recorder));
        if (recorder == NULL)
                PUMAS_RETURN(PUMAS_RETURN_MEMORY_ERROR, pumas_recorder_create);

        /* Configure the context in order to use the new recorder. */
        *recorder_ = (struct pumas_recorder *)recorder;
        if (context != NULL) context->recorder = *recorder_;

        /* configure the new recorder and return it. */
        (*recorder_)->period = 1;
        (*recorder_)->length = 0;
        (*recorder_)->first = NULL;
        recorder->last = NULL;
        recorder->stack = NULL;

        return PUMAS_RETURN_SUCCESS;
}

void pumas_recorder_destroy(struct pumas_recorder ** recorder)
{
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
enum pumas_return pumas_property_grammage(
    enum pumas_scheme scheme, int material, double kinetic, double * grammage)
{
        enum pumas_return rc;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        } else if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
            (scheme >= PUMAS_SCHEME_DETAILED)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        } else if ((material < 0) || (material >= s_shared->n_materials)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        }

        *grammage = cel_grammage(NULL, scheme, material, kinetic);
        return PUMAS_RETURN_SUCCESS;

error:
        *grammage = 0.;
        PUMAS_RETURN(rc, pumas_property_grammage);
}

enum pumas_return pumas_property_proper_time(
    enum pumas_scheme scheme, int material, double kinetic, double * time)
{
        enum pumas_return rc;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        } else if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
            (scheme >= PUMAS_SCHEME_DETAILED)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        } else if ((material < 0) || (material >= s_shared->n_materials)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        }

        *time = cel_proper_time(NULL, scheme, material, kinetic);
        return PUMAS_RETURN_SUCCESS;

error:
        *time = 0.;
        PUMAS_RETURN(rc, pumas_property_proper_time);
}

enum pumas_return pumas_property_magnetic_rotation(
    int material, double kinetic, double * angle)
{
        enum pumas_return rc;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        } else if ((material < 0) || (material >= s_shared->n_materials)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        }

        *angle = cel_magnetic_rotation(NULL, material, kinetic);
        return PUMAS_RETURN_SUCCESS;

error:
        *angle = 0.;
        PUMAS_RETURN(rc, pumas_property_magnetic_rotation);
}

enum pumas_return pumas_property_kinetic_energy(
    enum pumas_scheme scheme, int material, double grammage, double * kinetic)
{
        enum pumas_return rc;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        } else if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
            (scheme >= PUMAS_SCHEME_DETAILED)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        } else if ((material < 0) || (material >= s_shared->n_materials)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        }

        *kinetic = cel_kinetic_energy(NULL, scheme, material, grammage);
        return PUMAS_RETURN_SUCCESS;

error:
        *kinetic = 0.;
        PUMAS_RETURN(rc, pumas_property_kinetic_energy);
}

enum pumas_return pumas_property_energy_loss(
    enum pumas_scheme scheme, int material, double kinetic, double * dedx)
{
        enum pumas_return rc;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        } else if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
            (scheme >= PUMAS_SCHEME_DETAILED)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        } else if ((material < 0) || (material >= s_shared->n_materials)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        }

        *dedx = cel_energy_loss(NULL, scheme, material, kinetic);
        return PUMAS_RETURN_SUCCESS;

error:
        *dedx = 0.;
        PUMAS_RETURN(rc, pumas_property_energy_loss);
}

/* Public library function: elastic scattering 1st transport path length. */
enum pumas_return pumas_property_scattering_length(
    int material, double kinetic, double * length)
{
        *length = 0.;
        enum pumas_return rc;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        } else if ((material < 0) || (material >= s_shared->n_materials)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        }

        double path = 0.;
        int i = 0;
        for (; i < s_shared->elements_in[material]; i++) {
                const struct material_component * component =
                    &s_shared->composition[material][i];
                const struct atomic_element * element =
                    s_shared->element[component->element];
                double kinetic0, screening[3], coefficient[2], fCM[2];
                coulomb_frame_parameters(kinetic, element->A, &kinetic0, fCM);
                coulomb_screening_parameters(
                    NULL, kinetic0, component->element, screening);
                const double fspin = coulomb_spin_factor(kinetic0);
                double a[3], b[3];
                coulomb_pole_decomposition(screening, a, b);
                coulomb_transport_coefficients(
                    1., fspin, screening, a, b, coefficient);
                double d = 1. / (fCM[0] * (1. + fCM[1]));
                d *= d;
                coefficient[1] *= d;
                path += component->fraction /
                    coulomb_wentzel_path(
                        kinetic, element->Z, element->A, screening[0]) *
                    coefficient[1];
        }

        if (path == 0.) {
                rc = PUMAS_RETURN_VALUE_ERROR;
                *length = DBL_MAX;
                goto error;
        }

        *length = 1. / path;
        return PUMAS_RETURN_SUCCESS;

error:
        PUMAS_RETURN(rc, pumas_property_scattering_length);
}

/* Public library function: total inelastic cross section. */
enum pumas_return pumas_property_cross_section(
    int material, double kinetic, double * cross_section)
{
        enum pumas_return rc;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        } else if ((material < 0) || (material >= s_shared->n_materials)) {
                rc = PUMAS_RETURN_INDEX_ERROR;
                goto error;
        }

        *cross_section = (kinetic < *table_get_Kt(material)) ?
            0. :
            del_cross_section(NULL, material, kinetic);
        return PUMAS_RETURN_SUCCESS;

error:
        *cross_section = 0.;
        PUMAS_RETURN(rc, pumas_property_cross_section);
}

/* Public library function: the main transport routine. */
enum pumas_return pumas_transport(
    struct pumas_context * context, struct pumas_state * state)
{
        /* Check the initial state. */
        if (state->decayed) return PUMAS_RETURN_SUCCESS;

        /* Check the configuration. */
        if (s_shared == NULL)
                PUMAS_RETURN(
                    PUMAS_RETURN_INITIALISATION_ERROR, pumas_transport);
        if (context->medium == NULL)
                PUMAS_RETURN(PUMAS_RETURN_MEDIUM_ERROR, pumas_transport);
        if ((s_shared->particle == PUMAS_PARTICLE_TAU) && context->forward &&
            (context->decay == PUMAS_DECAY_WEIGHT))
                PUMAS_RETURN(PUMAS_RETURN_DECAY_ERROR, pumas_transport);

        /* Get the start medium. */
        struct pumas_medium * medium;
        double step_max_medium = context->medium(context, state, &medium);
        if (medium == NULL) {
                /* Register the start of the the track, if recording. */
                if (context->recorder != NULL)
                        record_state(context->recorder, medium, state);
                return PUMAS_RETURN_SUCCESS;
        } else if (step_max_medium > 0.)
                step_max_medium += 0.5 * STEP_MIN;
        struct medium_locals locals = { { 0., { 0., 0., 0. } } };
        const double step_max_locals =
            transport_set_locals(medium, state, &locals);
        if ((step_max_locals > 0.) && (step_max_locals < step_max_medium))
                step_max_medium = step_max_locals;
        if (locals.api.density <= 0.)
                PUMAS_RETURN(PUMAS_RETURN_DENSITY_ERROR, pumas_transport);

        /* Randomise the lifetime, if required. */
        if (context->decay == PUMAS_DECAY_PROCESS) {
                if (context->random == NULL)
                        PUMAS_RETURN(
                            PUMAS_RETURN_MISSING_RANDOM, pumas_transport);
                const double u = context->random(context);
                struct simulation_context * context_ =
                    (struct simulation_context *)context;
                context_->lifetime = state->time - s_shared->ctau * log(u);
        }

        /* Check the transport mode. */
        if ((step_max_medium <= 0.) && (step_max_locals <= 0.)) {
                /* This is an infinite and uniform medium. */
                if ((context->scheme == PUMAS_SCHEME_NO_LOSS) &&
                    (context->distance_max <= 0.) &&
                    (context->grammage_max <= 0.))
                        PUMAS_RETURN(
                            PUMAS_RETURN_MISSING_LIMIT, pumas_transport);
                else if ((context->longitudinal != 0) &&
                    (context->scheme == PUMAS_SCHEME_CSDA)) {
                        /* This is a purely deterministic case. */
                        PUMAS_RETURN(transport_with_csda(
                                         context, state, medium, &locals),
                            pumas_transport);
                }
        }

        /* Transport with a detailed stepping. */
        enum pumas_return rc = transport_with_stepping(
            context, state, medium, &locals, step_max_medium, step_max_locals);

        PUMAS_RETURN(rc, pumas_transport);
}

/* Public library function: transported particle info. */
enum pumas_return pumas_particle(
    enum pumas_particle * particle, double * lifetime, double * mass)
{
        if (s_shared == NULL)
                PUMAS_RETURN(PUMAS_RETURN_INITIALISATION_ERROR, pumas_particle);
        if (particle != NULL) *particle = s_shared->particle;
        if (lifetime != NULL) *lifetime = s_shared->ctau;
        if (mass != NULL) *mass = s_shared->mass;

        return PUMAS_RETURN_SUCCESS;
}

/* Public library functions: materials handling. */
enum pumas_return pumas_material_name(int index, const char ** material)
{
        if (s_shared == NULL)
                PUMAS_RETURN(
                    PUMAS_RETURN_INITIALISATION_ERROR, pumas_material_name);
        if ((index < 0) || (index >= s_shared->n_materials))
                PUMAS_RETURN(PUMAS_RETURN_INDEX_ERROR, pumas_material_name);
        *material = s_shared->material_name[index];

        return PUMAS_RETURN_SUCCESS;
}

enum pumas_return pumas_material_index(const char * material, int * index)
{
        if (s_shared == NULL)
                PUMAS_RETURN(
                    PUMAS_RETURN_INITIALISATION_ERROR, pumas_material_index);
        int i;
        for (i = 0; i < s_shared->n_materials; i++) {
                if (strcmp(s_shared->material_name[i], material) == 0) {
                        *index = i;
                        return PUMAS_RETURN_SUCCESS;
                }
        }

        PUMAS_RETURN(PUMAS_RETURN_UNKNOWN_MATERIAL, pumas_material_index);
}

int pumas_material_length(void)
{
        if (s_shared == NULL) return 0;
        return s_shared->n_materials;
}

/* Public library functions: accessing and modifying composite materials. */
int pumas_composite_length(void)
{
        if (s_shared == NULL) return 0;
        return s_shared->n_composites;
}

enum pumas_return pumas_composite_update(
    int material, const double * fractions, const double * densities)
{
        if (s_shared == NULL)
                PUMAS_RETURN(
                    PUMAS_RETURN_INITIALISATION_ERROR, pumas_composite_update);

        const int i0 = s_shared->n_materials - s_shared->n_composites;
        if ((material < i0) || (material > s_shared->n_materials - 1))
                PUMAS_RETURN(PUMAS_RETURN_INDEX_ERROR, pumas_composite_update);

        const int icomp =
            material - s_shared->n_materials + s_shared->n_composites;
        int i;
        for (i = 0; i < s_shared->composite[icomp]->n_components; i++) {
                struct composite_component * component =
                    s_shared->composite[icomp]->component + i;
                if (fractions != NULL) component->fraction = fractions[i];
                if (densities != NULL) component->density = densities[i];
        }

        enum pumas_return rc;
        rc = compute_composite_density(material);
        if ((rc != PUMAS_RETURN_SUCCESS) || (fractions == NULL))
                goto clean_and_exit;
        rc = compute_composite(material);

clean_and_exit:
        /* Free temporary workspace and return. */
        compute_coulomb_parameters(-1, -1);
        PUMAS_RETURN(rc, pumas_composite_update);
}

enum pumas_return pumas_composite_properties(int material, double * density,
    int * components, double * fractions, double * densities)
{
        if (s_shared == NULL)
                PUMAS_RETURN(
                    PUMAS_RETURN_INITIALISATION_ERROR, pumas_composite_update);

        const int i0 = s_shared->n_materials - s_shared->n_composites;
        if ((material < i0) || (material > s_shared->n_materials - 1))
                PUMAS_RETURN(
                    PUMAS_RETURN_INDEX_ERROR, pumas_composite_properties);

        const int icomp =
            material - s_shared->n_materials + s_shared->n_composites;
        if (density != NULL) *density = s_shared->composite[icomp]->density;
        if (components != NULL)
                *components = s_shared->composite[icomp]->n_components;
        int i;
        for (i = 0; i < s_shared->composite[icomp]->n_components; i++) {
                struct composite_component * component =
                    s_shared->composite[icomp]->component + i;
                if (fractions != NULL) fractions[i] = component->fraction;
                if (densities != NULL) densities[i] = component->density;
        }

        return PUMAS_RETURN_SUCCESS;
}

/* Public library functions: info on tabulations. */
enum pumas_return pumas_table_value(enum pumas_property property,
    enum pumas_scheme scheme, int material, int row, double * value)
{
        /* Check the input parameters. */
        enum pumas_return rc;
        *value = 0.;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        }

        rc = PUMAS_RETURN_INDEX_ERROR;
        if ((material < 0) || (material >= s_shared->n_materials))
                goto error;
        else if ((row < 0) || (row >= s_shared->n_kinetics))
                goto error;

        if (property == PUMAS_PROPERTY_KINETIC_ENERGY) {
                *value = *table_get_K(row);
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_GRAMMAGE) {
                if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
                    (scheme >= PUMAS_SCHEME_DETAILED))
                        goto error;
                *value = *table_get_X(scheme, material, row);
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_PROPER_TIME) {
                if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
                    (scheme >= PUMAS_SCHEME_DETAILED))
                        goto error;
                *value = *table_get_T(scheme, material, row);
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_ENERGY_LOSS) {
                if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
                    (scheme >= PUMAS_SCHEME_DETAILED))
                        goto error;
                *value = *table_get_dE(scheme, material, row);
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_MAGNETIC_ROTATION) {
                const double factor = LARMOR_FACTOR / s_shared->mass;
                *value =
                    *table_get_T(PUMAS_SCHEME_CSDA, material, row) * factor;
                return PUMAS_RETURN_SUCCESS;
        } else if (property == PUMAS_PROPERTY_CROSS_SECTION) {
                *value = *table_get_CS(material, row);
                return PUMAS_RETURN_SUCCESS;
        }

error:
        PUMAS_RETURN(rc, pumas_table_value);
}

int pumas_table_length(void)
{
        if (s_shared == NULL) return 0;
        return s_shared->n_kinetics;
}

enum pumas_return pumas_table_index(enum pumas_property property,
    enum pumas_scheme scheme, int material, double value, int * index)
{
        const int imax = s_shared->n_kinetics - 1;
        const double * table;
        enum pumas_return rc;

        /* Check some input parameters. */
        *index = 0;
        if (s_shared == NULL) {
                rc = PUMAS_RETURN_INITIALISATION_ERROR;
                goto error;
        }

        rc = PUMAS_RETURN_INDEX_ERROR;
        if ((material < 0) || (material >= s_shared->n_materials)) goto error;

        /* Get the tabulated value's index. */
        if (property == PUMAS_PROPERTY_KINETIC_ENERGY)
                table = table_get_K(0);
        else if (property == PUMAS_PROPERTY_GRAMMAGE) {
                if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
                    (scheme >= PUMAS_SCHEME_DETAILED))
                        goto error;
                table = table_get_X(scheme, material, 0);
        } else if (property == PUMAS_PROPERTY_PROPER_TIME) {
                if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
                    (scheme >= PUMAS_SCHEME_DETAILED))
                        goto error;
                table = table_get_T(scheme, material, 0);
        } else if (property == PUMAS_PROPERTY_ENERGY_LOSS) {
                if ((scheme <= PUMAS_SCHEME_NO_LOSS) ||
                    (scheme >= PUMAS_SCHEME_DETAILED))
                        goto error;
                table = table_get_dE(scheme, material, 0);
        } else if (property == PUMAS_PROPERTY_MAGNETIC_ROTATION) {
                value *= s_shared->mass / LARMOR_FACTOR;
                table = table_get_T(PUMAS_SCHEME_CSDA, material, 0);
        } else if (property == PUMAS_PROPERTY_CROSS_SECTION)
                table = table_get_CS(material, 0);
        else
                goto error;

        rc = PUMAS_RETURN_VALUE_ERROR;
        if (value < table[0]) {
                goto error;
        } else if (value >= table[imax]) {
                *index = imax;
                goto error;
        }

        int i1 = 0, i2 = imax;
        table_bracket(table, value, &i1, &i2);
        *index = i1;

        return PUMAS_RETURN_SUCCESS;
error:
        PUMAS_RETURN(rc, pumas_table_index);
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
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the total grammage in kg/m^2 otherwise `0`.
 */
double cel_grammage(struct pumas_context * context, enum pumas_scheme scheme,
    int material, double kinetic)
{
        const int imax = s_shared->n_kinetics - 1;
        if (kinetic < *table_get_K(0)) return 0.;

        if (kinetic >= *table_get_K(imax)) {
                /* Constant energy loss model. */
                const double a = *table_get_a_max(material);
                const double b = *table_get_b_max(scheme, material);
                const double K0 = *table_get_K(imax);
                const double K1 = a / b + s_shared->mass;

                return *table_get_X(scheme, material, imax) +
                    1. / b * log((kinetic + K1) / (K0 + K1));
        }

        /* Interpolation. */
        return table_interpolate(
            context, table_get_K(0), table_get_X(scheme, material, 0), kinetic);
}

/**
 * Total grammage for a deterministic CEL as function of total proper time.
 *
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param time     The total proper time variation.
 * @return On success, the total grammage in kg/m^2 otherwise `0`.
 */
double cel_grammage_as_time(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double time)
{
        const int imax = s_shared->n_kinetics - 1;
        if (time < *table_get_T(scheme, material, 0)) return 0.;

        if (time >= *table_get_T(scheme, material, imax)) {
                /* Constant energy loss model. */
                const double a = *table_get_a_max(material);
                const double b = *table_get_b_max(scheme, material);
                const double E0 = *table_get_K(imax) + s_shared->mass;
                const double tmax = *table_get_T(scheme, material, imax);
                const double r =
                    E0 / (a + b * E0) * exp(a * (time - tmax) / s_shared->mass);
                const double kinetic = a * r / (1. - b * r) - s_shared->mass;
                const double K1 = a / b + s_shared->mass;
                return *table_get_X(scheme, material, imax) +
                    1. / b * log((kinetic + K1) / (E0 - s_shared->mass + K1));
        }

        /* Interpolation. */
        return table_interpolate(context, table_get_T(scheme, material, 0),
            table_get_X(scheme, material, 0), time);
}

/**
 * Total proper time for a deterministic CEL in a homogeneous medium.
 *
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the normalised proper time in kg/m^2 otherwise `0`.
 */
double cel_proper_time(struct pumas_context * context, enum pumas_scheme scheme,
    int material, double kinetic)
{
        const int imax = s_shared->n_kinetics - 1;
        if (kinetic < *table_get_K(0)) return 0.;

        if (kinetic >= *table_get_K(imax)) {
                /* Constant energy loss model. */
                const double a = *table_get_a_max(material);
                const double b = *table_get_b_max(scheme, material);
                const double E0 = *table_get_K(imax) + s_shared->mass;
                const double E1 = kinetic + s_shared->mass;

                return *table_get_T(scheme, material, imax) +
                    s_shared->mass / a *
                    log((E1 / E0) * (a + b * E0) / (a + b * E1));
        }

        /* Interpolation. */
        return table_interpolate(
            context, table_get_K(0), table_get_T(scheme, material, 0), kinetic);
}

/**
 * The initial kinetic energy for a given total grammage assuming a determistic
 * CEL.
 *
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param grammage The total grammage depth.
 * @return On success, the initial kinetic energy in GeV otherwise `0`.
 */
double cel_kinetic_energy(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double grammage)
{
        const int imax = s_shared->n_kinetics - 1;
        if (grammage < *table_get_X(scheme, material, 0)) return 0.;

        if (grammage >= *table_get_X(scheme, material, imax)) {
                /* Constant energy loss model. */
                const double a = *table_get_a_max(material);
                const double b = *table_get_b_max(scheme, material);
                const double K0 = *table_get_K(imax);
                const double K1 = a / b + s_shared->mass;

                return (K0 + K1) *
                    exp(b * (grammage - *table_get_X(scheme, material, imax))) -
                    K1;
        }

        /* Interpolation. */
        return table_interpolate(context, table_get_X(scheme, material, 0),
            table_get_K(0), grammage);
}

/**
 * The average CEL.
 *
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the CEL in GeV/(kg/m^2) otherwise `0`.
 */
double cel_energy_loss(struct pumas_context * context, enum pumas_scheme scheme,
    int material, double kinetic)
{
        const int imax = s_shared->n_kinetics - 1;
        if (kinetic < *table_get_K(0)) return 0.;

        if (kinetic >= *table_get_K(imax)) {
                /* Constants energy loss model. */
                return *table_get_a_max(material) +
                    *table_get_b_max(scheme, material) *
                    (kinetic + s_shared->mass);
        }

        /* Interpolation. */
        return table_interpolate(context, table_get_K(0),
            table_get_dE(scheme, material, 0), kinetic);
}

/**
 * The normalised magnetic rotation angle for a deterministic CEL.
 *
 * @param context  The simulation context.
 * @param material The index of the material in which the
 *                 particle travels.
 * @param kinetic  The initial kinetic energy, in GeV.
 * @return The normalised rotation angle in kg/m^2/T.
 *
 * The magnetic rotation angle is proportional to the proper time integral.
 * Therefore it is computed from the proper time table.
 */
double cel_magnetic_rotation(
    struct pumas_context * context, int material, double kinetic)
{
        const int imax = s_shared->n_kinetics - 1;
        const double factor = LARMOR_FACTOR / s_shared->mass;
        double * const T = table_get_T(PUMAS_SCHEME_CSDA, material, 0);
        if (kinetic <= *table_get_K(0)) return T[imax] * factor;

        if (kinetic >= *table_get_K(imax)) {
                /*
                 * Neglect any magnetic rotation above the max tabulated value
                 * of the kinetic energy.
                 */
                return 0.;
        }

        /* Interpolation. */
        return (T[imax] -
                   table_interpolate(context, table_get_K(0), T, kinetic)) *
            factor;
}

/**
 * The total cross section for inelastic DELs.
 *
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the DEL total cross section, otherwise `0`.
 */
double del_cross_section(
    struct pumas_context * context, int material, double kinetic)
{
        const int imax = s_shared->n_kinetics - 1;
        if (kinetic < *table_get_K(0)) return 0.;

        if (kinetic >= *table_get_K(imax)) {
                /* Constant cross section model. */
                return *table_get_CS(material, imax);
        }

        /* Interpolation. */
        return table_interpolate(
            NULL, table_get_K(0), table_get_CS(material, 0), kinetic);
}

/**
 * Total number of interaction lengths for a deterministic CEL.
 *
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the number of interaction lenths, otherwise `0`.
 */
double del_interaction_length(
    struct pumas_context * context, int material, double kinetic)
{
        const int imax = s_shared->n_kinetics - 1;
        if (kinetic < *table_get_K(0)) return 0.;

        if (kinetic >= *table_get_K(imax)) {
                /* constant loss model. */
                const double k0 = *table_get_K(imax);
                const double a0 = *table_get_a_max(material);
                const double b0 =
                    *table_get_b_max(PUMAS_SCHEME_HYBRID, material);
                const double cs = *table_get_CS(material, imax);
                const double dZ = cs / b0 *
                    log((a0 + b0 * (kinetic + s_shared->mass)) /
                        (a0 + b0 * (k0 + s_shared->mass)));
                return *table_get_NI_in(material, imax) + dZ;
        }

        /* Interpolation. */
        return table_interpolate(
            context, table_get_K(0), table_get_NI_in(material, 0), kinetic);
}

/**
 * Initial kinetic energy for a given number of interaction lengths.
 *
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param nI       The number of interaction lengths.
 * @return On success, the initial kinetic energy, otherwise `0`.
 */
double del_kinetic_from_interaction_length(
    struct pumas_context * context, int material, double nI)
{
        const int imax = s_shared->n_kinetics - 1;
        if (nI < *table_get_NI_in(material, 0)) return 0.;

        if (nI >= *table_get_NI_in(material, imax)) {
                /* constant loss model. */
                const double k0 = *table_get_K(imax);
                const double a0 = *table_get_a_max(material);
                const double b0 =
                    *table_get_b_max(PUMAS_SCHEME_HYBRID, material);
                const double cs = *table_get_CS(material, imax);
                const double nI0 = *table_get_NI_in(material, imax);
                return ((a0 + b0 * (k0 + s_shared->mass)) *
                               exp(b0 * (nI - nI0) / cs) -
                           a0) /
                    b0 -
                    s_shared->mass;
        }

        /* Interpolation. */
        return table_interpolate(
            context, table_get_NI_in(material, 0), table_get_K(0), nI);
}

/**
 * Total number of EHS interaction lengths for a deterministic CEL.
 *
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param kinetic  The initial kinetic energy.
 * @return On success, the number of interaction lengths, otherwise `0`.
 */
double ehs_interaction_length(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double kinetic)
{
        const int imax = s_shared->n_kinetics - 1;
        if (kinetic < *table_get_K(0)) return 0.;

        if (kinetic >= *table_get_K(imax)) {
                /* linear extrapolation. */
                const double k0 = *table_get_K(imax - 1);
                const double k1 = *table_get_K(imax);
                const double nI0 = *table_get_NI_el(scheme, material, imax - 1);
                const double nI1 = *table_get_NI_el(scheme, material, imax);
                return nI1 + (kinetic - k1) * (nI1 - nI0) / (k1 - k0);
        }

        /* Interpolation. */
        return table_interpolate(context, table_get_K(0),
            table_get_NI_el(scheme, material, 0), kinetic);
}

/**
 * Initial kinetic energy for a total number of EHS interaction lengths.
 *
 * @param context  The simulation context.
 * @param scheme   The energy loss scheme.
 * @param material The index of the propagation material.
 * @param nI       The number of interaction lengths.
 * @return On success, the initial kinetic energy, otherwise `0`.
 */
double ehs_kinetic_from_interaction_length(struct pumas_context * context,
    enum pumas_scheme scheme, int material, double nI)
{
        const int imax = s_shared->n_kinetics - 1;
        if (nI < *table_get_NI_el(scheme, material, 0)) return 0.;

        if (nI >= *table_get_NI_el(scheme, material, imax)) {
                /* linear extrapolation. */
                const double k0 = *table_get_K(imax - 1);
                const double k1 = *table_get_K(imax);
                const double nI0 = *table_get_NI_el(scheme, material, imax - 1);
                const double nI1 = *table_get_NI_el(scheme, material, imax);
                return k1 + (nI - nI1) * (k1 - k0) / (nI1 - nI0);
        }

        /* Interpolation. */
        return table_interpolate(
            context, table_get_NI_el(scheme, material, 0), table_get_K(0), nI);
}

/*
 * Low level routines: generic and specific interpolations of various
 * tabulated data.
 */
/**
 * Interpolation of the Multiple SCattering (MSC) parameters.
 *
 * @param context  The simulation context.
 * @param material The propagation material.
 * @param kinetic  The kinetic energy.
 * @param mu0      The interpolated angular cuttof for Coulomb scattering.
 * @param invlb1   The interpolated 1st transport inverse grammage.
 *
 * Interpolate the cutt-off angle for Coulomb scattering and the total 1st
 * transport path length for MSC.
 */
void table_get_msc(struct pumas_context * context, int material, double kinetic,
    double * mu0, double * invlb1)
{
        const int imax = s_shared->n_kinetics - 1;
        if (kinetic < *table_get_K(1)) {
                *mu0 = *table_get_Mu0(material, 1);
                /* Use asymptotic limit as lb1 ~ sqrt(kinetic). */
                *invlb1 = *table_get_Ms1(material, 1) *
                    sqrt((*table_get_K(1)) / kinetic);
        } else if (kinetic >= *table_get_K(imax)) {
                *mu0 = *table_get_Mu0(material, imax);
                /* Use asymptotic limit as lb1 ~ kinetic. */
                *invlb1 = *table_get_Ms1(material, imax) *
                    (*table_get_K(imax)) / kinetic;
        } else {
                const int i1 = table_index(context, table_get_K(0), kinetic);
                const int i2 = i1 + 1;
                double h = (kinetic - *table_get_K(i1)) /
                    (*table_get_K(i2) - *table_get_K(i1));
                *mu0 = *table_get_Mu0(material, i1) +
                    h *
                        (*table_get_Mu0(material, i2) -
                            *table_get_Mu0(material, i1));
                *invlb1 = *table_get_Ms1(material, i1) +
                    h *
                        (*table_get_Ms1(material, i2) -
                            *table_get_Ms1(material, i1));
        }
}

/**
 * Interpolate an arbitrary tabulated property.
 *
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
double table_interpolate(struct pumas_context * context, const double * table_X,
    const double * table_Y, double x)
{
        const int i1 = table_index(context, table_X, x);
        const int i2 = i1 + 1;
        double h = (x - table_X[i1]) / (table_X[i2] - table_X[i1]);
        return table_Y[i1] + h * (table_Y[i2] - table_Y[i1]);
}

/**
 * Find the index closest to `value`, from below.
 *
 * @param context The simulation context.
 * @param table   The tabulated values.
 * @param value   The value to bracket.
 * @return In case of success the closest index from below is returned,
 * otherwise -1 in case of underflow or `shared_data::n_kinetics-1` in
 * case of overflow.
 *
 * Compute the table index for the given entry `value` using a dichotomy
 * algorithm. If a `context` is not `NULL`, `value` is checked against the
 * last used indices in the table, before doing the dichotomy search.
 */
int table_index(
    struct pumas_context * context, const double * table, double value)
{
        int * last = NULL;
        if (context != NULL) {
                /* Check if the last used indices are still relevant. */
                struct simulation_context * const ctx =
                    (struct simulation_context * const)context;
                if (table == table_get_K(0))
                        last = ctx->index_K_last;
                else
                        last = ctx->index_X_last;

                if ((value >= table[last[0]]) && (value < table[last[0] + 1]))
                        return last[0];
                if ((value >= table[last[1]]) && (value < table[last[1] + 1]))
                        return last[1];
        }

        /* Check the boundaries. */
        const int imax = s_shared->n_kinetics - 1;
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
 * @param row The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_K(int row) { return s_shared->table_K + row; }

/**
 * Encapsulation of the total grammage table.
 *
 * @param scheme   The energy loss scheme.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_X(int scheme, int material, int row)
{
        scheme = (scheme > PUMAS_SCHEME_HYBRID) ? PUMAS_SCHEME_HYBRID : scheme;
        return s_shared->table_X +
            (scheme * s_shared->n_materials + material) * s_shared->n_kinetics +
            row;
}

/**
 * Encapsulation of the proper time table.
 *
 * @param scheme   The energy loss scheme.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_T(int scheme, int material, int row)
{
        scheme = (scheme > PUMAS_SCHEME_HYBRID) ? PUMAS_SCHEME_HYBRID : scheme;
        return s_shared->table_T +
            (scheme * s_shared->n_materials + material) * s_shared->n_kinetics +
            row;
}

/**
 * Encapsulation of the dE/dX table.
 *
 * @param scheme   The energy loss scheme.
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_dE(int scheme, int material, int row)
{
        scheme = (scheme > PUMAS_SCHEME_HYBRID) ? PUMAS_SCHEME_HYBRID : scheme;
        return s_shared->table_dE +
            (scheme * s_shared->n_materials + material) * s_shared->n_kinetics +
            row;
}

/**
 * Encapsulation of the number of EHS interaction lengths.
 *
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_NI_el(int scheme, int material, int row)
{
        scheme = (scheme >= PUMAS_SCHEME_HYBRID) ? PUMAS_SCHEME_HYBRID :
                                                   PUMAS_SCHEME_CSDA;
        return s_shared->table_NI_el +
            (scheme * s_shared->n_materials + material) * s_shared->n_kinetics +
            row;
}

/**
 * Encapsulation of the number of interaction lengths for inelastic DELs.
 *
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_NI_in(int material, int row)
{
        return s_shared->table_NI_in + material * s_shared->n_kinetics + row;
}

/**
 * Encapsulation of the cross section table, for inelastic DELs.
 *
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_CS(int material, int row)
{
        return s_shared->table_CS + material * s_shared->n_kinetics + row;
}

/**
 * Encapsulation of the fractional cross-sections table, for inelastic DELs.
 *
 * @param process   The process index.
 * @param component The component index.
 * @param row       The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_CSf(int process, int component, int row)
{
        return s_shared->table_CSf +
            (process * s_shared->n_components + component) *
            s_shared->n_kinetics +
            row;
}

/**
 * Encapsulation of the cross-sections normalisation table.
 *
 * @param process The process index.
 * @param element The element index.
 * @param row     The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_CSn(int process, int element, int row)
{
        return s_shared->table_CSn +
            (process * s_shared->n_elements + element) * s_shared->n_kinetics +
            row;
}

/**
 * Encapsulation of the fractional lower threshold for inelastic DELs.
 *
 * @param process The process index.
 * @param element The element index.
 * @param row     The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_Xt(int process, int element, int row)
{
        return s_shared->table_Xt +
            (process * s_shared->n_elements + element) * s_shared->n_kinetics +
            row;
}

/**
 * Encapsulation of the lower kinetic threshold for inelatic DELs.
 *
 * @param material The material index.
 * @return A pointer to the table element.
 */
double * table_get_Kt(int material) { return s_shared->table_Kt + material; }

/**
 * Encapsulation of the temporary average CEL table, element wise.
 *
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
double * table_get_cel(int process, int element, int row, double * table)
{
        return table +
            (process * s_shared->n_elements + element) * s_shared->n_kinetics +
            row;
}

/**
 * Encapsulation of the magnetic deflection table, when using CSDA in a
 * homogeneous medium of infinite extension.
 *
 * @param material The material index.
 * @param order    The order of the Taylor expansion.
 * @param row      The kinetic row index.
 * @return A pointer to the table element.
 */
double * table_get_Li(int material, int order, int row)
{
        return s_shared->table_Li +
            (material * N_LARMOR_ORDERS + order) * s_shared->n_kinetics + row;
}

/**
 * Encapsulation of the last tabulated ionisation dE/dX.
 *
 * @param material The material index.
 * @return A pointer to the table element.
 */
double * table_get_a_max(int material)
{
        return s_shared->table_a_max + material;
}

/**
 * Encapsulation of the last tabulated radiative dE/dX parameter.
 *
 * @param scheme   The energy loss scheme.
 * @param material The material index.
 * @return A pointer to the table element.
 */
double * table_get_b_max(int scheme, int material)
{
        scheme = (scheme > PUMAS_SCHEME_HYBRID) ? PUMAS_SCHEME_HYBRID : scheme;
        return s_shared->table_b_max + scheme * s_shared->n_materials +
            material;
}

/**
 * Encapsulation of the angular cutoff for Coulomb scattering.
 *
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_Mu0(int material, int row)
{
        return s_shared->table_Mu0 + material * s_shared->n_kinetics + row;
}

/**
 * Encapsulation of the interaction length for EHS events.
 *
 * @param material The material index.
 * @param row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_Lb(int material, int row)
{
        return s_shared->table_Lb + material * s_shared->n_kinetics + row;
}

/*!
 * Encapsulation of the Multiple SCattering (MSC) 1st transport path length.
 *
 * @param [in] material The material index.
 * @param [in] row      The kinetic energy row index.
 * @return A pointer to the table element.
 */
double * table_get_Ms1(int material, int row)
{
        return s_shared->table_Ms1 + material * s_shared->n_kinetics + row;
}

/*!
 * Encapsulation of the temporary MSC 1st transport path length, element wise.
 *
 * @param [in] element  The material index.
 * @param [in] row      The kinetic energy row index.
 * @param [in] table    The temporary table.
 * @return A pointer to the table element.
 */
double * table_get_ms1(int element, int row, double * table)
{
        return table + element * s_shared->n_kinetics + row;
}

/**
 * Encapsulation of the polynomial coefficients of the DCS model.
 *
 * @param element The element index.
 * @param process The process index.
 * @param row     The kinetic energy row index.
 * @return A pointer to the table element.
 */
float * table_get_dcs_coeff(
    const struct atomic_element * element, int process, int kinetic)
{
        const int n =
            DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1 + DCS_SAMPLING_N;
        return element->dcs_data +
            (process * (s_shared->n_kinetics - s_shared->dcs_model_offset) +
                kinetic) *
            n;
}

/**
 * Encapsulation of the tabulated DCS values.
 *
 * @param element The element index.
 * @param process The process index.
 * @param row     The kinetic energy row index.
 * @return A pointer to the table element.
 */
float * table_get_dcs_value(
    const struct atomic_element * element, int process, int kinetic)
{
        const int n = DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1;
        return element->dcs_data +
            (process * (s_shared->n_kinetics - s_shared->dcs_model_offset) +
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
 * @param context      The simulation context.
 * @param state        The initial/final state.
 * @param medium_index The index of the propagation medium.
 * @param locals       Handle for the local properties of the uniform medium.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
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
enum pumas_return transport_with_csda(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    struct medium_locals * locals)
{
        /* Unpack and backup the initial state. */
        const int material = medium->material;
        const double density = locals->api.density;
        const double ki = state->kinetic;
        const double di = state->distance;
        const double ti = state->time;
        const double xi =
            cel_grammage(context, PUMAS_SCHEME_CSDA, material, ki);
        const double Ti =
            cel_proper_time(context, PUMAS_SCHEME_CSDA, material, ki);

        /* Register the start of the the track, if recording. */
        struct pumas_recorder * recorder = context->recorder;
        const int record = (recorder != NULL);
        if (record) record_state(recorder, 0, state);

        /* Get the end point with CSDA. */
        enum stepping_event event = EVENT_NONE;
        double xB;
        if (context->grammage_max > 0.) {
                xB = context->grammage_max - state->grammage;
                event = EVENT_GRAMMAGE;
        } else
                xB = DBL_MAX;
        if (context->distance_max > 0.) {
                const double xD = density * (context->distance_max - di);
                if (xD < xB) {
                        xB = xD;
                        event = EVENT_RANGE;
                }
        }
        int decayed = 0;
        double time_max = context->time_max;
        if (context->decay == PUMAS_DECAY_PROCESS) {
                struct simulation_context * c =
                    (struct simulation_context *)context;
                if ((time_max <= 0.) || (c->lifetime < time_max)) {
                        time_max = c->lifetime;
                        decayed = 1;
                }
        }
        if (time_max > 0.) {
                const double dt = time_max - ti;
                if (dt <= 0.)
                        xB = 0.;
                else {
                        const double sgn = context->forward ? 1. : -1.;
                        const double Tf = Ti - sgn * dt * density;
                        if (Tf > 0.) {
                                const double xT = fabs(xi -
                                    cel_grammage_as_time(context,
                                        PUMAS_SCHEME_CSDA, material, Tf));
                                if (xT < xB) {
                                        xB = xT;
                                        event = EVENT_TIME;
                                }
                        }
                }
        }
        if (xB <= 0.) return PUMAS_RETURN_SUCCESS;

        double xf, kf;
        if (context->forward != 0) {
                if (xi > xB) {
                        xf = xi - xB;
                        kf = cel_kinetic_energy(
                            context, PUMAS_SCHEME_CSDA, material, xf);
                } else {
                        xf = kf = 0.;
                }

                if (context->kinetic_limit > 0.) {
                        if (ki <= context->kinetic_limit)
                                return PUMAS_RETURN_SUCCESS;
                        if (kf < context->kinetic_limit) {
                                kf = context->kinetic_limit;
                                xf = cel_grammage(
                                    context, PUMAS_SCHEME_CSDA, material, kf);
                                event = EVENT_KINETIC;
                        }
                }
        } else {
                xf = xB + xi;
                kf = cel_kinetic_energy(
                    context, PUMAS_SCHEME_CSDA, material, xf);
                if (context->kinetic_limit > 0.) {
                        if (ki >= context->kinetic_limit)
                                return PUMAS_RETURN_SUCCESS;
                        if (kf > context->kinetic_limit) {
                                kf = context->kinetic_limit;
                                xf = cel_grammage(
                                    context, PUMAS_SCHEME_CSDA, material, kf);
                                event = EVENT_KINETIC;
                        }
                }
        }

        /* Update the end point statistics. */
        const double distance = fabs(xf - xi) / density;
        state->kinetic = kf;
        if (event == EVENT_RANGE)
                state->distance = context->distance_max;
        else
                state->distance += distance;
        if (event == EVENT_GRAMMAGE)
                state->grammage = context->grammage_max;
        else
                state->grammage += distance * density;
        if (event == EVENT_TIME) {
                state->time = time_max;
                state->decayed = decayed;
        } else
                state->time += fabs(Ti -
                                   cel_proper_time(context, PUMAS_SCHEME_CSDA,
                                       material, kf)) /
                    density;
        if (!context->forward)
                state->weight *=
                    cel_energy_loss(context, PUMAS_SCHEME_CSDA, material, kf) /
                    cel_energy_loss(context, PUMAS_SCHEME_CSDA, material, ki);
        if (context->decay == PUMAS_DECAY_WEIGHT)
                state->weight *= exp(-fabs(ti - state->time) / s_shared->ctau);

        /* Update the position and direction. */
        if ((locals->magnetized != 0)) {
                enum pumas_return rc;
                rc = transport_csda_deflect(
                    context, state, medium, locals, ki, distance);
                if (rc != PUMAS_RETURN_SUCCESS) return rc;
        } else {
                double path;
                if (context->forward != 0)
                        path = distance;
                else
                        path = -distance;
                state->position[0] += path * state->direction[0];
                state->position[1] += path * state->direction[1];
                state->position[2] += path * state->direction[2];
        }

        /* Register the end of the track, if recording. */
        if (record) record_state(recorder, 0, state);

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Apply the magnetic deflection in CSDA scheme and uniform medium.
 *
 * @param context      The simulation context.
 * @param state        The final state.
 * @param medium       The propagation medium.
 * @param locals       Handle for the local properties of the uniform medium.
 * @param ki           The initial kinetic energy.
 * @param distance     The travelled distance.
 * @return #PUMAS_RETURN_SUCCESS on success, #PUMAS_ERROR otherwise.
 *
 * Compute the magnetic deflection between two kinetic energies using the CSDA.
 * The final kinetic energy is read from the state. At return the final state
 * position and direction are updated.
 */
enum pumas_return transport_csda_deflect(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    struct medium_locals * locals, double ki, double distance)
{
        /* Unpack arguments */
        const double charge = state->charge;
        const double kf = state->kinetic;
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
                        double d = b0_i * b0_i *
                            (magnet[0] * direction[0] +
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
                double pi = a * cel_magnetic_rotation(context, material, ki);
                double pf = a * cel_magnetic_rotation(context, material, kf);
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
                        enum pumas_return rc = PUMAS_RETURN_SUCCESS;
                        rc = csda_magnetic_transport(context, material, density,
                            b0, charge, ki, pi, &xi, &yi, &zi);
                        if (rc != PUMAS_RETURN_SUCCESS) return rc;
                        rc = csda_magnetic_transport(context, material, density,
                            b0, charge, kf, pf, &xf, &yf, &zf);
                        if (rc != PUMAS_RETURN_SUCCESS) return rc;

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
 * @return `PUMAS_RETURN_SUCCESS` on success or `PUMAS_ERROR` otherwise, i.e. if
 * an
 * out of bound error ooccured.
 */
enum pumas_return csda_magnetic_transport(struct pumas_context * context,
    int material, double density, double magnet, double charge, double kinetic,
    double phase, double * x, double * y, double * z)
{
        const int imax = s_shared->n_kinetics - 1;
        double k0 = kinetic;
        int i1, i2;
        if (kinetic <= *table_get_K(0)) {
                k0 = *table_get_K(0);
                i1 = 0;
                i2 = 1;
        } else if (kinetic >= *table_get_K(imax)) {
                /* Neglect deflection at very high energy. */
                *y = 0.0;
                *x = *z = -(cel_grammage(
                                context, PUMAS_SCHEME_CSDA, material, kinetic) -
                    *table_get_X(0, material, imax));
                return PUMAS_RETURN_SUCCESS;
        } else {
                /* Interpolation of the step starting energy. */
                i1 = table_index(context, table_get_K(0), k0);
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

        if (fabs(phase) > max_phi) return PUMAS_RETURN_VALUE_ERROR;

        /* Compute the local x, y and z. */
        double x1, x2, y1, y2, z1, z2;
        {
                /* Negative charge convention. */
                double a = -magnet * charge / density;
                x1 = x2 = y1 = y2 = 0.;
                double u = 1.;
                int j;
                for (j = 0; j < N_LARMOR_ORDERS; j++) {
                        x1 += u * poly_x[j] * (*table_get_Li(material, j, i1));
                        y1 += u * poly_y[j] * (*table_get_Li(material, j, i1));
                        x2 += u * poly_x[j] * (*table_get_Li(material, j, i2));
                        y2 += u * poly_y[j] * (*table_get_Li(material, j, i2));
                        u *= a;
                }
                z1 = *table_get_Li(material, 0, i1);
                z2 = *table_get_Li(material, 0, i2);
        }

        /* Update the position. */
        {
                double h = (k0 - *table_get_K(i1)) /
                    (*table_get_K(i2) - *table_get_K(i1));
                *x = (x1 + h * (x2 - x1)) / density;
                *y = (y1 + h * (y2 - y1)) / density;
                *z = (z1 + h * (z2 - z1)) / density;
        }

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Propagation in arbitrary media.
 *
 * @param context         The simulation context.
 * @param state           The initial/final state.
 * @param medium          The initial propagation medium.
 * @param locals          Handle for the local properties of the starting
 *                        medium.
 * @param step_max_medium The step limitation from the medium.
 * @param step_max_locals The step limitation from the local properties.
 * @return `PUMAS_RETURN_SUCCESS` on success or `PUMAS_ERROR` otherwise.
 *
 * Transport through a set of media described by a medium callback. At output
 * the final kinetic energy, the total distance travelled and the total proper
 * time spent are updated.
 *
 * **Warning** : The initial state must have been initialized before the call.
 */
enum pumas_return transport_with_stepping(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    struct medium_locals * locals, double step_max_medium,
    double step_max_locals)
{
        /* Check the config. */
        if (context->random == NULL) return PUMAS_RETURN_MISSING_RANDOM;

        /* Unpack data. */
        int material = medium->material;

        /* Check for a straight path in a uniform medium of infinite
         * extension.
         */
        const enum pumas_scheme scheme = context->scheme;
        int straight =
            (context->longitudinal && (scheme <= PUMAS_SCHEME_HYBRID) &&
                (step_max_medium <= 0.) && (step_max_locals <= 0.) &&
                !locals->magnetized) ?
            1 :
            0;

        /* Register the start of the the track, if recording. */
        struct pumas_recorder * recorder = context->recorder;
        const int record = (recorder != NULL);
        if (record) record_state(recorder, medium, state);

        /* Initialise some temporary data for the propagation, weights, ect ...
         */
        double ki = state->kinetic;
        double ti = state->time;
        double wi = state->weight;
        double Xi = state->grammage;
        double dei, Xf;
        if (scheme > PUMAS_SCHEME_NO_LOSS) {
                Xf = cel_grammage(context, scheme, material, ki);
                dei = 1. / cel_energy_loss(context, scheme, material, ki);

        } else {
                Xf = dei = 0.;
        }

        /* Check for any initial violation of external limits. */
        if ((context->distance_max > 0.) &&
            (state->distance >= context->distance_max))
                return PUMAS_RETURN_SUCCESS;
        if ((context->grammage_max > 0.) &&
            (state->grammage >= context->grammage_max))
                return PUMAS_RETURN_SUCCESS;
        if ((context->time_max > 0.) && (state->time >= context->time_max))
                return PUMAS_RETURN_SUCCESS;
        if (state->weight <= 0.) return PUMAS_RETURN_SUCCESS;

        /* Initialise the stepping data. */
        struct simulation_context * const context_ =
            (struct simulation_context *)context;
        double grammage_max;
        context_->step_first = 1;
        context_->step_X_limit = (context->kinetic_limit <= 0.) ?
            0. :
            cel_grammage(context, scheme, material, context->kinetic_limit);
        context_->step_invlb1 = 0;
        context_->step_rLarmor = 0.;
        memset(context_->step_uT, 0x0, 3 * sizeof(*(context_->step_uT)));
        transport_limit(context, state, material, Xi, Xf, &grammage_max);
        if (context_->step_event) return PUMAS_RETURN_SUCCESS;

        /* Step through the media. */
        int step_index = 1;
        for (;;) {
                /* Do a transportation step. */
                struct pumas_medium * new_medium = NULL;
                enum pumas_return rc = step_transport(context, state, straight,
                    medium, locals, grammage_max, &step_max_medium,
                    &step_max_locals, &new_medium);
                if (rc != PUMAS_RETURN_SUCCESS) return rc;
                step_index++;

                /* Check for any event. */
                int record_step;
                if ((record != 0) && (recorder->period != 0) &&
                    (step_index != 1)) {
                        record_step = (recorder->period < 0) ?
                            1 :
                            (step_index % recorder->period) == 0;
                } else {
                        record_step = 0;
                }
                if (!context_->step_event && !record_step) continue;

                /* Update the weight if a boundary or hard energy loss
                 * occured.
                 */
                if (!(context_->step_event & EVENT_EHS)) {
                        const double w0 =
                            (context->decay == PUMAS_DECAY_WEIGHT) ?
                            wi * exp(-fabs(state->time - ti) / s_shared->ctau) :
                            wi;
                        state->weight = context->forward ? w0 :
                                                           w0 *
                                cel_energy_loss(
                                    context, scheme, material, state->kinetic) *
                                dei;
                }

                /* Process the event. */
                if (!context_->step_event) {
                        /* Register the current state. */
                        record_state(recorder, medium, state);
                } else {
                        if (context_->step_event &
                            (EVENT_KINETIC | EVENT_RANGE | EVENT_GRAMMAGE |
                                EVENT_TIME | EVENT_WEIGHT)) {
                                /* A boundary was reached. Let's stop the
                                 * simulation.
                                 */
                                break;
                        } else if (context_->step_event &
                            (EVENT_DEL | EVENT_EHS)) {
                                /* A discrete process occured. */
                                if (context_->step_event & EVENT_DEL) {
                                        /* Record the pre step point. */
                                        const int rec =
                                            record && (recorder->period != 0);
                                        if (rec != 0)
                                                record_state(
                                                    recorder, medium, state);

                                        /* Apply the inelastic DEL. */
                                        transport_do_del(
                                            context, state, material);

                                        /* Record the post step point. */
                                        if (rec != 0)
                                                record_state(
                                                    recorder, medium, state);

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
                                            context, state, material);
                                }

                                /* Recompute the geometric step if the
                                 * direction has changed.
                                 */
                                if (!context->longitudinal &&
                                    (step_max_medium > 0.)) {
                                        struct pumas_medium * tmp_medium;
                                        step_max_medium = context->medium(
                                            context, state, &tmp_medium);
                                        if (step_max_medium > 0.)
                                                step_max_medium +=
                                                    0.5 * STEP_MIN;
                                }

                                /* Update the locals if needed. */
                                if (step_max_locals > 0.) {
                                        step_max_locals = transport_set_locals(
                                            medium, state, locals);
                                        if (locals->api.density <= 0.)
                                                return PUMAS_RETURN_DENSITY_ERROR;
                                }
                        } else if (context_->step_event & EVENT_MEDIUM) {
                                /* A medium change occured. Let's update the
                                 * medium.
                                 */
                                medium = new_medium;
                                if (medium == NULL) break;
                                material = medium->material;
                                memset(locals, 0x0, sizeof(*locals));
                                step_max_locals =
                                    transport_set_locals(medium, state, locals);
                                if (locals->api.density <= 0.)
                                        return PUMAS_RETURN_DENSITY_ERROR;
                                straight =
                                    (context->longitudinal &&
                                        (scheme <= PUMAS_SCHEME_HYBRID) &&
                                        (step_max_medium <= 0.) &&
                                        (step_max_locals <= 0.) &&
                                        !locals->magnetized) ?
                                    1 :
                                    0;

                                /* Update the kinetic limit converted
                                 * to grammage for this material.
                                 */
                                context_->step_X_limit =
                                    (context->kinetic_limit <= 0.) ?
                                    0. :
                                    cel_grammage(context, scheme, material,
                                        context->kinetic_limit);

                                /* Reset the stepping data memory. */
                                context_->step_first = 1;
                                context_->step_invlb1 = 0;
                                context_->step_rLarmor = 0.;
                                memset(context_->step_uT, 0x0,
                                    3 * sizeof(*(context_->step_uT)));

                                /* Record the change of medium. */
                                if (record)
                                        record_state(recorder, medium, state);
                        } else {
                                /*  This should not happen. */
                                assert(0);
                        }

                        /* Update the initial conditions and the tracking of
                         * stepping events.
                         */
                        if (!(context_->step_event & EVENT_EHS)) {
                                ki = state->kinetic;
                                ti = state->time;
                                wi = state->weight;
                                Xi = state->grammage;
                                if (scheme > PUMAS_SCHEME_NO_LOSS) {
                                        Xf = cel_grammage(
                                            context, scheme, material, ki);
                                        dei = 1. /
                                            cel_energy_loss(
                                                context, scheme, material, ki);
                                }
                        }
                        transport_limit(
                            context, state, material, Xi, Xf, &grammage_max);
                        if (context_->step_event) break;
                }
        }

        /* Protect final kinetic value against rounding errors. */
        if (context->forward != 0) {
                const double kinetic_min =
                    (context->kinetic_limit < 0.) ? 0. : context->kinetic_limit;
                if (fabs(state->kinetic - kinetic_min) < FLT_EPSILON)
                        state->kinetic = kinetic_min;
        } else {
                if ((context->kinetic_limit > 0.) &&
                    (fabs(state->kinetic - context->kinetic_limit) <
                        FLT_EPSILON))
                        state->kinetic = context->kinetic_limit;
        }

        /* Register the end of the track, if recording. */
        if (record) record_state(recorder, medium, state);

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Set the local properties of a medium.
 *
 * @param state  The Monte-Carlo state.
 * @param locals The local properties.
 * @return The proposed max step length is returned. A null or negative value
 * indicates a uniform and infinite medium.
 *
 * This is an encapsulation of the API locals callback where internal data are
 * also initialised.
 */
double transport_set_locals(struct pumas_medium * medium,
    struct pumas_state * state, struct medium_locals * locals)
{
        struct pumas_locals * loc = (struct pumas_locals *)locals;
        const double step_max = medium->locals(medium, state, loc);
        const double * const b = loc->magnet;
        locals->magnetized =
            ((b[0] != 0.) || (b[1] != 0.) || (b[2] != 0.)) ? 1 : 0;
        return step_max;
}

/**
 * Prepare the various limits for a MC propagation.
 *
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
void transport_limit(struct pumas_context * context,
    const struct pumas_state * state, int material, double Xi, double Xf,
    double * grammage_max)
{
        /* Initialise the stepping event flags. */
        struct simulation_context * const context_ =
            (struct simulation_context *)context;
        context_->step_event = context_->step_foreseen = EVENT_NONE;

        /* Check the Monte-Carlo weight. */
        if (state->weight <= 0.) {
                context_->step_event = EVENT_WEIGHT;
                return;
        }

        /* Check the kinetic limits. */
        if (context->forward != 0) {
                const double kinetic_min =
                    (context->kinetic_limit < 0.) ? 0. : context->kinetic_limit;
                if (state->kinetic <= kinetic_min) {
                        context_->step_event = EVENT_KINETIC;
                        return;
                }
        } else if ((context->kinetic_limit > 0.) &&
            (state->kinetic >= context->kinetic_limit)) {
                context_->step_event = EVENT_KINETIC;
                return;
        };

        /* Initialise with the context grammage limit. */
        if (context->grammage_max <= 0.)
                *grammage_max = 0.;
        else {
                *grammage_max = context->grammage_max;
                context_->step_foreseen = EVENT_GRAMMAGE;
        }

        /* Check the NO LOSS case. */
        const enum pumas_scheme scheme = context->scheme;
        if (scheme == PUMAS_SCHEME_NO_LOSS) {
                if (context->longitudinal == 0) {
                        const double X = Xi -
                            coulomb_ehs_length(
                                context, material, state->kinetic) *
                                log(context->random(context));
                        if ((*grammage_max == 0.) || (X < *grammage_max)) {
                                *grammage_max = X;
                                context_->step_foreseen = EVENT_EHS;
                                return;
                        }
                }
                return;
        }

        /* Check for an inelastic DEL. */
        const double sgn = context->forward ? 1. : -1.;
        enum stepping_event foreseen = EVENT_NONE;
        double kinetic_limit = 0.;
        if (scheme == PUMAS_SCHEME_HYBRID) {
                const double nI =
                    del_interaction_length(context, material, state->kinetic) +
                    sgn * log(context->random(context));
                if (nI > 0.) {
                        const double k = del_kinetic_from_interaction_length(
                            context, material, nI);
                        if (!context->forward ||
                            (k > *table_get_Kt(material))) {
                                kinetic_limit = k;
                                foreseen = EVENT_DEL;
                        }
                }
        }

        /* Check for an EHS event. */
        if ((scheme < PUMAS_SCHEME_DETAILED) && !context->longitudinal) {
                const double nI = ehs_interaction_length(context, scheme,
                                      material, state->kinetic) +
                    sgn * log(context->random(context));
                if (nI > 0.) {
                        const double k = ehs_kinetic_from_interaction_length(
                            context, scheme, material, nI);
                        if ((kinetic_limit <= 0.) ||
                            ((context->forward != 0) && (k > kinetic_limit)) ||
                            ((context->forward == 0) && (k < kinetic_limit))) {
                                kinetic_limit = k;
                                foreseen = EVENT_EHS;
                        }
                }
        }

        /* Return if no discrete event might occur. */
        if (foreseen == EVENT_NONE) return;

        /* Convert the kinetic limit to a grammage one and update. */
        const double X = Xi +
            sgn * (Xf - cel_grammage(context, scheme, material, kinetic_limit));
        if ((*grammage_max <= 0) || (X < *grammage_max)) {
                *grammage_max = X;
                context_->step_foreseen = foreseen;
        }
}

/**
 * Apply an inelastic DEL during the Monte-Carlo propagation.
 *
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param state    The initial/final state.
 */
void transport_do_del(
    struct pumas_context * context, struct pumas_state * state, int material)
{
        /* Update the energy. */
        double ki, kf;
        polar_function_t * polar_func;
        if (context->forward != 0) {
                ki = state->kinetic;
                polar_func = del_randomise_forward(context, state, material);
                kf = state->kinetic;
        } else {
                kf = state->kinetic;
                polar_func = del_randomise_reverse(context, state, material);
                ki = state->kinetic;
        }

        /* Update the direction. */
        if ((context->longitudinal == 0) && (polar_func != NULL)) {
                const double ct = polar_func(context, ki, kf);
                step_rotate_direction(context, state, ct);
        }
}

/*!
 * Apply an Elastic Hard Scattering (EHS) event during the Monte-Carlo
 * propagation.
 *
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param state    The initial/final state.
 */
void transport_do_ehs(
    struct pumas_context * context, struct pumas_state * state, int material)
{
        /* Unpack data. */
        struct simulation_context * ctx = (struct simulation_context *)context;
        struct coulomb_workspace * workspace = ctx->workspace;
        const double kinetic = state->kinetic;

        /* Get the cut-off angle and the soft scattering 1st moment. */
        double mu0 = 0., invlb1 = 0.;
        table_get_msc(context, material, kinetic, &mu0, &invlb1);

        /* Compute the scattering parameters and the the total Coulomb
         * cross section.
         */
        int i;
        struct coulomb_data * data;
        double cs_tot = 0.;
        for (i = 0, data = workspace->data; i < s_shared->elements_in[material];
             i++, data++) {
                /* Compute the scattering parameters. */
                const struct material_component * const component =
                    &s_shared->composition[material][i];
                const struct atomic_element * const element =
                    s_shared->element[component->element];
                double kinetic0;
                coulomb_frame_parameters(
                    kinetic, element->A, &kinetic0, data->fCM);
                data->fspin = coulomb_spin_factor(kinetic0);
                coulomb_screening_parameters(
                    context, kinetic0, component->element, data->screening);
                coulomb_pole_decomposition(data->screening, data->a, data->b);
                data->invlambda = component->fraction /
                    coulomb_wentzel_path(
                        kinetic0, element->Z, element->A, data->screening[0]);

                /* Compute the restricted Coulomb cross-section in the
                 * CM frame.
                 */
                data->cs_hard = data->invlambda *
                    coulomb_restricted_cs(
                        mu0, data->fspin, data->screening, data->a, data->b);
                cs_tot += data->cs_hard;
        }

        /* Randomise the hard scatterer element. */
        const double cc0 = context->random(context) * cs_tot;
        int ihard = s_shared->elements_in[material] - 1;
        double cs = 0.;
        for (i = 0, data = workspace->data; i < s_shared->elements_in[material];
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
        const double f1 = transport_hard_coulomb_objective(mu1, workspace);
        if ((f1 < 0.) && (-f1 / workspace->cs_h >= 1E-06)) {
                /* We hit a form factor suppression. Let's call the
                 * root solver.
                 */
                const double f0 = data->cs_hard * zeta;
                math_find_root(transport_hard_coulomb_objective, mu0, mu1, &f0,
                    &f1, 1E-06 * mu0, 1E-06, 100, workspace, &mu1);
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
 * @param mu         The proposed Coulomb scattering angular parameter.
 * @param parameters A pointer to the scattering workspace.
 * @return The difference between the current and expected restricted cross-
 * section.
 *
 * This is a wrapper for the root solver. It provides the objective function
 * to resolve.
 */
double transport_hard_coulomb_objective(double mu, void * parameters)
{
        /* Unpack the workspace. */
        struct coulomb_workspace * workspace = parameters;
        struct coulomb_data * data = workspace->data + workspace->ihard;

        /* Compute the restricted cross section. */
        double cs_exp = data->invlambda *
            coulomb_restricted_cs(
                mu, data->fspin, data->screening, data->a, data->b);

        /* Return the difference with the expectation. */
        return cs_exp - workspace->cs_h;
}

/**
 * Randomise an inelastic DEL in forward MC.
 *
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param state    The initial/final state.
 * @return The polar function for the randomisation of the corresponding TT or
 * `NULL` if none.
 *
 * Below 10 GeV the DEL is randomised from a power law bias PDF. Above,
 * a ziggurat algorithm is used.
 */

polar_function_t * del_randomise_forward(
    struct pumas_context * context, struct pumas_state * state, int material)
{
        /* Check for a *do nothing* process. */
        if (state->kinetic <= *table_get_Kt(material)) return NULL;

        /* Randomise the target element and the DEL sub-process. */
        struct del_info info;
        del_randomise_target(context, state, material, &info);
        dcs_function_t * const dcs_func = dcs_get(info.process);
        const struct atomic_element * element = s_shared->element[info.element];

        if (info.process == 3) {
                const double m1 = s_shared->mass - ELECTRON_MASS;
                if (state->kinetic <= 0.5 * m1 * m1 / ELECTRON_MASS) {
                        state->kinetic = dcs_ionisation_randomise(
                            context, element, state->kinetic, X_FRACTION);
                        return polar_get(info.process);
                }
        }

        /* Interpolate the lower fractional threshold. */
        int row;
        double xmin, hK;
        if (state->kinetic >= *table_get_K(s_shared->n_kinetics - 1)) {
                row = s_shared->n_kinetics - 1;
                hK = 0.;
                xmin = *table_get_Xt(info.process, info.element, row);
        } else {
                row = table_index(context, table_get_K(0), state->kinetic);
                if (row <= 0) {
                        row = 0;
                        hK = 0.;
                        xmin = *table_get_Xt(info.process, info.element, row);
                } else {
                        const double * const k = table_get_K(row);
                        const double * const y =
                            table_get_Xt(info.process, info.element, row);
                        hK = (state->kinetic - k[0]) / (k[1] - k[0]);
                        xmin = hK * y[1] + (1. - hK) * y[0];
                }
        }

        /* Set the upper fractionnal threshold. */
        double xmax = 1.;
        if (info.process == 3) {
                /* Ionisation case with upper kinematic threshold. */
                const double P2 =
                    state->kinetic * (state->kinetic + 2. * s_shared->mass);
                const double E = state->kinetic + s_shared->mass;
                const double Wmax = 2. * ELECTRON_MASS * P2 /
                    (s_shared->mass * s_shared->mass +
                        ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
                xmax = Wmax / state->kinetic;
                if (xmax > 1.) xmax = 1.;
        }
        if (xmin >= xmax) return NULL;

        if (state->kinetic <= 1.E+01) {
                /* Randomise the final energy with a bias model. */
                const double alpha[N_DEL_PROCESSES] = { 1.5, 3., 1.6, 2.3 };
                double r, w_bias;
                del_randomise_power_law(
                    context, alpha[info.process], xmin, xmax, &r, &w_bias);

                /* Update the kinetic energy and the Monte-Carlo weight. */
                const double d = dcs_evaluate(context, dcs_func, element,
                                     state->kinetic, state->kinetic * (1 - r));
                state->kinetic *= r;
                state->weight *= info.reverse.weight * d * w_bias;
        } else {
                /* Above 10 GeV rely on a direct sampling. */
                float tmp[DCS_SAMPLING_N], *dcs_samples = NULL;
                if ((info.process < N_DEL_PROCESSES - 1) &&
                    (state->kinetic >= DCS_MODEL_MIN_KINETIC)) {
                        /*  Interpolate the tabulated values and pass them
                         * to the sampling routine.
                         */
                        dcs_samples = tmp;
                        const int n = DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q +
                            1 + DCS_SAMPLING_N;
                        const float * y = table_get_dcs_value(element,
                            info.process, row - s_shared->dcs_model_offset);
                        memcpy(tmp, y, DCS_SAMPLING_N * sizeof(float));
                        if (hK > 0.) {
                                y += n;
                                int i;
                                for (i = 0; i < DCS_SAMPLING_N; i++)
                                        tmp[i] = tmp[i] * (float)(1. - hK) +
                                            (float)(hK * y[i]);
                        }
                }
                del_randomise_ziggurat(
                    context, state, dcs_func, element, xmin, xmax, dcs_samples);
        }

        return polar_get(info.process);
}

/**
 * Randomise an inelastic DEL in backward MC.
 *
 * @param context  The simulation context.
 * @param material The index of the propagation material.
 * @param state    The initial/final state.
 * @return The polar function for the randomisation of the corresponding TT or
 * `NULL` if none.
 *
 * A mixture PDF with a weight function is used for the randomisation
 * of the ancestor's state. First we check for a *do nothing* event, then the
 * DEL is processed by randomising the initial kinetic energy over a power law
 * distribution and then randomising the target element.
 */
polar_function_t * del_randomise_reverse(
    struct pumas_context * context, struct pumas_state * state, int material)
{
        /* Check for a pure CEL event. */
        const double lnq0 = -log(X_FRACTION);
        const double kt = *table_get_Kt(material);
        const double kf = state->kinetic;
        const double pCEL =
            (kf >= kt) ? 0. : lnq0 / (lnq0 + log(kt / (kt - kf)));
        if (pCEL && (context->random(context) < pCEL)) {
                state->weight /= pCEL;
                return NULL;
        }
        const double xf =
            (kf >= kt * (1. - X_FRACTION)) ? X_FRACTION : 1. - kf / kt;

        /* Randomise the initial energy with a bias model. */
        double xmax, alpha;
        const double m1 = s_shared->mass - ELECTRON_MASS;
        if (state->kinetic < 0.5 * m1 * m1 / ELECTRON_MASS) {
                const double m2 = s_shared->mass + ELECTRON_MASS;
                xmax = 2. * ELECTRON_MASS * (state->kinetic +
                    2. * s_shared->mass) / (m2 * m2);
                if (xmax < xf) return NULL;
                alpha = RMC_ALPHA_LOW;
        } else {
                alpha = RMC_ALPHA_HIGH;
                xmax = 1.;
        }

        double r, w_bias;
        del_randomise_power_law(context, alpha, xf, xmax, &r, &w_bias);
        w_bias /= r; /* Jacobian factor. */
        state->kinetic /= r;

        /* Randomise the target element and the SEL process. */
        struct del_info info;
        info.reverse.Q = state->kinetic - kf;
        del_randomise_target(context, state, material, &info);
        dcs_function_t * const dcs_func = dcs_get(info.process);
        const struct atomic_element * element = s_shared->element[info.element];

        /* Update the MC weight. */
        const double f = dcs_evaluate(context, dcs_func, element,
                             state->kinetic, state->kinetic - kf) *
            info.reverse.weight;
        state->weight *= w_bias * f *
            del_cross_section(context, material, state->kinetic) /
            (del_cross_section(context, material, kf) * (1. - pCEL));
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
 * @param context The simulation context.
 * @param alpha   The power law exponent.
 * @param xmin    The minimum fractional energy transfer.
 * @param xmax    The maximum fractional energy transfer.
 *
 * Randomise the initial energy using a rejection sampling method. The DCS must
 * be a decreasing function of the energy loss. It is bounded by a piecewise
 * uniform pdf.
 */
void del_randomise_ziggurat(struct pumas_context * context,
    struct pumas_state * state, dcs_function_t * dcs_func,
    const struct atomic_element * element, double xmin, double xmax,
    float * dcs_sampling)
{
        double x_b[DCS_SAMPLING_N + 1], dcs_b[DCS_SAMPLING_N];
        int n_b = DCS_SAMPLING_N, ib;
        if (dcs_sampling == NULL) {
                /* The DCS is uniformly sampled. */
                const double dx =
                    state->kinetic * (xmax - xmin) / DCS_SAMPLING_N;
                x_b[0] = xmin * state->kinetic;
                dcs_b[0] = dcs_evaluate(
                    context, dcs_func, element, state->kinetic, x_b[0]);
                for (ib = 1; ib < DCS_SAMPLING_N; ib++) {
                        x_b[ib] = x_b[ib - 1] + dx;
                        dcs_b[ib] = dcs_evaluate(context, dcs_func, element,
                            state->kinetic, x_b[ib]);
                        if ((dcs_b[ib] <= 0.) || (dcs_b[ib] > dcs_b[ib - 1])) {
                                /* Protect against numeric errors. */
                                x_b[ib] = xmax * state->kinetic;
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
                        x_b[ib] = nu * state->kinetic;
                        dcs_b[ib] = (double)dcs_sampling[ib];
                }
                x_b[DCS_SAMPLING_N] =
                    x_b[DCS_SAMPLING_N - 1] + dnu * state->kinetic;
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
                    context, dcs_func, element, state->kinetic, xsol);
                if (context->random(context) * dcs_b[i] <= dd) break;
        }

        /* Update the kinetic energy. */
        state->kinetic -= xsol;
}

/**
 * Randomise the target element and the sub-process for a DEL.
 *
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
void del_randomise_target(struct pumas_context * context,
    struct pumas_state * state, int material, struct del_info * info)
{
        /* Interpolate the kinetic table. */
        int i1, i2;
        double h;
        i1 = table_index(context, table_get_K(0), state->kinetic);
        if (i1 < 0) {
                i1 = i2 = 0;
                h = 0.;
        } else if (i1 >= s_shared->n_kinetics - 1) {
                i1 = i2 = s_shared->n_kinetics - 1;
                h = 0.;
        } else {
                i2 = i1 + 1;
                const double K1 = *table_get_K(i1);
                h = (state->kinetic - K1) / (*table_get_K(i2) - K1);
        }

        /* Randomise the target element and the DEL process. */
        const struct material_component * component;
        int ic, ic0 = 0, ip;
        for (ic = 0; ic < material; ic++) ic0 += s_shared->elements_in[ic];
        if (context->forward != 0) {
                /* Randomise according to the total cross section. */
                double zeta = context->random(context);
                component = s_shared->composition[material];
                for (ic = ic0; ic < ic0 + s_shared->elements_in[material];
                     ic++, component++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                const double f1 = *table_get_CSf(ip, ic, i1);
                                const double f2 = *table_get_CSf(ip, ic, i2);
                                const double csf = (f2 - f1) * h + f1;
                                if (!(zeta > csf)) {
                                        double csn;
                                        double csn1 = *table_get_CSn(
                                            ip, component->element, i1);
                                        if (csn1 == 0.) {
                                                /* Linear interpolation. */
                                                const double csn2 =
                                                    *table_get_CSn(ip,
                                                        component->element, i2);
                                                csn = (csn2 - csn1) * h + csn1;
                                        } else {
                                                /* Log interpolation. */
                                                csn1 = log(csn1);
                                                const double csn2 =
                                                    log(*table_get_CSn(ip,
                                                        component->element,
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
                        component = s_shared->composition[material];
                        for (ic = ic0;
                             ic < ic0 + s_shared->elements_in[material];
                             ic++, component++) {
                                const struct atomic_element * element =
                                    s_shared->element[component->element];
                                const double d =
                                    dcs_evaluate(context, dcs_get(ip), element,
                                        state->kinetic, info->reverse.Q);
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
                component = s_shared->composition[material];
                for (ic = ic0; ic < ic0 + s_shared->elements_in[material];
                     ic++, component++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                const int iel = component->element;
                                const struct atomic_element * element =
                                    s_shared->element[iel];
                                const double si =
                                    dcs_evaluate(context, dcs_get(ip), element,
                                        state->kinetic, info->reverse.Q) *
                                    component->fraction;
                                s += si;
                                const double csf1 = *table_get_CSf(ip, ic, i1);
                                const double csf2 = *table_get_CSf(ip, ic, i2);
                                const double csf = (csf2 - csf1) * h + csf1;
                                if (!(zeta > s)) {
                                        double csn;
                                        double csn1 =
                                            *table_get_CSn(ip, iel, i1);
                                        if (csn1 == 0.) {
                                                /* Linear interpolation. */
                                                const double csn2 =
                                                    *table_get_CSn(ip, iel, i2);
                                                csn = (csn2 - csn1) * h + csn1;
                                        } else {
                                                /* Log interpolation. */
                                                csn1 = log(csn1);
                                                const double csn2 =
                                                    log(*table_get_CSn(
                                                        ip, iel, i2));
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
 * @param context            The simulation context.
 * @param state              The initial/final state.
 * @param straight           Flag for a straight step.
 * @param medium_index       The index of the propagation medium.
 * @param locals             Handle for the local properties of the medium.
 * @param grammage_max       The maximum grammage until a limit is reached.
 * @param step_max_medium    The stepping limitation from medium boundaries.
 * @param step_max_locals    The stepping limitation from a non uniform medium.
 * @param out_index          The index of the end step medium or a negative
 *                           value if a boundary condition was reached.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise
 * `PUMAS_ERROR`.
 *
 * At return the state kinetic energy, position, direction and distance are
 * updated. In case of a stochastic CEL the proper time is also updated.
 */
enum pumas_return step_transport(struct pumas_context * context,
    struct pumas_state * state, int straight, struct pumas_medium * medium,
    struct medium_locals * locals, double grammage_max,
    double * step_max_medium, double * step_max_locals,
    struct pumas_medium ** out_medium)
{
        /* Unpack the data. */
        struct simulation_context * const context_ =
            (struct simulation_context *)context;
        const enum pumas_scheme scheme = context->scheme;
        double * const direction = state->direction;
        double * const position = state->position;
        const int material = medium->material;
        const double density = locals->api.density;
        const double density_i = 1. / density;
        const double momentum =
            sqrt(state->kinetic * (state->kinetic + 2. * s_shared->mass));

        /* Total grammage for the initial kinetic energy.  */
        const int tmp_scheme =
            (scheme == PUMAS_SCHEME_NO_LOSS) ? PUMAS_SCHEME_CSDA : scheme;
        const double Xtot =
            cel_grammage(context, tmp_scheme, material, state->kinetic);

        /* Compute the local step length. */
        double step_loc, rLarmor = 0., uT[3];
        double invlb1 = 0.;
        if (!straight) {
                /* Compute the kinetic step length. */
                double r = RATIO_ENERGY_LOSS;
                const double k_threshold = 1E+09;
                if ((context->scheme == PUMAS_SCHEME_DETAILED) &&
                    (state->kinetic > k_threshold)) {
                        /* In detailed mode, at very high energies shorter
                         * steps are needed.
                         */
                        double f = k_threshold / state->kinetic;
                        if (f < 0.1) f = 0.1;
                        r *= f;
                }
                step_loc = r * density_i * Xtot;

                if (context->longitudinal == 0) {
                        /* Compute the soft scattering path length. */
                        if (context_->step_first != 0) {
                                double mu0;
                                table_get_msc(context, material, state->kinetic,
                                    &mu0, &invlb1);
                                invlb1 *= density;
                        } else
                                invlb1 = context_->step_invlb1;
                        const double stepT = RATIO_MSC / invlb1;
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
                        const double stepM = RATIO_MAGNETIC * rLarmor;
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
        double step = (*step_max_medium <= 0) ?
            step_loc :
            (*step_max_medium < step_loc ? *step_max_medium : step_loc);
        if ((*step_max_locals > 0.) && (step > *step_max_locals))
                step = *step_max_locals;

        /* Check the total distance limitation. */
        *out_medium = medium;
        enum stepping_event event = EVENT_NONE;
        if (context->distance_max > 0.) {
                const double d = context->distance_max - state->distance;
                if (d <= step) {
                        step = d;
                        event = EVENT_RANGE;
                }
        }

        /* Update the position. */
        const double sgn = context->forward ? 1. : -1.;
        position[0] += step * sgn * direction[0];
        position[1] += step * sgn * direction[1];
        position[2] += step * sgn * direction[2];

        /* Check for a change of medium. */
        struct pumas_medium * end_medium;
        *step_max_medium = context->medium(context, state, &end_medium);
        if (*step_max_medium > 0.) *step_max_medium += 0.5 * STEP_MIN;
        double ds = 0.;
        if (end_medium != medium) {
                /* Check for an exact boundary. */
                double pi[3];
                memcpy(pi, position, sizeof(pi));
                double s1 = 0., s2 = -STEP_MIN;
                position[0] = pi[0] + s2 * sgn * direction[0];
                position[1] = pi[1] + s2 * sgn * direction[1];
                position[2] = pi[2] + s2 * sgn * direction[2];
                struct pumas_medium * tmp_medium;
                const double smax =
                    context->medium(context, state, &tmp_medium);
                if (tmp_medium != medium) {
                        /* Locate the medium change by dichotomy. */
                        *step_max_medium = smax;
                        if (*step_max_medium > 0.)
                                *step_max_medium += 0.5 * STEP_MIN;
                        if (tmp_medium != end_medium) end_medium = tmp_medium;
                        s1 = s2;
                        s2 = -step;
                        while (fabs(s1 - s2) > STEP_MIN) {
                                double s3 = 0.5 * (s1 + s2);
                                position[0] = pi[0] + s3 * sgn * direction[0];
                                position[1] = pi[1] + s3 * sgn * direction[1];
                                position[2] = pi[2] + s3 * sgn * direction[2];
                                const double smax = context->medium(
                                    context, state, &tmp_medium);
                                if (tmp_medium == medium) {
                                        s2 = s3;
                                } else {
                                        *step_max_medium = smax;
                                        if (*step_max_medium > 0.)
                                                *step_max_medium +=
                                                    0.5 * STEP_MIN;
                                        s1 = s3;
                                        /* Update the end medium if required. */
                                        if (tmp_medium != end_medium)
                                                end_medium = tmp_medium;
                                }
                        }
                        position[0] = pi[0] + s2 * sgn * direction[0];
                        position[1] = pi[1] + s2 * sgn * direction[1];
                        position[2] = pi[2] + s2 * sgn * direction[2];
                        step += s1;
                }
                ds = s1 - s2;
                event = EVENT_MEDIUM;
                *out_medium = end_medium;
        }

        /*  Get the end step locals. */
        double Bi[3] = { 0., 0., 0. };
        if (locals->magnetized != 0) {
                Bi[0] = locals->api.magnet[0];
                Bi[1] = locals->api.magnet[1];
                Bi[2] = locals->api.magnet[2];
        }
        if (*step_max_locals > 0.) {
                /* Update the locals. */
                *step_max_locals = transport_set_locals(medium, state, locals);
                if (locals->api.density <= 0.)
                        return PUMAS_RETURN_DENSITY_ERROR;
        }

        /* Offset the end step position for a boundary crossing. */
        if (ds != 0.) {
                position[0] += ds * sgn * direction[0];
                position[1] += ds * sgn * direction[1];
                position[2] += ds * sgn * direction[2];
        }

        /* Set the end step kinetic energy. */
        double k1 = state->kinetic, dk = 0.;
        const double dX = 0.5 * step * (density + locals->api.density);
        if ((scheme >= PUMAS_SCHEME_CSDA) && (scheme <= PUMAS_SCHEME_HYBRID)) {
                /* Deterministic CEL with check for any kinetic limit. */
                const double X = Xtot - sgn * dX;
                if (context->forward && (X <= context_->step_X_limit)) {
                        k1 = context->kinetic_limit;
                        const double grammage =
                            state->grammage + Xtot - context_->step_X_limit;
                        if ((grammage_max <= 0.) || (grammage < grammage_max)) {
                                grammage_max = grammage;
                                event = EVENT_KINETIC;
                        }
                } else if (!context->forward && (context_->step_X_limit > 0.) &&
                    X > context_->step_X_limit) {
                        k1 = context->kinetic_limit;
                        const double grammage =
                            state->grammage + context_->step_X_limit - Xtot;
                        if ((grammage_max <= 0.) || (grammage < grammage_max)) {
                                grammage_max = grammage;
                                event = EVENT_KINETIC;
                        }
                } else
                        k1 = cel_kinetic_energy(context, scheme, material, X);
        } else if (scheme == PUMAS_SCHEME_DETAILED) {
                /* Fluctuate the CEL around its average value. */
                step_fluctuate(context, state, material, Xtot, dX, &k1, &dk);

                /* Check for a kinetic limit. */
                double kinetic_limit = -1.;
                if (context->forward) {
                        const double kinetic_min =
                            (context->kinetic_limit < 0.) ?
                            0. :
                            context->kinetic_limit;
                        if (k1 <= kinetic_min) kinetic_limit = kinetic_min;
                } else {
                        if ((context->kinetic_limit > 0.) &&
                            (k1 >= context->kinetic_limit))
                                kinetic_limit = context->kinetic_limit;
                }
                if (kinetic_limit >= 0.) {
                        const double grammage = state->grammage +
                            fabs(state->kinetic - kinetic_limit) / dk * dX;
                        k1 = kinetic_limit;
                        if ((grammage_max <= 0.) || (grammage < grammage_max)) {
                                grammage_max = grammage;
                                event = EVENT_KINETIC;
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
                double xs_del =
                    del_cross_section(context, material, state->kinetic);
                const double tmp_del = del_cross_section(context, material, k1);
                if (tmp_del > xs_del) xs_del = tmp_del;
                if (xs_del <= 0.) goto no_del_event;
                const double X_del = -log(context->random(context)) / xs_del;
                if (X_del >= Xmax) goto no_del_event;
                const double k_del = state->kinetic - sgn * X_del * dk / dX;
                if (k_del <= 0.) goto no_del_event;
                const double r =
                    del_cross_section(context, material, k_del) / xs_del;
                if (context->random(context) > r) goto no_del_event;
                k_x = k_del;
                Xmax = X_del;
                event = EVENT_DEL;
        no_del_event:

                if (!context->longitudinal) {
                        /* Randomise an EHS. */
                        double kmin, kmax;
                        if (k1 <= state->kinetic)
                                kmin = k1, kmax = state->kinetic;
                        else
                                kmax = k1, kmin = state->kinetic;
                        if (kmin < s_shared->table_K[1]) {
                                const double tmp = kmin;
                                kmin = kmax, kmax = tmp;
                        }
                        double lb_ehs =
                            coulomb_ehs_length(context, material, kmin);
                        if (lb_ehs <= 0.) goto no_ehs_event;
                        const double X_ehs =
                            -log(context->random(context)) * lb_ehs;
                        if (X_ehs >= Xmax) goto no_ehs_event;
                        const double k_ehs =
                            state->kinetic - sgn * X_ehs * dk / dX;
                        if (k_ehs <= 0.) goto no_ehs_event;
                        const double r =
                            coulomb_ehs_length(context, material, k_ehs) /
                            lb_ehs;
                        if ((r <= 0.) || (context->random(context) > 1. / r))
                                goto no_ehs_event;
                        k_x = k_ehs;
                        Xmax = X_ehs;
                        event = EVENT_EHS;
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
                        const double drho = density - locals->api.density;
                        const double tmp =
                            density * density + 2. * dX_ * drho / sf0;
                        h_int = (sqrt(tmp) - density) / drho;
                        step *= h_int;
                }
                state->grammage = grammage_max;
                if (!(event & EVENT_KINETIC)) {
                        /*  Update the kinetic energy. */
                        if (scheme <= PUMAS_SCHEME_HYBRID) {
                                if (scheme != PUMAS_SCHEME_NO_LOSS)
                                        k1 = cel_kinetic_energy(context, scheme,
                                            material,
                                            Xtot -
                                                sgn * (state->grammage - Xi));
                                event = context_->step_foreseen;
                        } else if (!(event & (EVENT_EHS | EVENT_DEL))) {
                                /*
                                 * This is a grammage limit in the detailed
                                 * scheme.
                                 */
                                k1 = state->kinetic -
                                    sgn * (state->grammage - Xi) * dk / dX;
                                if (k1 < 0.) k1 = 0.;
                                event = EVENT_GRAMMAGE;
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
        double time_max = context->time_max;
        if (context->decay == PUMAS_DECAY_PROCESS) {
                if ((time_max <= 0.) || (context_->lifetime < time_max)) {
                        time_max = context_->lifetime;
                        decayed = 1;
                }
        }

        const double sf1 = step;
        if (straight && (scheme != PUMAS_SCHEME_NO_LOSS)) {
                const double Ti =
                    cel_proper_time(context, scheme, material, state->kinetic);
                if (time_max > 0.) {
                        const double Tf =
                            Ti - sgn * (time_max - state->time) * density;
                        if (Tf > 0.) {
                                const double dxT = fabs(Xtot -
                                    cel_grammage_as_time(
                                        context, scheme, material, Tf));
                                if (Xi + dxT < state->grammage) {
                                        /* A proper time limit is reached. */
                                        event = EVENT_TIME;
                                        state->time = time_max;
                                        state->decayed = decayed;
                                        step = dxT * density_i;
                                        if (step > sf1) step = sf1;
                                        state->grammage = Xi + dxT;
                                        const double xf = Xtot - sgn * dxT;
                                        k1 = cel_kinetic_energy(
                                            context, scheme, material, xf);
                                }
                        }
                }

                if (event != EVENT_TIME) {
                        const double Tf =
                            cel_proper_time(context, scheme, material, k1);
                        state->time += fabs(Tf - Ti) * density_i;
                }
        } else {
                const double p_f = (k1 <= 0) ?
                    momentum :
                    sqrt(k1 * (k1 + 2. * s_shared->mass));
                const double ti = state->time;
                state->time +=
                    0.5 * step * s_shared->mass * (1. / momentum + 1. / p_f);
                if ((time_max > 0.) && (state->time >= time_max)) {
                        /* A proper time limit is reached. */
                        event = EVENT_TIME;
                        state->time = time_max;
                        state->decayed = decayed;

                        /* Interpolate the step length. */
                        const double a = (momentum / p_f - 1.) / sf1;
                        const double c =
                            -2. * (time_max - ti) * momentum / s_shared->mass;
                        if (a != 0.) {
                                double delta = 1. - a * c;
                                if (delta < 0.) delta = 0.;
                                step = 1. / a * (sqrt(delta) - 1.);
                        } else {
                                step = -0.5 * c;
                        }
                        if (step < 0.) step = 0.;

                        /*  Update the kinetic energy. */
                        if (scheme != PUMAS_SCHEME_NO_LOSS) {
                                const double p1 = momentum / (1 + a * step);
                                k1 = sqrt(p1 * p1 -
                                         s_shared->mass * s_shared->mass) -
                                    s_shared->mass;
                                if (k1 < 0.) k1 = 0.;
                        }

                        /* Correct the grammage. */
                        const double Xf = Xi + dX;
                        if (density == locals->api.density)
                                state->grammage = Xi + step * density;
                        else
                                state->grammage = Xi +
                                    step *
                                        (density +
                                            0.5 * step / sf0 *
                                                (locals->api.density -
                                                    density));
                        if (state->grammage > Xf) state->grammage = Xf;
                }
        }

        if (event == EVENT_TIME) {
                /* Correct the position. */
                const double ds_ = step - sf1;
                position[0] += ds_ * sgn * direction[0];
                position[1] += ds_ * sgn * direction[1];
                position[2] += ds_ * sgn * direction[2];
        }

        /* Update the event flag. */
        context_->step_event = event;

        /* Update the kinetic energy. */
        state->kinetic = k1;

        /* Update the travelled distance. */
        state->distance += step;

        /* Compute the multiple scattering path length. */
        if (context->longitudinal == 0) {
                double mu0, invlb1_;
                if (state->kinetic <= 0.) {
                        context_->step_invlb1 = invlb1;
                } else {
                        table_get_msc(
                            context, material, state->kinetic, &mu0, &invlb1_);
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
                                const double p = (state->kinetic <= 0) ?
                                    momentum :
                                    sqrt(state->kinetic *
                                        (state->kinetic + 2. * s_shared->mass));
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
        int rotated = 0;
        if ((rLarmor > 0.) || (context_->step_rLarmor > 0.)) {
                const double u[3] = { direction[0], direction[1],
                        direction[2] };
                double theta_i =
                    (rLarmor > 0.) ? state->charge * step / rLarmor : 0.;
                double theta_f = (context_->step_rLarmor > 0.) ?
                    state->charge * step / context_->step_rLarmor :
                    0.;
                if (context->forward == 0) {
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
                rotated = 1;
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
                rotated = 1;
        }

        /* Update the geometric step length if the particle was rotated. */
        if (rotated && (*step_max_medium > 0.)) {
                struct pumas_medium * tmp_medium;
                *step_max_medium = context->medium(context, state, &tmp_medium);
                if (*step_max_medium > 0.) *step_max_medium += 0.5 * STEP_MIN;
        }

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Apply a stochastic CEL on a MC step.
 *
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
static void step_fluctuate(struct pumas_context * context,
    struct pumas_state * state, int material, double Xtot, double dX,
    double * kf, double * dE)
{
        const enum pumas_scheme scheme = context->scheme;
        const double sgn = context->forward ? 1. : -1.;
        double k1, dk = 0.;
        k1 = cel_kinetic_energy(context, scheme, material, Xtot - sgn * dX);
        if (k1 > 0.) {
                const double dk1 = sqrt(0.5 * dX *
                    (step_fluctuations2(material, state->kinetic) +
                        step_fluctuations2(material, k1)));
                const double dk0 = fabs(state->kinetic - k1);
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
                                k1 = state->kinetic - sgn * b * u;
                        } else {
                                k1 = state->kinetic;
                        }
                }
                dk = fabs(state->kinetic - k1);
        }
        if (dk == 0.)
                dk =
                    cel_energy_loss(context, scheme, material, state->kinetic) *
                    dX;

        /* Copy back the result. */
        *kf = k1;
        *dE = dk;
}

/**
 * Squared standard deviation of the stochastic CEL.
 *
 * @param material The target material.
 * @param kinetic  The projectile kinetic energy.
 * @return The squared standard deviation.
 *
 * The Gaussian thick aborber approximation of ICRU Report 49 is assumed.
 */
static double step_fluctuations2(int material, double kinetic)
{
        const double r = ELECTRON_MASS / s_shared->mass;
        const double g = 1. + kinetic / s_shared->mass;
        const double b2 = 1. - 1. / (g * g);
        double qmax = 2. * ELECTRON_MASS * b2 * g * g / (1. + r * (2. * g + r));
        const double qcut = X_FRACTION * kinetic;
        if (qmax > qcut) qmax = qcut;
        return 1.535375E-05 * s_shared->material_ZoA[material] *
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
                        const double nrm = 1. /
                            sqrt(direction[0] * direction[0] +
                                direction[2] * direction[2]);
                        u0x = -direction[2] * nrm, u0z = direction[0] * nrm;
                } else {
                        const double nrm = 1. /
                            sqrt(direction[1] * direction[1] +
                                direction[2] * direction[2]);
                        u0y = direction[2] * nrm, u0z = -direction[1] * nrm;
                }
        } else {
                if (a1 > a2) {
                        const double nrm = 1. /
                            sqrt(direction[0] * direction[0] +
                                direction[1] * direction[1]);
                        u0x = direction[1] * nrm, u0y = -direction[0] * nrm;
                } else {
                        const double nrm = 1. /
                            sqrt(direction[1] * direction[1] +
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
 * @param name The element name.
 * @return The element index or `-1` if it wasn't found.
 */
int element_index(const char * name)
{
        int index = 0;
        for (; index < s_shared->n_elements; index++)
                if (strcmp(s_shared->element[index]->name, name) == 0)
                        return index;
        return -1;
}

/*
 * Low level routines: various functions related to Coulomb scattering.
 */
/**
 * Compute the atomic and nuclear screening parameters.
 *
 * @param context   The simulation context.
 * @param kinetic   The kinetic energy.
 * @param element   The scatterer atomic element.
 * @param screening The computed screening parameters.
 *
 * The atomic screening parameter is computed according to Kuraev's
 * parameterisation for relativistic particles or using Moliere's one at low
 * energies.
 */
void coulomb_screening_parameters(struct pumas_context * context,
    double kinetic, int element, double * screening)
{
        /* Nuclear screening. */
        const double third = 1. / 3;
        const double A13 = pow(s_shared->element[element]->A, third);
        const double R1 = 1.02934 * A13 + 0.435;
        const double R2 = 2.;
        const double p2 = kinetic * (kinetic + 2. * s_shared->mass);
        const double d = 5.8406E-02 / p2;
        screening[1] = d / (R1 * R1);
        screening[2] = d / (R2 * R2);

        /* Atomic Moliere screening with Coulomb correction from Kuraev et al.
         * Phys. Rev. D 89, 116016 (2014). Valid for ultra-relativistic
         * particles only.
         */
        struct atomic_element * e = s_shared->element[element];
        const double etot = kinetic + s_shared->mass;
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
 * @param kinetic   The kinetic energy.
 * @param Z         The atomic number of the target element.
 * @param A         The atomic mass of the target element.
 * @param screening The atomic screening factor.
 * @return The mean free path in kg/m^2.
 */
double coulomb_wentzel_path(
    double kinetic, double Z, double A, double screening)
{
        const double d = kinetic * (kinetic + 2. * s_shared->mass) /
            (Z * (kinetic + s_shared->mass));
        return A * 2.54910918E+08 * screening * (1. + screening) * d * d;
}

/**
 * Encapsulation of the interaction length for EHS.
 *
 * @param context  The simulation context.
 * @param material The propagation material.
 * @param kinetic  The kinetic energy.
 * @return The grammage length for hard scattering events.
 */
double coulomb_ehs_length(
    struct pumas_context * context, int material, double kinetic)
{
        const double p2 = kinetic * (kinetic + 2. * s_shared->mass);
        const int imax = s_shared->n_kinetics - 1;
        if (kinetic < *table_get_K(1)) {
                return *table_get_Lb(material, 1) / p2;
        } else if (kinetic >= *table_get_K(imax)) {
                return *table_get_Lb(material, imax) / p2;
        } else {
                const int i1 = table_index(context, table_get_K(0), kinetic);
                const int i2 = i1 + 1;
                double h = (kinetic - *table_get_K(i1)) /
                    (*table_get_K(i2) - *table_get_K(i1));
                return (*table_get_Lb(material, i1) +
                           h *
                               (*table_get_Lb(material, i2) -
                                   *table_get_Lb(material, i1))) /
                    p2;
        }
}

/**
 * Compute the spin factor to the Coulomb DCS.
 *
 * @param kinetic The kinetic energy.
 * @return The spin factor.
 */
double coulomb_spin_factor(double kinetic)
{
        const double e = kinetic + s_shared->mass;
        return kinetic * (e + s_shared->mass) / (e * e);
}

/**
 * Compute the parameters of the CM to the Lab frame transform for the
 * Coulomb DCS.
 *
 * @param kinetic    The kinetic energy.
 * @param Ma         The target mass, in atomic unit.
 * @param kinetic0   The CM kinetic energy.
 * @param parameters The vector of Lorentz parameters (gamma, tau).
 */
void coulomb_frame_parameters(
    double kinetic, double Ma, double * kinetic0, double * parameters)
{
        Ma *= 0.931494; /* Convert from atomic unit to GeV. */
        double M2 = s_shared->mass + Ma;
        M2 *= M2;
        const double sCM12i = 1. / sqrt(M2 + 2. * Ma * kinetic);
        parameters[0] = (kinetic + s_shared->mass + Ma) * sCM12i;
        *kinetic0 =
            (kinetic * Ma + s_shared->mass * (s_shared->mass + Ma)) * sCM12i -
            s_shared->mass;
        if (*kinetic0 < 1E-09) *kinetic0 = 1E-09;
        const double etot = kinetic + s_shared->mass + Ma;
        const double betaCM2 =
            kinetic * (kinetic + 2. * s_shared->mass) / (etot * etot);
        double rM2 = s_shared->mass / Ma;
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
 * @param element The target atomic element.
 * @param kinetic The projectile initiale kinetic energy.
 * @return The inverse of the electronic 1st transport path length in kg/m^2.
 *
 * The contribution from atomic electronic shells is computed following
 * Salvat et al., NIMB316 (2013) 144-159, considering only close interactions
 * and approximating the electronic structure by a single shell of energy I
 * with occupancy Z.
 */
double transverse_transport_ionisation(
    const struct atomic_element * element, double kinetic)
{
        /* Soft close interactions, restricted to X_FRACTION. */
        const double momentum2 = kinetic * (kinetic + 2. * s_shared->mass);
        const double E = kinetic + s_shared->mass;
        const double Wmax = 2. * ELECTRON_MASS * momentum2 /
            (s_shared->mass * s_shared->mass +
                ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
        const double W0 = 2. * momentum2 / ELECTRON_MASS;
        const double mu_max = Wmax / W0;
        double mu3 = kinetic * X_FRACTION / W0;
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
double transverse_transport_photonuclear(
    const struct atomic_element * element, double kinetic)
{
        /* Integration over the kinetic transfer, q, done with a log sampling.
         */
        const double E = kinetic + s_shared->mass;
        double x0 = log(1E-06);
        double x1 = 0.;
        math_gauss_quad(100, &x0, &x1); /* Initialisation. */

        double xi, wi, lbipn = 0.;
        while (math_gauss_quad(0, &xi, &wi) == 0) { /* Stepping. */
                const double nu = X_FRACTION * exp(xi);
                const double q = nu * kinetic;

                /* Analytical integration over mu. */
                const double m02 = 0.4;
                const double q2 = q * q;
                const double tmax = 1.876544 * q;
                const double tmin =
                    q2 * s_shared->mass * s_shared->mass / (E * (E - q));
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
                const double ratio = (I1 * tmax - I2) /
                    ((I0 * tmax - I1) * kinetic *
                        (kinetic + 2. * s_shared->mass));

                /* Update the double integral value.  */
                lbipn += ratio * nu *
                    dcs_photonuclear(element, kinetic, nu * kinetic) * wi;
        }
        return 2. * lbipn;
}

/* Low level routine: helper function for recording a MC state. */
/**
 * Register a Monte-Carlo state.
 *
 * @param recorder     The recorder handle.
 * @param medium       The medium in which the particle is located.
 * @param state        The Monte-Carlo state to record.
 *
 * This routine adds the given state to the recorder's stack.
 */
void record_state(struct pumas_recorder * recorder,
    struct pumas_medium * medium, struct pumas_state * state)
{
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
 * Utility function for encapsulating `returns` with an error handler.
 *
 * @param context    The simulation context or `NULL`.
 * @param rc         The return code.
 * @param caller     The calling function from which to return.
 * @return The return code is forwarded.
 *
 * Note that depending on the value of *context* the function operates with
 * the library global error handler or the context's one.
 */
static enum pumas_return pumas_return(
    enum pumas_return rc, pumas_function_t * caller)
{
        if (rc == PUMAS_RETURN_SUCCESS) return rc;

        if (s_error.catch) {
                if (s_error.catch_rc == PUMAS_RETURN_SUCCESS) {
                        s_error.catch_rc = rc;
                        s_error.catch_function = caller;
                }
                return rc;
        }

        if (s_error.handler != NULL) {
                struct pumas_error error, *ptr;
                if (caller == (pumas_function_t *)pumas_initialise) {
                        ptr = &error;
                        error.file = s_error.file;
                        error.line = s_error.line;
                } else
                        ptr = NULL;
                s_error.handler(rc, caller, ptr);
        }

        return rc;
}

/*
 * Low level routines: I/O and parsing.
 */
/**
 * Parse a dE/dX file.
 *
 * @param fid         The file handle.
 * @param material    The material index.
 * @param line        The line where any error occured.
 * @return On succces `PUMAS_RETURN_SUCCESS`, otherwise `PUMAS_ERROR`.
 *
 * Parse a dE/dX data table in PDG text file format.
 */
enum pumas_return io_parse_dedx_file(FILE * fid, int material, int * line)
{
        enum pumas_return rc = PUMAS_RETURN_SUCCESS;
        char * buffer = NULL;

        /* Skip the header lines. */
        *line = 0;
        int i;
        for (i = 0; i < s_shared->n_energy_loss_header; i++) {
                rc = io_read_line(fid, &buffer);
                (*line)++;
                if (rc != PUMAS_RETURN_SUCCESS)
                        return PUMAS_RETURN_FORMAT_ERROR;
        }

        /* Initialise the new table. */
        int row = 0;
        *table_get_T(PUMAS_SCHEME_CSDA, material, row) = 0.;
        *table_get_T(PUMAS_SCHEME_HYBRID, material, row) = 0.;
        *table_get_K(row) = 0.;
        *table_get_dE(PUMAS_SCHEME_CSDA, material, row) = 0.;
        *table_get_dE(PUMAS_SCHEME_HYBRID, material, row) = 0.;
        *table_get_NI_el(PUMAS_SCHEME_CSDA, material, row) = 0.;
        *table_get_NI_el(PUMAS_SCHEME_HYBRID, material, row) = 0.;
        *table_get_NI_in(material, row) = 0.;
        *table_get_CS(material, row) = 0.;
        int ip;
        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                *table_get_CSf(ip, material, row) = 0.;
                *table_get_CSn(ip, material, row) = 0.;
        }
        row++;

        /* Scan the new table. */
        while (rc == PUMAS_RETURN_SUCCESS) {
                rc = io_read_line(fid, &buffer);
                (*line)++;
                if (rc != PUMAS_RETURN_SUCCESS) break;
                rc = io_parse_dedx_row(buffer, material, &row);
        }

        if (rc != PUMAS_RETURN_SUCCESS) {
                if ((rc == PUMAS_RETURN_END_OF_FILE) || feof(fid))
                        rc = PUMAS_RETURN_SUCCESS;
        }
        if (rc != PUMAS_RETURN_SUCCESS) return rc;

        if (row != s_shared->n_kinetics) return PUMAS_RETURN_FORMAT_ERROR;

        compute_regularise_del(material);

        return rc;
}

/**
 * Parse a row of a dE/dX file.
 *
 * @param buffer   The read buffer containing the row data.
 * @param material The index of the material.
 * @param row      The index of the parsed row.
 * @return On succees `PUMAS_RETURN_SUCCESS`, or `PUMAS_ERROR` otherwise.
 *
 * Parse a row of a dE/dX data table formated in PDG text file format.
 */
enum pumas_return io_parse_dedx_row(char * buffer, int material, int * row)
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
        if ((count != 7) || (*row >= s_shared->n_kinetics))
                return PUMAS_RETURN_FORMAT_ERROR;

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
                *table_get_K(*row) = k;
        } else if (*table_get_K(*row) != k)
                return PUMAS_RETURN_VALUE_ERROR;

        /* Compute the fractional contributions of the energy loss
         * processes to the CEL and to DELs.
         */
        static double * cel_table = NULL;
        if (material == 0) {
                /* Precompute the per element terms. */
                if ((cel_table = compute_cel_and_del(*row)) == NULL)
                        return PUMAS_RETURN_MEMORY_ERROR;
        }

        struct material_component * component = s_shared->composition[material];
        double frct_cel[] = { 0., 0., 0., 0. };
        double frct_cs[] = { 0., 0., 0., 0. };
        int ic, ic0 = 0;
        for (ic = 0; ic < material; ic++) ic0 += s_shared->elements_in[ic];
        for (ic = ic0; ic < ic0 + s_shared->elements_in[material];
             ic++, component++) {
                int iel = component->element;
                const double w = component->fraction;
                /* Loop over processes. */
                int ip;
                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                        const double f = (*table_get_CSn(ip, iel, *row)) * w;
                        *table_get_CSf(ip, ic, *row) = f;
                        frct_cs[ip] += f;
                        frct_cel[ip] +=
                            *table_get_cel(ip, iel, *row, cel_table) * w;
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
                for (ic = ic0; ic < ic0 + s_shared->elements_in[material]; ic++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++)
                                *table_get_CSf(ip, ic, *row) = 0.;
        } else {
                double sum_tot = 0.;
                for (ic = ic0; ic < ic0 + s_shared->elements_in[material]; ic++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++)
                                sum_tot += *table_get_CSf(ip, ic, *row);
                double sum = 0.;
                for (ic = ic0; ic < ic0 + s_shared->elements_in[material]; ic++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                sum += *table_get_CSf(ip, ic, *row);
                                *table_get_CSf(ip, ic, *row) = sum / sum_tot;
                        }

                /* Protect against rounding errors. */
                ic = ic0 + s_shared->elements_in[material] - 1;
                ip = N_DEL_PROCESSES - 1;
                *table_get_CSf(ip, ic, *row) = 1.;
        }

        /* Update the table values */
        const double de_cel = de - be_cel;

        /* If this is the last entry, save the energy loss values. */
        if (*row == s_shared->n_kinetics - 1) {
                const double etot = k + s_shared->mass;
                *table_get_a_max(material) = a;
                *table_get_b_max(PUMAS_SCHEME_CSDA, material) = be / etot;
                *table_get_b_max(PUMAS_SCHEME_HYBRID, material) =
                    (be - be_cel) / etot;
        }

        /* End point statistics */
        *table_get_dE(PUMAS_SCHEME_CSDA, material, *row) = de;
        *table_get_dE(PUMAS_SCHEME_HYBRID, material, *row) = de_cel;
        *table_get_CS(material, *row) = frct_cs_del;

        /* Weighted integrands */
        const double dei = 1. / de_cel;
        *table_get_X(PUMAS_SCHEME_HYBRID, material, *row) = dei;
        *table_get_NI_in(material, *row) = frct_cs_del * dei;

        (*row)++;
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Read a line from a file.
 *
 * @param fid The file handle.
 * @param buf A pointer to the new line or `NULL` in case of faillure.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise an error code.
 *
 * This routine manages a dynamic buffer. If *fid* is not `NULL`, this routine
 * reads a new line and fills in a pointer to a string holding the read data.
 * Otherwise, if *fid* is `NULL` any memory previously allocated by the routine
 * is released.
 */
enum pumas_return io_read_line(FILE * fid, char ** buf)
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
                if (buffer == NULL) return PUMAS_RETURN_MEMORY_ERROR;
        }

        /* Get a full line. */
        char * s = buffer;
        int to_read = size;
        for (;;) {
                const char endline1 = '\n';
                const char endline2 = '\r';
                const char * r = fgets(s, to_read, fid);
                if (r == NULL) return PUMAS_RETURN_END_OF_FILE;
                int n_read = strlen(s);
                if ((n_read >= to_read - 1) && (s[to_read - 2] != endline1) &&
                    (s[to_read - 2] != endline2)) {
                        size += 2048;
                        char * new_buffer =
                            reallocate(buffer, size * sizeof(*buffer));
                        if (new_buffer == NULL)
                                return PUMAS_RETURN_MEMORY_ERROR;
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
 * @param mdf          The MDF handle.
 * @param dedx_path    The path to the energy loss table(s).
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise an error code.
 */
enum pumas_return mdf_parse_settings(
    struct mdf_buffer * mdf, const char * dedx_path)
{
        /* Initialisation of settings. */
        mdf->n_kinetics = 0;
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
        enum pumas_return rc;
        char * full_path = NULL;
        char * filename = NULL;
        int offset_dir, size_path, size_name = 0;
        if (!mdf->dry_mode) {
                mdf->size_dedx_path = strlen(dedx_path) + 1;
                if (mdf_format_path(dedx_path, mdf->mdf_path, &full_path,
                        &offset_dir, &size_path) != PUMAS_RETURN_SUCCESS) {
                        rc = PUMAS_RETURN_MEMORY_ERROR;
                        goto clean_and_exit;
                }
        }

        /* Prepare the index tables. */
        if (mdf_settings_index(MDF_INDEX_INITIALISE, 0) !=
            PUMAS_RETURN_SUCCESS) {
                rc = PUMAS_RETURN_MEMORY_ERROR;
                goto clean_and_exit;
        }

        /* Loop on the XML nodes. */
        const int pad_size = sizeof(*(s_shared->data));
        mdf->left = 0;
        mdf->line = 1;
        mdf->depth = MDF_DEPTH_EXTERN;
        struct mdf_node node;
        for (;;) {
                /* Get the next node. */
                rc = mdf_get_node(mdf, &node);
                /* Check the termination. */
                if (rc == PUMAS_RETURN_END_OF_FILE) {
                        if (mdf->depth == MDF_DEPTH_EXTERN) {
                                /* This is a normal terminations. */
                                rc = PUMAS_RETURN_SUCCESS;
                                break;
                        } else {
                                goto clean_and_exit;
                        }
                } else if (rc != PUMAS_RETURN_SUCCESS)
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
                                                rc = PUMAS_RETURN_MEMORY_ERROR;
                                                goto clean_and_exit;
                                        }
                                        strcpy(filename, node.at2.file);
                                }

                                /* Update the names table and book a new index
                                 * table.
                                 */
                                if ((mdf_settings_name(
                                         size, 'M', node.at1.name) !=
                                        PUMAS_RETURN_SUCCESS) ||
                                    (mdf_settings_index(
                                         MDF_INDEX_WRITE_MATERIAL, 0) !=
                                        PUMAS_RETURN_SUCCESS)) {
                                        rc = PUMAS_RETURN_MEMORY_ERROR;
                                        goto clean_and_exit;
                                }
                        } else {
                                if (mdf->elements_in > mdf->max_components)
                                        mdf->max_components = mdf->elements_in;
                                /* Finalise the index table. */
                                mdf_settings_index(MDF_INDEX_FINALISE_MATERIAL,
                                    mdf->elements_in);
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
                                        mdf->n_elements) !=
                                    PUMAS_RETURN_SUCCESS) {
                                        rc = PUMAS_RETURN_MEMORY_ERROR;
                                        goto clean_and_exit;
                                }
                        } else {
                                /* Analyse the index table. */
                                int elements_in = mdf_settings_index(
                                    MDF_INDEX_FINALISE_COMPOSITE, 0);
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
                        if (mdf_settings_name(size, 'E', node.at1.name) !=
                            PUMAS_RETURN_SUCCESS) {
                                rc = PUMAS_RETURN_MEMORY_ERROR;
                                goto clean_and_exit;
                        }
                } else if (node.key == MDF_KEY_ATOMIC_COMPONENT) {
                        mdf->n_components++;
                        int index = mdf_settings_name(0, 'E', node.at1.name);
                        if (index < 0) {
                                rc = PUMAS_RETURN_UNKNOWN_ELEMENT;
                                goto clean_and_exit;
                        }
                        if (mdf_settings_index(MDF_INDEX_WRITE_MATERIAL,
                                index) != PUMAS_RETURN_SUCCESS) {
                                rc = PUMAS_RETURN_FORMAT_ERROR;
                                goto clean_and_exit;
                        }
                } else if (node.key == MDF_KEY_COMPOSITE_COMPONENT) {
                        /* Update the components count. */
                        mdf->materials_in++;

                        /* Update the index table. */
                        int index = mdf_settings_name(0, 'M', node.at1.name);
                        if (index < 0) {
                                rc = PUMAS_RETURN_UNKNOWN_MATERIAL;
                                goto clean_and_exit;
                        }
                        mdf_settings_index(MDF_INDEX_UPDATE_COMPOSITE, index);
                }
        }
        mdf->line = 0;

        /* Check the content .*/
        if ((mdf->n_elements == 0) || (mdf->n_materials == 0)) {
                /* There are no elements or materials. */
                rc = PUMAS_RETURN_INCOMPLETE_FILE;
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
                                rc = PUMAS_RETURN_MEMORY_ERROR;
                                goto clean_and_exit;
                        }
                        full_path = new_name;
                }
                strcpy(full_path + offset_dir, filename);

                rc = mdf_parse_kinetic(mdf, full_path);
        } else
                mdf->n_kinetics = 1;

clean_and_exit:
        /* Copy any info on the error. */
        if ((rc != PUMAS_RETURN_SUCCESS) && (mdf->line > 0)) {
                s_error.line = mdf->line;
                if (s_error.file[0] == '\0')
                        strncpy(s_error.file, mdf->mdf_path, ERROR_FILE_LENGTH);
                s_error.file[ERROR_FILE_LENGTH - 1] = '\0';
        }

        /* Free the temporary memory and return. */
        mdf_settings_name(-1, 0x0, NULL);
        mdf_settings_index(MDF_INDEX_FREE, 0);
        deallocate(filename);
        deallocate(full_path);
        return rc;
}

/**
 * Manage a temporary mapping of material indices.
 *
 * @param operation The operation to perform.
 * @param value     An optional value for the operation.
 * @return The return value depends on the operation.
 */
int mdf_settings_index(int operation, int value)
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
                for (i = 0; i < value; i++) { table += (*table) + 1; }
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
                        return PUMAS_RETURN_MEMORY_ERROR;
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

        return PUMAS_RETURN_INDEX_ERROR;
}

/**
 * Manage a temporary table of material names.
 *
 * @param size   The operation to perform or the size of the new name.
 * @param prefix A prefix for the category of the name (base, composite, ...).
 * @param name   The name to handle.
 * @return The return value depends on the operation performed.
 */
int mdf_settings_name(int size, char prefix, const char * name)
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
                        return PUMAS_RETURN_MEMORY_ERROR;
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
 * @param mdf  The MDF handle.
 * @param path The full path to the dE/dX file.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
 */
enum pumas_return mdf_parse_kinetic(struct mdf_buffer * mdf, const char * path)
{
        /* Initialise the settings. */
        mdf->n_kinetics = 0;
        mdf->dcs_model_offset = 0;
        mdf->line = 0;

        /* Open the dedx file. */
        FILE * fid = fopen(path, "r");
        if (fid == NULL) {
                mdf->line = 0;
                strncpy(s_error.file, path, ERROR_FILE_LENGTH);
                return PUMAS_RETURN_PATH_ERROR;
        }

        /* Skip the header lines. */
        char * buffer = NULL;
        enum pumas_return rc;
        for (mdf->n_energy_loss_header = 0;; mdf->n_energy_loss_header++) {
                rc = io_read_line(fid, &buffer);
                mdf->line++;
                if (rc != PUMAS_RETURN_SUCCESS) {
                        fclose(fid);
                        return rc;
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
                rc = io_read_line(fid, &buffer);
                mdf->line++;
                if (rc != PUMAS_RETURN_SUCCESS) {
                        if ((rc == PUMAS_RETURN_END_OF_FILE) || feof(fid)) {
                                rc = PUMAS_RETURN_SUCCESS;
                                break;
                        }
                }
        }
        io_read_line(NULL, NULL);

        /*  Update the settings. */
        if (rc == PUMAS_RETURN_SUCCESS)
                mdf->line = 0;
        else
                strncpy(s_error.file, path, ERROR_FILE_LENGTH);
        mdf->n_kinetics = nk;
        mdf->dcs_model_offset = offset;

        fclose(fid);
        return rc;
}

/**
 * Parse the atomic elements from a MDF.
 *
 * @param mdf The MDF handle.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
 */
enum pumas_return mdf_parse_elements(struct mdf_buffer * mdf)
{
        /* Loop on the XML nodes. */
        const int m = s_shared->n_kinetics - s_shared->dcs_model_offset;
        const int n =
            DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1 + DCS_SAMPLING_N;
        const int pad_size = sizeof(*(s_shared->data));
        const int c_mem = memory_padded_size(
            m * n * (N_DEL_PROCESSES - 1) * sizeof(float), pad_size);
        rewind(mdf->fid);
        mdf->left = 0;
        mdf->line = 1;
        mdf->depth = MDF_DEPTH_EXTERN;
        struct mdf_node node;
        enum pumas_return rc;
        int iel = 0;
        for (;;) {
                /* Get the next node. */
                rc = mdf_get_node(mdf, &node);
                if (rc != PUMAS_RETURN_SUCCESS) break;

                if (node.key != MDF_KEY_ELEMENT)
                        continue;
                else if ((node.head == 0) && (node.tail == 1)) {
                        if (iel >= s_shared->n_elements)
                                break;
                        else
                                continue;
                }

                /* Set the element data. */
                if (iel == 0) {
                        char * tmp = (char *)s_shared->element +
                            memory_padded_size(s_shared->n_elements *
                                    sizeof(s_shared->element[0]),
                                pad_size);
                        s_shared->element[0] = (struct atomic_element *)tmp;
                } else {
                        const struct atomic_element * e =
                            s_shared->element[iel - 1];
                        s_shared->element[iel] =
                            (struct atomic_element *)((char *)(e) +
                                sizeof(struct atomic_element) + c_mem +
                                memory_padded_size(
                                    strlen(e->name) + 1, pad_size));
                }
                struct atomic_element * e = s_shared->element[iel];
                e->dcs_data = (float *)e->data;
                e->name = (char *)(e->data) + c_mem;
                memset(e->dcs_data, 0x0, c_mem);
                strcpy(e->name, node.at1.name);
                rc = PUMAS_RETURN_VALUE_ERROR;
                if (sscanf(node.at2.Z, "%lf", &(e->Z)) != 1) break;
                if (sscanf(node.at3.A, "%lf", &(e->A)) != 1) break;
                if (sscanf(node.at4.I, "%lf", &(e->I)) != 1) break;
                if ((e->Z <= 0.) || (e->A <= 0.) || (e->I <= 0.)) break;
                e->I *= 1E-09;
                rc = PUMAS_RETURN_SUCCESS;

                /* Increment. */
                iel++;
                if ((node.tail == 1) && (iel >= s_shared->n_elements)) break;
        }

        if (rc != PUMAS_RETURN_SUCCESS) {
                s_error.line = mdf->line;
                strncpy(s_error.file, mdf->mdf_path, ERROR_FILE_LENGTH);
                s_error.file[ERROR_FILE_LENGTH - 1] = '\0';
        }

        return rc;
}

/**
 * Parse the base materials from a MDF.
 *
 * @param mdf The MDF handle.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
 *
 * The materials dE/dX data are loaded according to the path provided
 * by the `<table>` XML node. If path is not absolute, it specifies a relative
 * path to the directory where the MDF is located.
 */
enum pumas_return mdf_parse_materials(struct mdf_buffer * mdf)
{
        /* Format the path directory name. */
        enum pumas_return rc;
        char * filename = NULL;
        int offset_dir, size_name;
        if (!mdf->dry_mode) {
                if ((rc = mdf_format_path(s_shared->dedx_path,
                         s_shared->mdf_path, &filename, &offset_dir,
                         &size_name)) != PUMAS_RETURN_SUCCESS)
                        return rc;
        }

        /* Loop on the XML nodes. */
        rewind(mdf->fid);
        mdf->left = 0;
        mdf->line = 1;
        mdf->depth = 0;
        struct mdf_node node;
        int imat = 0;
        const int pad_size = sizeof(*(s_shared->data));
        void * tmp_ptr = ((char *)s_shared->composition) +
            memory_padded_size(
                s_shared->n_materials * sizeof(struct material_component *),
                pad_size);
        struct material_component * data = tmp_ptr;

        for (;;) {
                /* Get the next node. */
                rc = mdf_get_node(mdf, &node);
                if (rc != PUMAS_RETURN_SUCCESS) break;

                /* Set the material data. */
                if (node.key == MDF_KEY_MATERIAL) {
                        if (node.head == 0) {
                                /* This is a material closing. */
                                s_shared->elements_in[imat] = mdf->elements_in;

                                /* Compute the relative electron density. */
                                compute_ZoA(imat);

                                /* Check for dry mode. */
                                if (mdf->dry_mode) goto update_count;

                                /* Read the energy loss data. */
                                FILE * fid = fopen(filename, "r");
                                if (fid == NULL) {
                                        mdf->line = 0;
                                        strncpy(s_error.file, filename,
                                            ERROR_FILE_LENGTH);
                                        rc = PUMAS_RETURN_PATH_ERROR;
                                        break;
                                }
                                int line;
                                rc = io_parse_dedx_file(fid, imat, &line);
                                fclose(fid);
                                if (rc != PUMAS_RETURN_SUCCESS) {
                                        mdf->line = line;
                                        strncpy(s_error.file, filename,
                                            ERROR_FILE_LENGTH);
                                        break;
                                }

                        /* Update the material count. */
                        update_count:
                                imat++;
                                if (imat >= s_shared->n_materials -
                                        s_shared->n_composites)
                                        break;
                                continue;
                        }

                        /* We have a new material opening. */
                        if (imat == 0) {
                                s_shared->material_name[0] =
                                    (char *)(s_shared->material_name +
                                        s_shared->n_materials);
                        } else {
                                s_shared->material_name[imat] =
                                    s_shared->material_name[imat - 1] +
                                    strlen(s_shared->material_name[imat - 1]) +
                                    1;
                        }
                        strcpy(s_shared->material_name[imat], node.at1.name);
                        s_shared->composition[imat] = data;

                        /* Format the energy loss filename. */
                        if (!mdf->dry_mode) {
                                const int size_new =
                                    offset_dir + strlen(node.at2.file) + 1;
                                if (size_new > size_name) {
                                        /* Get enough memory. */
                                        char * new_name =
                                            reallocate(filename, size_new);
                                        if (new_name == NULL) {
                                                rc = PUMAS_RETURN_MEMORY_ERROR;
                                                break;
                                        }
                                        filename = new_name;
                                        size_name = size_new;
                                }
                                strcpy(filename + offset_dir, node.at2.file);
                        } else {
                                int n = strlen(node.at2.file) + 1;
                                s_shared->dedx_filename[imat] = allocate(n);
                                if (s_shared->dedx_filename[imat] == NULL) {
                                        rc = PUMAS_RETURN_MEMORY_ERROR;
                                        break;
                                }
                                memcpy(s_shared->dedx_filename[imat],
                                    node.at2.file, n);
                        }
                }

                /* Skip other closings. */
                if (node.head == 0) continue;

                /* Set the composition data. */
                if (node.key == MDF_KEY_ATOMIC_COMPONENT) {
                        int i = mdf->elements_in - 1;
                        int iel = element_index(node.at1.name);
                        if (iel < 0) return PUMAS_RETURN_UNKNOWN_ELEMENT;
                        s_shared->composition[imat][i].element = iel;

                        double f;
                        if ((sscanf(node.at2.fraction, "%lf", &f) != 1) ||
                            (f <= 0.)) {
                                rc = PUMAS_RETURN_VALUE_ERROR;
                                break;
                        }
                        s_shared->composition[imat][i].fraction = f;
                        data++;
                }
        }

        if (rc != PUMAS_RETURN_SUCCESS) {
                /* Clear the filenames, whenever allocated. */
                int i;
                for (i = 0; i < s_shared->n_materials - s_shared->n_composites;
                     i++) {
                        deallocate(s_shared->dedx_filename[imat]);
                        s_shared->dedx_filename[imat] = NULL;
                }

                /* Register the error data. */
                s_error.line = mdf->line;
                if (s_error.file[0] == '\0')
                        strncpy(s_error.file, mdf->mdf_path, ERROR_FILE_LENGTH);
                s_error.file[ERROR_FILE_LENGTH - 1] = '\0';
        }

        deallocate(filename);
        return rc;
}

/**
 * Parse the composite materials from a MDF.
 *
 * @param mdf The MDF handle.
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise an error code.
 *
 * The composite materials properties are given as a linear combination
 * of the base material ones.
 */
enum pumas_return mdf_parse_composites(struct mdf_buffer * mdf)
{
        /* Loop on the XML nodes. */
        rewind(mdf->fid);
        mdf->left = 0;
        mdf->line = 1;
        mdf->depth = 0;
        struct mdf_node node;
        int icomp = 0;
        int elements_in = 0;
        const int imat = s_shared->n_materials - s_shared->n_composites - 1;
        struct material_component * data_el =
            (struct material_component *)((char *)(s_shared
                                                       ->composition[imat]) +
                s_shared->elements_in[imat] *
                    sizeof(struct material_component));
        const int pad_size = sizeof(*(s_shared->data));
        void * tmp_ptr = ((char *)s_shared->composite) +
            memory_padded_size(
                s_shared->n_composites * sizeof(struct composite_material *),
                pad_size);
        struct composite_material * data_ma = tmp_ptr;

        enum pumas_return rc;
        for (;;) {
                /* Get the next node. */
                if ((rc = mdf_get_node(mdf, &node)) != PUMAS_RETURN_SUCCESS)
                        break;

                /* Set the material data. */
                if (node.key == MDF_KEY_COMPOSITE) {
                        const int i0 = imat + icomp + 1;

                        if (node.head == 0) {
                                /* This is a composite closing. */
                                data_ma->n_components = mdf->materials_in;
                                s_shared->elements_in[i0] = elements_in;
                                data_ma =
                                    (struct composite_material *)((char *)
                                                                      data_ma +
                                        sizeof(struct composite_material) +
                                        mdf->materials_in *
                                            sizeof(struct composite_component));
                                data_el += elements_in;

                                /* update the composite material count. */
                                icomp++;
                                if (icomp >= s_shared->n_composites) break;
                                continue;
                        }

                        /* We have a new material opening. */
                        s_shared->material_name[i0] =
                            s_shared->material_name[i0 - 1] +
                            strlen(s_shared->material_name[i0 - 1]) + 1;
                        strcpy(s_shared->material_name[i0], node.at1.name);
                        s_shared->composition[i0] = data_el;
                        s_shared->composite[icomp] = data_ma;
                        elements_in = 0;
                }

                /* Skip other closings. */
                if (node.head == 0) continue;

                /* Set the composition data. */
                if (node.key == MDF_KEY_COMPOSITE_COMPONENT) {
                        int i = mdf->materials_in - 1;
                        int ima;
                        if ((rc = pumas_material_index(node.at1.name, &ima)) !=
                            PUMAS_RETURN_SUCCESS)
                                break;
                        data_ma->component[i].material = ima;
                        double f, rho;
                        rc = PUMAS_RETURN_VALUE_ERROR;
                        if (sscanf(node.at2.fraction, "%lf", &f) != 1) break;
                        if (sscanf(node.at3.density, "%lf", &rho) != 1) break;
                        if ((f < 0.) || (rho <= 0.)) break;
                        rho *= 1E+03; /* g/cm^3 -> kg/m^3.*/
                        rc = PUMAS_RETURN_SUCCESS;
                        data_ma->component[i].fraction = f;
                        data_ma->component[i].density = rho;

                        const int i0 = imat + icomp + 1;
                        for (i = 0; i < s_shared->elements_in[ima]; i++) {
                                int iel = s_shared->composition[ima][i].element;
                                int j, n = elements_in;
                                int already_listed = 0;
                                for (j = 0; j < n; j++) {
                                        if (iel ==
                                            s_shared->composition[i0][j]
                                                .element) {
                                                already_listed = 1;
                                                break;
                                        }
                                }
                                if (already_listed != 0) continue;
                                s_shared->composition[i0][elements_in++]
                                    .element = iel;
                        }
                }
        }

        if (rc == PUMAS_RETURN_END_OF_FILE) rc = PUMAS_RETURN_SUCCESS;
        if (rc != PUMAS_RETURN_SUCCESS) {
                s_error.line = mdf->line;
                strncpy(s_error.file, mdf->mdf_path, ERROR_FILE_LENGTH);
                s_error.file[ERROR_FILE_LENGTH - 1] = '\0';
        }

        return rc;
}

/**
 * Get the next XML node in the MDF.
 *
 * @param mdf  The MDF handle.
 * @param node The node handle.
 * @return On success `PUMAS_RETURN_SUCCES` otherwise the corresponding
 * error code.
 *
 * Search the file for the next valid XML node. On success, at return the *node*
 * attributes are updated with links to the XML buffer. The parsing enforces
 * consistency of the `<pumas>` node and children. Openings and closures
 * above are not checked, neither the tag names.
 */
enum pumas_return mdf_get_node(struct mdf_buffer * mdf, struct mdf_node * node)
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
                        if (mdf->left <= 0) return PUMAS_RETURN_END_OF_FILE;
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
                if (mdf->left <= 1) return PUMAS_RETURN_END_OF_FILE;
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
                                if (*tail == '\0') return PUMAS_RETURN_TOO_LONG;
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
                        enum pumas_return rc = mdf_skip_pattern(mdf, "-->");
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
                                return PUMAS_RETURN_FORMAT_ERROR;
                } else if (mdf->depth == MDF_DEPTH_ELEMENT) {
                        if (strcmp(key, "element") == 0)
                                node->key = MDF_KEY_ELEMENT;
                        else
                                return PUMAS_RETURN_FORMAT_ERROR;
                } else if (mdf->depth == MDF_DEPTH_MATERIAL) {
                        if (strcmp(key, "material") == 0)
                                node->key = MDF_KEY_MATERIAL;
                        else if (strcmp(key, "component") == 0)
                                node->key = MDF_KEY_ATOMIC_COMPONENT;
                        else
                                return PUMAS_RETURN_FORMAT_ERROR;
                } else if (mdf->depth == MDF_DEPTH_COMPOSITE) {
                        if (strcmp(key, "composite") == 0)
                                node->key = MDF_KEY_COMPOSITE;
                        else if (strcmp(key, "component") == 0)
                                node->key = MDF_KEY_COMPOSITE_COMPONENT;
                        else
                                return PUMAS_RETURN_FORMAT_ERROR;
                } else
                        return PUMAS_RETURN_FORMAT_ERROR;

                /* Check if we have an empty node. */
                if (tailler == '>')
                        goto consistency_check;
                else if ((node->head == 0) && (node->tail == 1)) {
                        /* We have a badly finished pure tailler. */
                        return PUMAS_RETURN_FORMAT_ERROR;
                }

                /* Parse any attributes. */
                for (;;) {
                        char * sep = mdf->pos;
                        while ((*sep != '=') && (*sep != '>')) {
                                if (*sep == '\0') return PUMAS_RETURN_TOO_LONG;
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
                        if (value == NULL) return PUMAS_RETURN_FORMAT_ERROR;
                        value++;
                        end = strchr(value, '"');
                        if (end == NULL) return PUMAS_RETURN_FORMAT_ERROR;
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
                                        return PUMAS_RETURN_FORMAT_ERROR;
                        } else if (node->key == MDF_KEY_MATERIAL) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else if (strcmp(attr, "file") == 0)
                                        node->at2.file = value;
                                else
                                        return PUMAS_RETURN_FORMAT_ERROR;
                        } else if (node->key == MDF_KEY_COMPOSITE) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else
                                        return PUMAS_RETURN_FORMAT_ERROR;
                        } else if (node->key == MDF_KEY_ATOMIC_COMPONENT) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else if (strcmp(attr, "fraction") == 0)
                                        node->at2.fraction = value;
                                else
                                        return PUMAS_RETURN_FORMAT_ERROR;
                        } else if (node->key == MDF_KEY_COMPOSITE_COMPONENT) {
                                if (strcmp(attr, "name") == 0)
                                        node->at1.name = value;
                                else if (strcmp(attr, "fraction") == 0)
                                        node->at2.fraction = value;
                                else if (strcmp(attr, "density") == 0)
                                        node->at3.density = value;
                                else
                                        return PUMAS_RETURN_FORMAT_ERROR;
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
                        if (node->tail == 1) return PUMAS_RETURN_FORMAT_ERROR;
                        mdf->depth = MDF_DEPTH_ROOT;
                }
        } else if (mdf->depth == MDF_DEPTH_ROOT) {
                if (node->key == MDF_KEY_PUMAS) {
                        if (node->head == 1) return PUMAS_RETURN_FORMAT_ERROR;
                        mdf->depth = MDF_DEPTH_EXTERN;
                } else if (node->key == MDF_KEY_ELEMENT) {
                        if (node->tail == 0) mdf->depth = MDF_DEPTH_ELEMENT;
                } else if (node->key == MDF_KEY_MATERIAL) {
                        if (node->tail == 1) return PUMAS_RETURN_FORMAT_ERROR;
                        mdf->depth = MDF_DEPTH_MATERIAL;
                } else if (node->key == MDF_KEY_COMPOSITE) {
                        if (node->tail == 1) return PUMAS_RETURN_FORMAT_ERROR;
                        mdf->depth = MDF_DEPTH_COMPOSITE;
                }
        } else if (mdf->depth == MDF_DEPTH_ELEMENT) {
                if (node->key == MDF_KEY_ELEMENT) {
                        if (node->head == 1) return PUMAS_RETURN_FORMAT_ERROR;
                        mdf->depth = MDF_DEPTH_ROOT;
                }
        } else if (mdf->depth == MDF_DEPTH_MATERIAL) {
                if (node->key == MDF_KEY_MATERIAL) {
                        if (node->head == 1) return PUMAS_RETURN_FORMAT_ERROR;
                        mdf->depth = MDF_DEPTH_ROOT;
                }
        } else if (mdf->depth == MDF_DEPTH_COMPOSITE) {
                if (node->key == MDF_KEY_COMPOSITE) {
                        if (node->head == 1) return PUMAS_RETURN_FORMAT_ERROR;
                        mdf->depth = MDF_DEPTH_ROOT;
                }
        } else
                return PUMAS_RETURN_FORMAT_ERROR;

        /* Check that the node has all its attributes defined. */
        if (node->head == 0) {
                if ((node->key == MDF_KEY_MATERIAL) && (mdf->elements_in == 0))
                        return PUMAS_RETURN_FORMAT_ERROR;
                else if ((node->key == MDF_KEY_COMPOSITE) &&
                    (mdf->materials_in == 0))
                        return PUMAS_RETURN_FORMAT_ERROR;
                return PUMAS_RETURN_SUCCESS;
        }

        if (node->key == MDF_KEY_COMPOSITE) {
                if (node->at1.name == NULL) return PUMAS_RETURN_FORMAT_ERROR;
        } else if ((node->key == MDF_KEY_MATERIAL) ||
            (node->key == MDF_KEY_ATOMIC_COMPONENT)) {
                if ((node->at1.name == NULL) || (node->at2.fraction == NULL))
                        return PUMAS_RETURN_FORMAT_ERROR;
        } else if (node->key == MDF_KEY_ELEMENT) {
                if ((node->at1.name == NULL) || (node->at2.Z == NULL) ||
                    (node->at3.A == NULL) || (node->at4.I == NULL))
                        return PUMAS_RETURN_FORMAT_ERROR;
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
 * @return On success `PUMAS_RETURN_SUCCES` otherwise the corresponding
 * error code.
 *
 * Search the file for *pattern* and skip it. This routine is used in
 * order to skip XML comments.
 */
enum pumas_return mdf_skip_pattern(
    struct mdf_buffer * mdf, const char * pattern)
{
        const int pattern_size = strlen(pattern);
        for (;;) {
                if (mdf->left < pattern_size) {
                        memmove(mdf->data, mdf->pos, mdf->left);
                        mdf->left += fread(mdf->data + mdf->left, sizeof(char),
                            mdf->size - 1 - mdf->left, mdf->fid);
                        if (mdf->left < pattern_size)
                                return PUMAS_RETURN_END_OF_FILE;
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
 * @return On success `PUMAS_RETURN_SUCCESS` otherwise `PUMAS_ERROR`.
 *
 * The full path to the filename string is formated by prepending the
 * relative or absolute directory path. If directory is a relative path,
 * the directory of the MDF is used as root.
 */
enum pumas_return mdf_format_path(const char * directory, const char * mdf_path,
    char ** filename, int * offset_dir, int * size_name)
{
        /* Format the path directory name. */
        const char sep =
#if (defined _WIN32) || (defined WIN32) || (defined __CYGWIN__)
            '\\';
#else
            '/';
#endif
        const int dir_size = strlen(directory);
        const int initial_size = 128;
        if (directory[0] == sep ||
            ((dir_size >= 3) && (directory[1] == ':') && directory[2] == sep)) {
                /* We have an absolute path name. */
                *offset_dir = dir_size + 1;
                *size_name = (*offset_dir) + initial_size;
                *filename = allocate(*size_name);
                if (*filename == NULL) return PUMAS_RETURN_MEMORY_ERROR;
                strcpy(*filename, directory);
                (*filename)[(*offset_dir) - 1] = sep;
        } else {
                /* We have a relative path name. */
                int n1 = strlen(mdf_path) - 1;
                while ((n1 >= 0) && (mdf_path[n1] != sep)) n1--;
                *offset_dir = n1 + dir_size + 2;
                *size_name = (*offset_dir) + initial_size;
                *filename = allocate(*size_name);
                if (*filename == NULL) return PUMAS_RETURN_MEMORY_ERROR;
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
 * @return `PUMAS_ERROR` if any computation failled, `PUMAS_RETURN_SUCCESS`
 * otherwise.
 */
enum pumas_return compute_composite(int material)
{
        compute_composite_weights(material);
        compute_composite_tables(material);
        compute_regularise_del(material);
        int ikin;
        for (ikin = 0; ikin < s_shared->n_kinetics; ikin++) {
                const enum pumas_return rc =
                    compute_coulomb_parameters(material, ikin);
                if (rc != PUMAS_RETURN_SUCCESS) return rc;
        }
        compute_cel_integrals(material);
        compute_csda_magnetic_transport(material);
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Compute various cumulative integrals for a deterministic CEL.
 *
 * @param material The index of the material to tabulate.
 */
void compute_cel_integrals(int material)
{
        compute_cel_grammage_integral(0, material);
        compute_cel_grammage_integral(1, material);
        compute_time_integrals(material);
        compute_kinetic_integral(
            table_get_NI_el(PUMAS_SCHEME_CSDA, material, 0));
        compute_kinetic_integral(
            table_get_NI_el(PUMAS_SCHEME_HYBRID, material, 0));
        compute_kinetic_integral(table_get_NI_in(material, 0));
}

/**
 * Computation of mixture weights for composite materials.
 *
 * @param material The material index of the composite to compute.
 */
void compute_composite_weights(int material)
{
        /* Compute the effective mass number (normalisation). */
        const int icomp =
            material - s_shared->n_materials + s_shared->n_composites;
        int i;
        double M = 0.;
        for (i = 0; i < s_shared->composite[icomp]->n_components; i++) {
                const struct composite_component * component =
                    s_shared->composite[icomp]->component + i;
                const int imat = component->material;
                int j;
                double Mi_inv = 0.;
                for (j = 0; j < s_shared->elements_in[imat]; j++) {
                        const struct material_component * component =
                            s_shared->composition[imat] + j;
                        const struct atomic_element * element =
                            s_shared->element[component->element];
                        Mi_inv += component->fraction / element->A;
                }
                M += component->fraction * Mi_inv;
        }
        M = 1. / M;

        /* Compute the materials and atomic elements weights. */
        for (i = 0; i < s_shared->elements_in[material]; i++) {
                struct material_component * component =
                    s_shared->composition[material] + i;
                component->fraction = 0.;
        }

        for (i = 0; i < s_shared->composite[icomp]->n_components; i++) {
                struct composite_component * component =
                    s_shared->composite[icomp]->component + i;
                const int imat = component->material;
                int j;
                double Mi_inv = 0.;
                for (j = 0; j < s_shared->elements_in[imat]; j++) {
                        const struct material_component * component =
                            s_shared->composition[imat] + j;
                        const struct atomic_element * element =
                            s_shared->element[component->element];
                        Mi_inv += component->fraction / element->A;
                }
                component->weight = component->fraction * M * Mi_inv;

                for (j = 0; j < s_shared->elements_in[imat]; j++) {
                        const struct material_component * cij =
                            s_shared->composition[imat] + j;
                        int k;
                        for (k = 0; k < s_shared->elements_in[material]; k++) {
                                struct material_component * c =
                                    s_shared->composition[material] + k;
                                if (c->element == cij->element) {
                                        c->fraction +=
                                            component->weight * cij->fraction;
                                        break;
                                }
                        }
                }
        }

        /* Compute the relative electron density. */
        compute_ZoA(material);
}

/**
 * Compute the density of a composite material.
 *
 * @param material The material index of the composite material.
 * @return `PUMAS_ERROR` if any density is not strictly positive,
 * `PUMAS_RETURN_SUCCESS` otherwise.
 */
enum pumas_return compute_composite_density(int material)
{
        const int icomp =
            material - s_shared->n_materials + s_shared->n_composites;
        int i;
        double rho_inv = 0., nrm = 0.;
        for (i = 0; i < s_shared->composite[icomp]->n_components; i++) {
                const struct composite_component * component =
                    s_shared->composite[icomp]->component + i;
                if (component->density <= 0.) return PUMAS_RETURN_DENSITY_ERROR;
                rho_inv += component->fraction / component->density;
                nrm += component->fraction;
        }

        if (rho_inv <= 0.) return PUMAS_RETURN_DENSITY_ERROR;
        s_shared->composite[icomp]->density = nrm / rho_inv;
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Computation of tabulated properties for composite materials.
 *
 * @param material The matrial index of the composite to compute.
 */
void compute_composite_tables(int material)
{
        const int icomp =
            material - s_shared->n_materials + s_shared->n_composites;

        /* Initialise to zero. */
        int i, row, k0 = 0;
        for (i = 0; i < material; i++) k0 += s_shared->elements_in[i];
        row = 0;
        *table_get_T(PUMAS_SCHEME_CSDA, material, row) = 0.;
        *table_get_T(PUMAS_SCHEME_HYBRID, material, row) = 0.;
        *table_get_NI_in(material, row) = 0.;
        for (row = 0; row < s_shared->n_kinetics; row++) {
                *table_get_dE(PUMAS_SCHEME_CSDA, material, row) = 0.;
                *table_get_dE(PUMAS_SCHEME_HYBRID, material, row) = 0.;
                *table_get_CS(material, row) = 0.;
                int k, ip;
                for (k = 0; k < s_shared->elements_in[material]; k++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++)
                                *table_get_CSf(ip, k0 + k, row) = 0.;
        }
        *table_get_Kt(material) = 0.;
        *table_get_a_max(material) = 0.;
        *table_get_b_max(PUMAS_SCHEME_CSDA, material) = 0.;
        *table_get_b_max(PUMAS_SCHEME_HYBRID, material) = 0.;

        /* End point statistics */
        for (i = 0; i < s_shared->composite[icomp]->n_components; i++) {
                struct composite_component * component =
                    s_shared->composite[icomp]->component + i;
                const int imat = component->material;
                const double kt = *table_get_Kt(imat);

                /* Total cross section and energy loss. */
                for (row = 0; row < s_shared->n_kinetics; row++) {
                        *table_get_dE(0, material, row) +=
                            *table_get_dE(0, imat, row) * component->weight;
                        *table_get_dE(1, material, row) +=
                            *table_get_dE(1, imat, row) * component->weight;
                        const double k = *table_get_K(row);
                        const double cs = (k < kt) ?
                            0. :
                            *table_get_CS(imat, row) * component->weight;
                        *table_get_CS(material, row) += cs;
                }

                /* Fractional contribution to the cross section. */
                int j, j0 = 0;
                for (j = 0; j < imat; j++) j0 += s_shared->elements_in[j];

                for (row = 0; row < s_shared->n_kinetics; row++) {
                        const double kinetic = *table_get_K(row);
                        if (kinetic < kt) continue;

                        const double cs_tot = *table_get_CS(imat, row);
                        double csf_last = 0.;
                        for (j = 0; j < s_shared->elements_in[imat]; j++) {
                                /* Locate the element. */
                                const struct material_component * cij =
                                    s_shared->composition[imat] + j;
                                int k;
                                for (k = 0; k < s_shared->elements_in[material];
                                     k++) {
                                        struct material_component * c =
                                            k + s_shared->composition[material];
                                        if (c->element == cij->element) break;
                                }

                                /* Update the fractional cross-section. */
                                int ip;
                                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                        const double csf =
                                            *table_get_CSf(ip, j0 + j, row);
                                        *table_get_CSf(ip, k0 + k, row) +=
                                            cs_tot * (csf - csf_last) *
                                            component->weight;
                                        csf_last = csf;
                                }
                        }
                }

                /* Maximum tabulated energy loss parameters. */
                *table_get_a_max(material) +=
                    *table_get_a_max(imat) * component->weight;
                *table_get_b_max(0, material) +=
                    *table_get_b_max(0, imat) * component->weight;
                *table_get_b_max(1, material) +=
                    *table_get_b_max(1, imat) * component->weight;
        }

        /* Normalise the fractional contributions to the cross section. */
        for (row = 0; row < s_shared->n_kinetics; row++) {
                int k, ip;
                const double cs = *table_get_CS(material, row);
                if (cs <= 0.) {
                        for (k = 0; k < s_shared->elements_in[material]; k++)
                                for (ip = 0; ip < N_DEL_PROCESSES; ip++)
                                        *table_get_CSf(ip, k0 + k, row) = 0.;
                        continue;
                }

                double sum = 0.;
                for (k = 0; k < s_shared->elements_in[material]; k++)
                        for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                                sum += *table_get_CSf(ip, k0 + k, row);
                                const double f = sum / cs;
                                *table_get_CSf(ip, k0 + k, row) =
                                    (f > 1.) ? 1. : f;
                        }

                /* Protect against rounding errors. */
                k = s_shared->elements_in[material] - 1;
                *table_get_CSf(N_DEL_PROCESSES - 1, k0 + k, row) = 1.;
        }

        /* Weighted integrands */
        for (row = 1; row < s_shared->n_kinetics; row++) {
                *table_get_NI_in(material, row) = *table_get_CS(material, row) /
                    *table_get_dE(PUMAS_SCHEME_HYBRID, material, row);
        }
}

/**
 * Compute the cumulative integral of a table column.
 *
 * @param table The table column to process.
 *
 * Compute a cumulative integral of the column by integrating over the kinetic
 * energy with a linear interpolation and a trapezoidal rule.
 */
void compute_kinetic_integral(double * table)
{
        double value = 0., dv;
        int i;
        for (i = 1; i < s_shared->n_kinetics; i++) {
                const double x0 = *table_get_K(i - 1);
                const double x1 = *table_get_K(i);
                dv = 0.5 * (x1 - x0) * (table[i - 1] + table[i]);
                table[i - 1] = value;
                value += dv;
        }
        table[s_shared->n_kinetics - 1] = value;
}

/**
 * Compute cumulative integrals for the proper time.
 *
 * @param table The table column to process.
 *
 * Compute the proper time cumulative integrals over path length with a
 * trapezoidal rule.
 */
void compute_time_integrals(int material)
{
        static double I0 = 0.;
        if (I0 == 0.) {
                /* Compute the integral of 1/momemtum for the lowest energy bin
                 * using trapezes. */
                const int n = 101;
                int i;
                const double dK = (*table_get_K(1)) / (n - 1);
                double Ki = dK;
                I0 = 0.5 / sqrt(Ki * (Ki + 2. * s_shared->mass));
                for (i = 2; i < n - 1; i++) {
                        Ki += dK;
                        const double pi = sqrt(Ki * (Ki + 2. * s_shared->mass));
                        I0 += 1. / pi;
                }
                Ki += dK;
                I0 += 0.5 / sqrt(Ki * (Ki + 2. * s_shared->mass));
                I0 /= n - 1;
        }

        /* Compute the cumulative path integrals . */
        double * const K = table_get_K(0);
        double * const T0 = table_get_T(PUMAS_SCHEME_CSDA, material, 0);
        double * const T1 = table_get_T(PUMAS_SCHEME_HYBRID, material, 0);
        double * const X0 = table_get_X(PUMAS_SCHEME_CSDA, material, 0);
        double * const X1 = table_get_X(PUMAS_SCHEME_HYBRID, material, 0);

        T0[0] = T1[0] = 0.;
        T0[1] = I0 * X0[1] * s_shared->mass;
        T1[1] = I0 * X1[1] * s_shared->mass;
        int i;
        for (i = 2; i < s_shared->n_kinetics; i++) {
                const double p0 =
                    sqrt(K[i - 1] * (K[i - 1] + 2. * s_shared->mass));
                const double p1 = sqrt(K[i] * (K[i] + 2. * s_shared->mass));
                const double psi = 1. / p0 + 1. / p1;
                const double dy0 = 0.5 * (X0[i] - X0[i - 1]) * psi;
                const double dy1 = 0.5 * (X1[i] - X1[i - 1]) * psi;
                T0[i] = T0[i - 1] + dy0 * s_shared->mass;
                T1[i] = T1[i - 1] + dy1 * s_shared->mass;
        }
}

/**
 * Compute the CSDA grammage range from the energy loss.
 *
 *  @param scheme   The index of the simulation scheme.
 *  @param material The index of the material to process.
 *
 * Compute the cumulative CSDA grammage integral from the energy loss
 * using a trapezoidal rule.
 */
void compute_cel_grammage_integral(int scheme, int material)
{
        const double * const kinetic = table_get_K(0);
        const double * const dEdX = table_get_dE(scheme, material, 0);
        double * const table = table_get_X(scheme, material, 0);

        /* Compute the cumulative integral. */
        int i;
        table[0] = 0.;
        double y0 = 1. / dEdX[1];
        for (i = 1; i < s_shared->n_kinetics; i++) {
                const double y1 = 1. / dEdX[i];
                table[i] = table[i - 1] +
                    0.5 * (kinetic[i] - kinetic[i - 1]) * (y0 + y1);
                y0 = y1;
        }
}

/**
 * Compute the cumulative integrals of the magnetic transport.
 *
 * @param imed The index of the material to tabulate.
 *
 * Compute the cumulative integrals for the momenta of the magnetic deflection
 * using a trapezoidal rule.
 */
void compute_csda_magnetic_transport(int material)
{
        double x[N_LARMOR_ORDERS];
        double dx[N_LARMOR_ORDERS];
        memset(x, 0x0, sizeof(x));

        /* The magnetic phase shift is proportional to the proper time integral.
         * We refer to this table. */
        const double factor = LARMOR_FACTOR / s_shared->mass;
        double * const T = table_get_T(PUMAS_SCHEME_CSDA, material, 0);

        /* Compute the deflection starting from max energy down to 0 */
        double * const X0 = table_get_X(PUMAS_SCHEME_CSDA, material, 0);
        const int imax = s_shared->n_kinetics - 1;
        int i, j;
        for (i = s_shared->n_kinetics - 2; i >= 1; i--) {
                double dX0 = 0.5 * (X0[i + 1] - X0[i]);
                double p1 = (T[imax] - T[i]) * factor;
                double p2 = (T[imax] - T[i + 1]) * factor;

                double f1 = 1., f2 = 1.;
                for (j = 0; j < N_LARMOR_ORDERS; j++) {
                        *table_get_Li(material, j, i + 1) = x[j];
                        dx[j] = dX0 * (f1 + f2);
                        x[j] += dx[j];
                        f1 *= p1;
                        f2 *= p2;
                }
        }

        /* Extrapolate the end points */
        for (j = 0; j < N_LARMOR_ORDERS; j++) {
                *table_get_Li(material, j, 1) = x[j];
                double hx = (X0[1] - X0[0]) / (X0[2] - X0[1]);
                *table_get_Li(material, j, 0) = x[j] + hx * dx[j];
        }
}

/**
 * Compute and tabulate the multiple scattering parameters.
 *
 * @param material The target material.
 * @param row The kinetic energy index to compute for.
 *
 * At output the *row* of `s_shared::table_Mu0` and `s_shared::table_Ms1` are
 * updated. This routine manages a temporary workspace memory buffer. Calling
 * it with a negative *material* index causes the workspace memory to be freed.
 */
enum pumas_return compute_coulomb_parameters(int material, int row)
{
        /* Handle the memory for the temporary workspace. */
        static struct coulomb_workspace * workspace = NULL;
        if (material < 0) {
                deallocate(workspace);
                workspace = NULL;
                return PUMAS_RETURN_SUCCESS;
        } else if (workspace == NULL) {
                const int work_size = sizeof(struct coulomb_workspace) +
                    s_shared->max_components * sizeof(struct coulomb_data);
                workspace = allocate(work_size);
                if (workspace == NULL) return PUMAS_RETURN_MEMORY_ERROR;
        }

        /* Check the kinetic energy. */
        const double kinetic = *table_get_K(row);
        if (kinetic <= 0.) {
                *table_get_Mu0(material, row) = 0.;
                *table_get_Ms1(material, row) = 0.;
                return PUMAS_RETURN_SUCCESS;
        }

        /* Compute the mean free paths, the screening factors and various
         * averages used for the hard scattering.
         */
        double invlb_m = 0., invlb1_m = 0.;
        double s_m_l = 0., s_m_h = 0.;
        int i;
        struct coulomb_data * data;
        for (i = 0, data = workspace->data; i < s_shared->elements_in[material];
             i++, data++) {
                double G[2];
                const struct material_component * const component =
                    &s_shared->composition[material][i];
                const struct atomic_element * const element =
                    s_shared->element[component->element];
                double kinetic0;
                coulomb_frame_parameters(
                    kinetic, element->A, &kinetic0, data->fCM);
                data->fspin = coulomb_spin_factor(kinetic0);
                coulomb_screening_parameters(
                    NULL, kinetic0, component->element, data->screening);
                coulomb_pole_decomposition(data->screening, data->a, data->b);
                coulomb_transport_coefficients(
                    1., data->fspin, data->screening, data->a, data->b, G);
                const double invlb = component->fraction /
                    coulomb_wentzel_path(
                        kinetic0, element->Z, element->A, data->screening[0]);

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
                double fmax = compute_cutoff_objective(mu_max, workspace);
                double fmin;
                if (fmax > 0.) {
                        /* This shouldn't occur, but let's be safe and handle
                         * this case. */
                        mu_min = mu_max;
                        fmin = fmax;
                        mu_max = 1.;
                        fmax = -workspace->cs_h;
                } else {
                        fmin = compute_cutoff_objective(mu_min, workspace);
                        if (fmin < 0.) {
                                /* This might occur at high energies when the
                                 * nuclear screening becomes significant. */
                                mu_max = mu_min;
                                fmax = fmin;
                                mu_min = 0.;
                                fmin =
                                    compute_cutoff_objective(mu_min, workspace);
                        }
                }
                if (mu_min < max_mu0) {
                        if (mu_max > max_mu0) mu_max = max_mu0;
                        if (math_find_root(compute_cutoff_objective, mu_min,
                                mu_max, &fmin, &fmax, 1E-06 * mu0, 1E-06, 100,
                                workspace, &mubest) == 0)
                                mu0 = mubest;
                }
                if (mu0 > max_mu0) mu0 = max_mu0;
                lb_h =
                    compute_cutoff_objective(mu0, workspace) + workspace->cs_h;
                if (lb_h <= 1. / EHS_PATH_MAX)
                        lb_h = EHS_PATH_MAX;
                else
                        lb_h = 1. / lb_h;
        } else
                lb_h = lb_m;
        *table_get_Mu0(material, row) = mu0;
        *table_get_Lb(material, row) =
            lb_h * kinetic * (kinetic + 2. * s_shared->mass);
        *table_get_NI_el(PUMAS_SCHEME_CSDA, material, row) =
            1. / (*table_get_dE(PUMAS_SCHEME_CSDA, material, row) * lb_h);
        *table_get_NI_el(PUMAS_SCHEME_HYBRID, material, row) =
            1. / (*table_get_dE(PUMAS_SCHEME_HYBRID, material, row) * lb_h);

        /* Compute the 1st moment of the soft scattering. */
        const int n0 = s_shared->n_materials - s_shared->n_composites;
        if (material < n0) {
                /* We have a base material. */
                static double * ms1_table = NULL;
                if (material == 0) {
                        /* Precompute the per element soft scattering terms. */
                        enum pumas_return rc;
                        if ((rc = compute_coulomb_soft(row, &ms1_table)) !=
                            PUMAS_RETURN_SUCCESS)
                                return rc;
                }

                double invlb1 = 0.;
                struct material_component * component =
                    s_shared->composition[material];
                for (i = 0, data = workspace->data;
                     i < s_shared->elements_in[material];
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
                        invlb1 += *table_get_ms1(iel, row, ms1_table) *
                            component->fraction;
                }
                *table_get_Ms1(material, row) = invlb1;
        } else {
                /* We have a composite material. Let's loop on the base
                 * material components.
                 */
                const struct composite_material * composite =
                    s_shared->composite[material - n0];
                double invlb1 = 0.;
                int icomp;
                for (icomp = 0; icomp < composite->n_components; icomp++) {
                        const struct composite_component * c =
                            composite->component + icomp;
                        const int imat = c->material;
                        struct material_component * component =
                            s_shared->composition[imat];

                        /* Compute the variation of the Coulomb scattering. */
                        const double mu0b = *table_get_Mu0(imat, row);
                        double delta_invlb1 = 0.;
                        for (i = 0, data = workspace->data;
                             i < s_shared->elements_in[imat];
                             i++, data++, component++) {
                                /* Compute the scattering parameters for the
                                 * base material.
                                 */
                                double G[2];
                                const struct atomic_element * const element =
                                    s_shared->element[component->element];
                                double kinetic0;
                                coulomb_frame_parameters(
                                    kinetic, element->A, &kinetic0, data->fCM);
                                data->fspin = coulomb_spin_factor(kinetic0);
                                coulomb_screening_parameters(NULL, kinetic0,
                                    component->element, data->screening);
                                coulomb_pole_decomposition(
                                    data->screening, data->a, data->b);
                                coulomb_transport_coefficients(1., data->fspin,
                                    data->screening, data->a, data->b, G);
                                const double invlb = component->fraction /
                                    coulomb_wentzel_path(kinetic0, element->Z,
                                        element->A, data->screening[0]);

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
                        invlb1 += (*table_get_Ms1(imat, row) + delta_invlb1) *
                            c->weight;
                }
                *table_get_Ms1(material, row) = invlb1;
        }

        return PUMAS_RETURN_SUCCESS;
}

/**
 * Tabulate the multiple scattering per element.
 *
 * @param row     The row index for the kinetic value.
 * @param data    The tabulated data.
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
enum pumas_return compute_coulomb_soft(int row, double ** data)
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
                ms1_table = allocate(s_shared->n_elements *
                    s_shared->n_kinetics * sizeof(double));
                if (ms1_table == NULL) return PUMAS_RETURN_MEMORY_ERROR;
        }

        /* Loop over atomic elements. */
        const double kinetic = *table_get_K(row);
        int iel;
        for (iel = 0; iel < s_shared->n_elements; iel++) {
                struct atomic_element * element = s_shared->element[iel];
                double invlb1 = 0.;

                /* Electron shells contribution to the transverse transport. */
                invlb1 = transverse_transport_ionisation(element, kinetic);

                /* Photonuclear contribution to the transverse transport. */
                invlb1 += transverse_transport_photonuclear(element, kinetic);

                *table_get_ms1(iel, row, ms1_table) = invlb1;
        }

        *data = ms1_table;
        return PUMAS_RETURN_SUCCESS;
}

/**
 * Objective function for solving the soft scattering cut-off angle.
 *
 * @param mu         The proposed cutoff angular parameter.
 * @param parameters A pointer to the temporary workspace.
 * @return The difference between the current and expected restricted cross
 * section.
 *
 * This is a wrapper for the root solver. It provides the objective function
 * to solve for.
 */
double compute_cutoff_objective(double mu, void * parameters)
{
        /* Unpack the workspace. */
        struct coulomb_workspace * workspace = parameters;

        /* Compute the restricted cross section. */
        double cs_tot = 0.;
        int i, n = s_shared->elements_in[workspace->material];
        struct coulomb_data * data;
        for (i = 0, data = workspace->data; i < n; i++, data++) {
                cs_tot += data->invlambda *
                    coulomb_restricted_cs(
                        mu, data->fspin, data->screening, data->a, data->b);
        }

        /* Return the difference with the expectation. */
        return cs_tot - workspace->cs_h;
}

/**
 * Tabulate the deterministic CEL and DEL cross-sections per element.
 *
 * @param row the row index for the kinetic value.
 * @return The tabulated data.
 *
 * Because this step is time consuming, the CEL integration is precomputed per
 * element at initialisation using a temporary buffer. These temporary tables
 * are used for computing per material CEL corrections and DEL cross-sections.
 *
 * **Note** This routine handles a static dynamically allocated table. If the
 * *row* index is negative the table is freed.
 */
double * compute_cel_and_del(int row)
{
        static double * cel_table = NULL;

        if (row < 0) {
                deallocate(cel_table);
                return (cel_table = NULL);
        }

        /* Allocate the temporary table. */
        if (cel_table == NULL) {
                cel_table = allocate(N_DEL_PROCESSES * s_shared->n_elements *
                    s_shared->n_kinetics * sizeof(double));
                if (cel_table == NULL) return NULL;
        }

        /* Loop over atomic elements. */
        const double kinetic = *table_get_K(row);
        int iel;
        for (iel = 0; iel < s_shared->n_elements; iel++) {
                const struct atomic_element * element = s_shared->element[iel];

                /* Loop over processes. */
                int ip;
                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                        *table_get_CSn(ip, iel, row) = compute_dcs_integral(
                            0, element, kinetic, dcs_get(ip), X_FRACTION, 180);
                        *table_get_cel(ip, iel, row, cel_table) =
                            compute_dcs_integral(1, element, kinetic,
                                dcs_get(ip), X_FRACTION, 180);
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
void compute_regularise_del(int material)
{
        /* Find the kinetic threshold above which the total cross-section is
         * not null. */
        int it;
        double cs0;
        for (it = 1; it < s_shared->n_kinetics; it++)
                if ((cs0 = *table_get_CS(material, it)) != 0) break;
        *table_get_Kt(material) = *table_get_K(it);

        /*
         * Regularise the cross section and the total number of interaction
         * lengths.
         */
        int row;
        for (row = 1; row < it; row++) {
                *table_get_CS(material, row) = cs0;
                const double dEdX =
                    *table_get_dE(PUMAS_SCHEME_HYBRID, material, row);
                *table_get_NI_in(material, row) = cs0 / dEdX;
        }

        if (material > 0)
                return; /* The computations below need to be done only
once. */

        /* Find the kinetic threshold per process and atomic element. */
        int iel, ip;
        for (iel = 0; iel < s_shared->n_elements; iel++) {
                const struct atomic_element * element = s_shared->element[iel];
                for (ip = 0; ip < N_DEL_PROCESSES; ip++) {
                        for (row = 0; row < it; row++)
                                *table_get_Xt(ip, iel, row) = 1.;
                        dcs_function_t * dcs_func = dcs_get(ip);
                        for (row = it; row < s_shared->n_kinetics; row++) {
                                const double k = *table_get_K(row);
                                double x = X_FRACTION;
                                while ((x < 1.) &&
                                    (dcs_func(element, k, k * x) <= 0.))
                                        x *= 2;
                                if (x >= 1.)
                                        x = 1.;
                                else if (x > X_FRACTION) {
                                        const double eps = 1E-02 * X_FRACTION;
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
                                                dcs = dcs_func(
                                                    element, k, k * x0);
                                        }
                                }
                                *table_get_Xt(ip, iel, row) = x;
                        }
                }
        }
}

/**
 * Compute integrals of DCSs.
 *
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
double compute_dcs_integral(int mode, const struct atomic_element * element,
    double kinetic, dcs_function_t * dcs, double xlow, int nint)
{

        /* Let us use the analytical form for ionisation when radiative
         * corrections can be neglected.
         */
        const double m1 = s_shared->mass - ELECTRON_MASS;
        if ((dcs == &dcs_ionisation) &&
            (kinetic <= 0.5 * m1 * m1 / ELECTRON_MASS))
                return dcs_ionisation_integrate(mode, element, kinetic, xlow);

        /* We integrate over the recoil energy using a logarithmic sampling. */
        double dcsint = 0.;
        double x0 = log(kinetic * xlow), x1 = log(kinetic);
        math_gauss_quad(nint, &x0, &x1); /* Initialisation. */

        double xi, wi;
        while (math_gauss_quad(0, &xi, &wi) == 0) { /* Iterations. */
                const double qi = exp(xi);
                double y = dcs(element, kinetic, qi) * qi;
                if (mode != 0) y *= qi;
                dcsint += y * wi;
        }
        dcsint /= kinetic + s_shared->mass;

        return dcsint;
}

/**
 * Compute or update the relative electron density of a material.
 *
 * @param material The material index.
 */
void compute_ZoA(int material)
{
        int i;
        double ZoA = 0.;
        for (i = 0; i < s_shared->elements_in[material]; i++) {

                const struct material_component * const c =
                    s_shared->composition[material] + i;
                const struct atomic_element * const e =
                    s_shared->element[c->element];
                ZoA += e->Z / e->A * c->fraction;
        }
        s_shared->material_ZoA[material] = ZoA;
}

/**
 * Compute the coefficients of the polynomial approximation to an
 * inelastic DCS.
 *
 * @param [in] dcs_func The DCS function.
 * @param [in] element  The target atomic element.
 *
 * The DCS is also tabulated at specific values for the ziggurat algorithm.
 */
enum pumas_return compute_dcs_model(
    dcs_function_t * dcs_func, struct atomic_element * element)
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
                m = (int)(100. * log10(DCS_MODEL_MAX_FRACTION / X_FRACTION)) +
                    1 + DCS_SAMPLING_N;
                if (m < 0) return PUMAS_RETURN_INTERNAL_ERROR;
                tmp = allocate((3 * m + n + DCS_SAMPLING_N) * sizeof(double));
                if (tmp == NULL) return PUMAS_RETURN_MEMORY_ERROR;
        }
        double * x = tmp;
        double * y = x + m;
        double * w = y + m;
        double * c_ = w + m;

        /*  Loop over the kinetic energy values. */
        const int index = dcs_get_index(dcs_func);
        const int nkeff = s_shared->n_kinetics - s_shared->dcs_model_offset;
        float * const coeff = table_get_dcs_coeff(element, index, 0);
        const double x0 = log(X_FRACTION);
        const double dx =
            log(DCS_MODEL_MAX_FRACTION / X_FRACTION) / (m - 1 - DCS_SAMPLING_N);
        int i;
        float * c;
        for (i = 0, c = coeff; i < nkeff; i++, c += n + DCS_SAMPLING_N) {
                /* Prepare the fit values using a log sampling. */
                const double K = *table_get_K(i + s_shared->dcs_model_offset);
                int j, first = 1, j0 = -1, j1 = -1;
                for (j = 0; j < m - DCS_SAMPLING_N; j++) {
                        const double nu = exp(x0 + j * dx);
                        x[j] = nu;
                        const double dcs = dcs_func(element, K, nu * K) * K /
                            (K + s_shared->mass);
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
                if ((j0 < 0) || (j1 < 0)) return PUMAS_RETURN_INTERNAL_ERROR;

                /* Add the tabulated values in linear scale. */
                const double dnu = (1. - x[j0]) / DCS_SAMPLING_N;
                double nu;
                for (j = m - DCS_SAMPLING_N, nu = x[j0]; j - m < 0;
                     j++, nu += dnu)
                /* Patch a gcc warning here. With -O2 j < m != j-m < 0 ... */
                {
                        x[j] = nu;
                        const double tmp = dcs_func(element, K, nu * K) * K /
                            (K + s_shared->mass);
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
 * The Bremsstrahlung differential cross section.
 *
 * @param element The target atomic element.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The differential cross section in m^2/kg.
 *
 * The differential cross-section is computed following the PDG:
 * http://pdg.lbl.gov/2014/AtomicNuclearProperties/adndt.pdf.
 */
double dcs_bremsstrahlung(
    const struct atomic_element * element, double K, double q)
{
        const double Z = element->Z;
        const double A = element->A;
        const double mu = s_shared->mass;
        const double me = ELECTRON_MASS;
        const double sqrte = 1.648721271;
        const double phie_factor = mu / (me * me * sqrte);
        const double rem = 5.63588E-13 * me / mu;

        const double BZ_n = (Z == 1.) ? 202.4 : 182.7 * pow(Z, -1. / 3.);
        const double BZ_e = (Z == 1.) ? 446. : 1429. * pow(Z, -2. / 3.);
        const double D_n = 1.54 * pow(A, 0.27);
        const double dcs_factor = 4.394466E+20 * rem * rem * Z / A;

        const double E = K + mu;
        const double delta_factor = 0.5 * mu * mu / E;
        const double qe_max = E / (1. + 0.5 * mu * mu / (me * E));

        const double nu = q / E;
        const double delta = delta_factor * nu / (1. - nu);
        double Phi_n, Phi_e;
        Phi_n = log(BZ_n * (mu + delta * (D_n * sqrte - 2.)) /
            (D_n * (me + delta * sqrte * BZ_n)));
        if (Phi_n < 0.) Phi_n = 0.;
        if (q < qe_max) {
                Phi_e = log(BZ_e * mu /
                    ((1. + delta * phie_factor) * (me + delta * sqrte * BZ_e)));
                if (Phi_e < 0.) Phi_e = 0.;
        } else
                Phi_e = 0.;

        const double dcs =
            dcs_factor * (Z * Phi_n + Phi_e) * (4. / 3. * (1. / nu - 1.) + nu);
        return (dcs < 0.) ? 0. : dcs;
}

/**
 * The Pair production differential cross section.
 *
 * @param element The target atomic element.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the e+e- pair.
 * @return The differential cross section in m^2/kg
 *
 * The differential cross section is computed following R.P. Kokoulin's
 * formulae taken from the Geant4 Physics Reference Manual.
 */
double dcs_pair_production(
    const struct atomic_element * element, double K, double q)
{
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
        const double Z13 = pow(element->Z, 1. / 3.);
        if (q >= K + s_shared->mass * (1. - 0.75 * sqrte * Z13)) return 0.;

        /*  Precompute some constant factors for the integration. */
        const double nu = q / (K + s_shared->mass);
        const double r = s_shared->mass / ELECTRON_MASS;
        const double beta = 0.5 * nu * nu / (1. - nu);
        const double xi_factor = 0.5 * r * r * beta;
        const double A = (element->Z == 1.) ? 202.4 : 183.;
        const double AZ13 = A / Z13;
        const double cL = 2. * sqrte * ELECTRON_MASS * AZ13;
        const double cLe = 2.25 * Z13 * Z13 / (r * r);

        /*  Compute the bound for the integral. */
        const double gamma = 1. + K / s_shared->mass;
        const double x0 = 4. * ELECTRON_MASS / q;
        const double x1 = 6. / (gamma * (gamma - q / s_shared->mass));
        const double argmin =
            (x0 + 2. * (1. - x0) * x1) / (1. + (1. - x1) * sqrt(1. - x0));
        if ((argmin >= 1.) || (argmin <= 0.)) return 0.;
        const double tmin = log(argmin);

        /*  Compute the integral over t = ln(1-rho). */
        double I = 0.;
        int i;
        for (i = 0; i < 8; i++) {
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
                if (element->Z == 1.) {
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
        const double dcs = 1.0807718E-07 * element->Z * (element->Z + zeta) *
            (K + s_shared->mass - q) * I / (element->A * q);
        return (dcs < 0.) ? 0. : dcs;

#undef N_GQ
}

/** ALLM97 parameterisation of the proton structure function, F2.
 *
 * @param x       The fractional kinetic energy lost to the photon.
 * @param Q2      The negative four momentum squared.
 * @return The corresponding value of the proton structure function, F2.
 *
 * References:
 *      DESY 97-251 [arXiv:hep-ph/9712415].
 */
double dcs_photonuclear_f2_allm(double x, double Q2)
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

        const double M2 = 0.8803505929;
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

/* The F2 structure function for atomic weight A.
 *
 * @param x       The fractional kinetic energy lost to the photon.
 * @param F2p     The proton structure function, F2.
 * @param A       The atomic weight.
 * @return The corresponding value of the structure function, F2.
 *
 * The F2 structure function for a nucleus of atomic weight A is computed
 * according to DRSS, including a Shadowing factor.
 *
 * References:
 *      Dutta et al., Phys.Rev. D63 (2001) 094020 [arXiv:hep-ph/0012350].
 */
double dcs_photonuclear_f2a_drss(double x, double F2p, double A)
{
        double a = 1.0;
        if (x < 0.0014)
                a = exp(-0.1 * log(A));
        else if (x < 0.04)
                a = exp((0.069 * log10(x) + 0.097) * log(A));

        return (0.5 * A * a *
            (2.0 + x * (-1.85 + x * (2.45 + x * (-2.35 + x)))) * F2p);
}

/* The R ratio of longitudinal to transverse structure functions.
 *
 * @param x       The fractional kinetic energy lost to the photon.
 * @param Q2      The negative four momentum squared.
 *
 * References:
 *      Whitlow, SLAC-PUB-5284.
 */
double dcs_photonuclear_r_whitlow(double x, double Q2)
{
        double q2 = Q2;
        if (Q2 < 0.3) q2 = 0.3;

        const double theta =
            1 + 12.0 * q2 / (1.0 + q2) * 0.015625 / (0.015625 + x * x);

        return (0.635 / log(q2 / 0.04) * theta + 0.5747 / q2 -
            0.3534 / (0.09 + q2 * q2));
}

/** The doubly differential cross sections d^2S/(dq*dQ2) for photonuclear
 * interactions.
 *
 * @param A       The target atomic weight.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @param Q2      The negative four momentum squared.
 * @return The doubly differential cross section in m^2/kg/GeV^3.
 *
 * References:
 *      Dutta et al., Phys.Rev. D63 (2001) 094020 [arXiv:hep-ph/0012350].
 */
double dcs_photonuclear_d2(double A, double K, double q, double Q2)
{
        const double ml = s_shared->mass;
        const double cf = 1.5676209E-08 / A;
        const double M = 0.931494;
        const double E = K + ml;

        const double y = q / E;
        const double x = 0.5 * Q2 / (M * q);
        const double F2p = dcs_photonuclear_f2_allm(x, Q2);
        const double F2A = dcs_photonuclear_f2a_drss(x, F2p, A);
        const double R = dcs_photonuclear_r_whitlow(x, Q2);

        const double dds = (1 - y +
                               0.5 * (1 - 2 * ml * ml / Q2) *
                                   (y * y + Q2 / (E * E)) / (1 + R)) /
                (Q2 * Q2) -
            0.25 / (E * E * Q2);

        return cf * F2A * dds / q;
}

/**
 * The Photonuclear differential cross section.
 *
 * @param element The target atomic element.
 * @param K       The projectile initial kinetic energy.
 * @param q       The kinetic energy lost to the photon.
 * @return The differential cross section in m^2/kg.
 *
 * The differential cross-section is computed following DRSS, with ALLM97
 * parameterisation of the structure function F2.
 *
 * References:
 *      Dutta et al., Phys.Rev. D63 (2001) 094020 [arXiv:hep-ph/0012350].
 */
double dcs_photonuclear(
    const struct atomic_element * element, double K, double q)
{
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

        /* Unpack and check the kinematic. */
        if (dcs_photonuclear_check(K, q)) return 0.;

        const double A = element->A;
        const double ml = s_shared->mass;
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
                ds += dcs_photonuclear_d2(A, K, q, Q2) * Q2 * wGQ[i];
        }

        if (ds < 0.) ds = 0.;
        return 0.5 * ds * dpQ2 * E;

#undef N_GQ
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
double dcs_ionisation(const struct atomic_element * element, double K, double q)
{
        const double Z = element->Z;
        const double A = element->A;
        const double P2 = K * (K + 2. * s_shared->mass);
        const double E = K + s_shared->mass;
        const double Wmax = 2. * ELECTRON_MASS * P2 /
            (s_shared->mass * s_shared->mass +
                ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
        if ((Wmax < X_FRACTION * K) || (q > Wmax))
                return 0.;
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
        const double m1 = s_shared->mass - ELECTRON_MASS;
        if (K >= 0.5 * m1 * m1 / ELECTRON_MASS) {
                const double L1 = log(1. + 2. * q / ELECTRON_MASS);
                Delta = 1.16141E-03 * L1 * (log(4. * E * (E - q) /
                    (s_shared->mass * s_shared->mass)) - L1);
        }

        return cs * (1. + Delta);
}

/**
 * The analytical form for the partial integral of the ionisation DCS.
 *
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
double dcs_ionisation_integrate(
    int mode, const struct atomic_element * element, double K, double xlow)
{
        const double Z = element->Z;
        const double A = element->A;
        const double P2 = K * (K + 2. * s_shared->mass);
        const double E = K + s_shared->mass;
        const double Wmax = 2. * ELECTRON_MASS * P2 /
            (s_shared->mass * s_shared->mass +
                ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
        if (Wmax < X_FRACTION * K)
                return 0.;
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

double dcs_ionisation_randomise(struct pumas_context * context,
    const struct atomic_element * element, double K, double xlow)
{
        const double P2 = K * (K + 2. * s_shared->mass);
        const double E = K + s_shared->mass;
        const double Wmax = 2. * ELECTRON_MASS * P2 /
            (s_shared->mass * s_shared->mass +
                ELECTRON_MASS * (ELECTRON_MASS + 2. * E));
        if (Wmax < X_FRACTION * K)
                return K;
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
                        q = Wmin / (1. - context->random(context) *
                            (1. - Wmin / Wmax));
                }

                /* Rejection sampling. */
                const double r0 = a0 + a2 / (q * q);
                const double r1 = r0 + a1 / q;
                if (context->random(context) * r0 < r1)
                        break;
        }

        return K - q;
}

/**
 * Encapsulation for the evaluation of DCS.
 *
 * @param context  The simulation context.
 * @param dcs_func The DCS function to evaluate.
 * @param element  The target atomic element.
 * @param K        The projectile initial kinetic energy.
 * @param q        The transfered energy.
 * @return The DCS value or `0`.
 *
 * This routine encapsulate the evluation of DCS during the MC. It takes care of
 * checking whether an approximate model can be used or not. In addition it
 * applies a Jacobian weight factor for changing from nu = q / E to x = q / K.
 */
double dcs_evaluate(struct pumas_context * context, dcs_function_t * dcs_func,
    const struct atomic_element * element, double K, double q)
{
        /* Compute the Jacobian factor. */
        const double wj = K / (K + s_shared->mass);

        /* Check if the process has a valid tabulated model. */
        if (dcs_func == dcs_ionisation)
                return dcs_func(element, K, q) * wj;
        else if ((dcs_func == dcs_photonuclear) && dcs_photonuclear_check(K, q))
                return 0.;

        /* Check if the exact computation should be used. */
        const double min_k = s_shared ?
            *table_get_K(s_shared->dcs_model_offset) :
            DCS_MODEL_MIN_KINETIC;
        if (!s_shared || (K <= min_k) || (q < X_FRACTION * K) ||
            (q > DCS_MODEL_MAX_FRACTION * K))
                return dcs_func(element, K, q) * wj;

        /* Get the coefficients of the polynomial approximation. */
        const int index = dcs_get_index(dcs_func);
        const int offset = s_shared->dcs_model_offset;
        const int n = DCS_MODEL_ORDER_P + DCS_MODEL_ORDER_Q + 1;
        const float * const coeff = table_get_dcs_coeff(element, index, 0);
        const int imax = s_shared->n_kinetics - 1;
        const float *c0, *c1;
        int i0;
        float h1;
        const double Kmax = *table_get_K(imax);
        if ((K >= Kmax) ||
            ((i0 = table_index(context, table_get_K(0), K)) >= imax)) {
                /*  Rescale and use the last tabulated value. */
                h1 = 0.;
                c0 = c1 = coeff + (imax - offset) * (n + DCS_SAMPLING_N);
        } else {
                const float K0 = (float)(*table_get_K(i0));
                const float K1 = (float)(*table_get_K(i0 + 1));
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
 * @param context The simulation context.
 * @param ki      The initial kinetic energy.
 * @param kf      The final kinetic energy.
 * @return The cosine of the polar angle.
 */
double polar_bremsstrahlung(
    struct pumas_context * context, double ki, double kf)
{
        const double e = ki + s_shared->mass;
        const double q = ki - kf;
        double rmax2 = e / q - 1.;
        rmax2 = (rmax2 > 1.) ? 1. : rmax2 * rmax2;
        double r = context->random(context) * rmax2 / (1. + rmax2);
        r = sqrt(r / (1. - r));

        return cos(r * s_shared->mass * q / (e * (kf + s_shared->mass)));
}

/**
 * Sample the polar angle in a Pair Production event.
 *
 * @param context The simulation context.
 * @param ki      The initial kinetic energy.
 * @param kf      The final kinetic energy.
 * @return The cosine of the polar angle.
 *
 * The polar angle is sampled assuming a virtual Bremsstrahlung event.
 */
double polar_pair_production(
    struct pumas_context * context, double ki, double kf)
{
        return polar_bremsstrahlung(context, ki, kf);
}

/**
 * Sample the polar angle in a photonuclear event.
 *
 * @param context The simulation context.
 * @param ki      The initial kinetic energy.
 * @param kf      The final kinetic energy.
 * @return The cosine of the polar angle.
 */
double polar_photonuclear(struct pumas_context * context, double ki, double kf)
{
        const double q = ki - kf;
        const double e = ki + s_shared->mass;
        const double y = q / e;
        const double tmin = s_shared->mass * s_shared->mass * y * y / (1. - y);
        const double tmax = 1.87914 * q;
        const double m02 = 0.4;
        const double q2 = q * q;
        const double t1 = (q2 < m02) ? q2 : m02;

        const double p = context->random(context);
        const double r = tmax * (tmin + t1) / (tmin * (tmax + t1));
        const double tp = tmax * t1 / ((tmax + t1) * pow(r, p) - tmax);
        const double ct = 1. -
            (tp - tmin) /
                (2. *
                        (e * (kf + s_shared->mass) -
                            s_shared->mass * s_shared->mass) -
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
 * @param context The simulation context.
 * @param ki      The initial kinetic energy.
 * @param kf      The final kinetic energy.
 * @return The cosine of the polar angle.
 *
 * The polar angle is set from energy-momentum conservation assuming that the
 * electron is initially at rest. See for example appendix A of Fernandez-Varea
 * et al., NIMB 229 (2005) 185-218.
 */
double polar_ionisation(struct pumas_context * context, double ki, double kf)
{
        if (kf > ki) {
                const double tmp = ki;
                ki = kf, kf = tmp;
        }
        const double p02 = ki * (ki + 2. * s_shared->mass);
        const double p12 = kf * (kf + 2. * s_shared->mass);
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
 * @param f      The objective function to resolve.
 * @param xa     The lower bound of the search interval.
 * @param xb     The upper bound of the search interval.
 * @param fa_p   The initial value at *a* if already computed.
 * @param fb_p   The initial value at *b* if already computed.
 * @param xtol   The absolute tolerance on the root value.
 * @param rtol   The absolute tolerance on the root value.
 * @param params A handle for passing additional parameters to the
 *               objective function.
 * @param x0     An estimate of the root over `[xa; xb]`.
 * @return On success `0` is returned. Otherwise a negative number.
 *
 * The root is searched for over `[xa; xb]` using Ridder's method
 * (https://en.wikipedia.org/wiki/Ridders%27_method). If initial values of the
 * objective have already been computed they can be passed over as *fa_p* or
 * *fb_p*. Otherwise these values must be `NULL`. The relative tolerance is
 * relative to min(xa, xb).
 */
int math_find_root(double (*f)(double x, void * params), double xa, double xb,
    const double * fa_p, const double * fb_p, double xtol, double rtol,
    int max_iter, void * params, double * x0)
{
        /*  Check the initial values. */
        double fa, fb;
        if (fa_p == NULL)
                fa = (*f)(xa, params);
        else
                fa = *fa_p;
        if (fb_p == NULL)
                fb = (*f)(xb, params);
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
                const double fm = (*f)(xm, params);
                double sgn = (fb > fa) ? 1. : -1.;
                double dn = sgn * dm * fm / sqrt(fm * fm - fa * fb);
                sgn = (dn > 0.) ? 1. : -1.;
                dn = fabs(dn);
                dm = fabs(dm) - 0.5 * tol;
                if (dn < dm) dm = dn;
                xn = xm - sgn * dm;
                const double fn = (*f)(xn, params);
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
                        for (k = l; k < n; k++) { scale += fabs(ai[k]); }
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
                                for (k = l; k < n; k++) { work[k] = ai[k] / h; }
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
                                for (k = l; k < n; k++) { ai[k] *= scale; }
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
                        for (j = l; j < n; j++) { ai[j] = 0.0; }
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
                                if (an + fabs(w[l1]) == an) { break; }
                        }

                        /* Cancellation of work[l] if l>0. */
                        if (cancel) {
                                c = 0.0;
                                s = 1.0;
                                for (i = l; i <= k; i++) {
                                        f = s * work[i];
                                        if (an + fabs(f) == an) { break; }
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
                                if (its >= ITMAX) { jstat = k + 1; }

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
                for (jj = 0; jj < n; jj++) { s += vj[jj] * work[jj]; }
                x[j] = s;
        }
}

#ifdef _BUILD_TABULATE
/*
 * Low level routines for computing energy loss tabulations.
 */
#include "pumas-tabulate.h"

/* Dry mode initialisation, without loading the energy loss tables. */
enum pumas_return _pumas_initialise_dry(enum pumas_particle particle,
    const char * mdf_path, struct pumas_error * error)
{
        return _initialise(particle, mdf_path, NULL, error, 1);
}

/* The density effect for the electronic energy loss. */
static void electronic_density_effect(
    struct tabulation_material * m, double kinetic)
{
        const double c = 2. * log(10.);
        const double r = kinetic / s_shared->mass;
        const double x = 0.5 * log10(r * (r + 2.));
        if (x < m->x_0)
                m->delta = m->delta0 > 0. ?
                    m->delta0 * pow(10., 2. * (x - m->x_0)) :
                    0.;
        else if (x < m->x_1)
                m->delta = c * x - m->Cbar + m->a * pow(m->x_1 - x, m->k);
        else
                m->delta = c * x - m->Cbar;
}

/* The average energy loss from atomic electrons. */
static double electronic_energy_loss(
    struct tabulation_material * m, double kinetic)
{
        /* Kinematic factors. */
        const double E = kinetic + s_shared->mass;
        const double P2 = kinetic * (kinetic + 2. * s_shared->mass);
        const double beta2 = P2 / (E * E);

        /* Electronic Bremsstrahlung correction. */
        const double r = ELECTRON_MASS / s_shared->mass;
        const double Qmax =
            2. * r * P2 / (s_shared->mass * (1. + r * r) + 2. * r * E);
        const double lQ = log(2. * Qmax / ELECTRON_MASS);
        const double Delta =
            5.8070487E-04 * (log(2. * E / s_shared->mass) - lQ / 3.) * lQ * lQ;

        /* Density effect. */
        electronic_density_effect(m, kinetic);

        /* Bethe-Bloch equation. */
        return 0.307075E-04 * s_shared->material_ZoA[m->index] *
            (0.5 / beta2 *
                    (log(2. * ELECTRON_MASS * P2 * Qmax /
                         (s_shared->mass * s_shared->mass * m->I * m->I)) -
                        m->delta) -
                1. + 0.125 * Qmax * Qmax / P2 + Delta);
}

static void tabulate_element(
    struct tabulation_element * data, int n_kinetics, double * kinetic)
{
        const struct atomic_element * element = s_shared->element[data->index];

        /* Loop over the kinetic energy values. */
        double * v;
        int ik, ip;
        for (ik = 0, v = data->data; ik < n_kinetics; ik++) {
                const double k = kinetic[ik];
                double x = 1E-06 / k;
                if (x < 1E-05) x = 1E-05;
                const int n = (int)(-1E+02 * log10(x)); /* 100 pts per decade. */
                for (ip = 0; ip < N_DEL_PROCESSES - 1; ip++, v++)
                        *v = compute_dcs_integral(
                            1, element, k, dcs_get(ip), x, n);
        }
}

/*
 * Create a new energy loss table for an element and add it to the stack of
 * temporary data.
 */
static struct tabulation_element * tabulation_element_create(
    struct tabulation_data * data, int element)
{
        /* Allocate memory for the new element. */
        struct tabulation_element * e =
            allocate(sizeof(*e) + 3 * data->n_kinetics * sizeof(*e->data));
        if (e == NULL) return NULL;
        e->index = element;
        e->fraction = 0.;

        /* Add the element's data on top of the stack. */
        if (data->e_stack != NULL) data->e_stack->next = e;
        e->prev = data->e_stack;
        e->next = NULL;
        data->e_stack = e;

        return e;
}

/*
 * Get the energy loss table for an element from the temporary data and
 * put it on top of the stack.
 */
static struct tabulation_element * tabulation_element_get(
    struct tabulation_data * data, int element)
{
        struct tabulation_element * e;
        for (e = data->e_stack; e != NULL; e = e->next) {
                if (e->index == element) {
                        struct tabulation_element * next = e->next;
                        if (next != NULL) {
                                /* Put the element on top of the stack. */
                                struct tabulation_element * prev = e->prev;
                                if (prev != NULL) prev->next = next;
                                next->prev = prev;
                                data->e_stack->next = e;
                                e->prev = data->e_stack;
                                data->e_stack = e;
                        }
                        return e;
                }
        }

        /* The element wasn't found, return `NULL`. */
        return NULL;
}

/**
 * Tabulate the energy loss for the material on top of the temporary
 * material's stack.
 *
 * @param data    Handle for the temporary data.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * On success the corresponding material data are dropped and the stack is
 * updated.
 *
 * __Warnings__
 *
 * This function is **not** thread safe.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_IO_ERROR        The output file already exists.
 *
 *     PUMAS_RETURN_MEMORY_ERROR    Some memory couldn't be allocated.
 *
 *     PUMAS_RETURN_PATH_ERROR      The output file could not be created.
 */
enum pumas_return _pumas_tabulate(struct tabulation_data * data)
{
        struct tabulation_material * m = data->m_stack;

        /* Compute the mean excitation energy, if not provided. */
        if (m->I <= 0.) {
                double lnI = 0., Z = 0.;
                struct material_component * component;
                int iel;
                for (iel = 0, component = s_shared->composition[m->index];
                     iel < s_shared->elements_in[m->index];
                     iel++, component++) {
                        struct atomic_element * e =
                            s_shared->element[component->element];
                        const double nZ = component->fraction * e->Z / e->A;
                        lnI += nZ * log(e->I);
                        Z += nZ;
                }
                m->I = exp(lnI / Z);
        }

        /* Compute the density effect coefficients, if not provided. */
        if (m->a <= 0.) {
                /* Use the Sternheimer and Peierls recipee. */
                m->k = 3.;
                const double hwp = 28.816E-09 *
                    sqrt(m->density * 1E-03 * s_shared->material_ZoA[m->index]);
                m->Cbar = 2. * log(m->I / hwp);
                if (m->state == TABULATION_STATE_GAZ) {
                        if (m->Cbar < 10.) {
                                m->x_0 = 1.6, m->x_1 = 4.;
                        } else if (m->Cbar < 10.5) {
                                m->x_0 = 1.7, m->x_1 = 4.;
                        } else if (m->Cbar < 11.) {
                                m->x_0 = 1.8, m->x_1 = 4.;
                        } else if (m->Cbar < 11.5) {
                                m->x_0 = 1.9, m->x_1 = 4.;
                        } else if (m->Cbar < 12.25) {
                                m->x_0 = 2.0, m->x_1 = 4.;
                        } else if (m->Cbar < 13.804) {
                                m->x_0 = 2.0, m->x_1 = 5.;
                        } else {
                                m->x_0 = 0.326 * m->Cbar - 1.5;
                                m->x_1 = 5.0;
                        }
                } else {
                        if (m->I < 100.E-09) {
                                m->x_1 = 2.;
                                if (m->Cbar < 3.681)
                                        m->x_0 = 0.2;
                                else
                                        m->x_0 = 0.326 * m->Cbar - 1.;
                        } else {
                                m->x_1 = 3.;
                                if (m->Cbar < 5.215)
                                        m->x_0 = 0.2;
                                else
                                        m->x_0 = 0.326 * m->Cbar - 1.5;
                        }
                }
                const double dx = m->x_1 - m->x_0;
                m->a = (m->Cbar - 2. * log(10.) * m->x_0) / (dx * dx * dx);
        }

        /*
         * Tabulate the radiative energy losses for the constitutive atomic
         * elements, if not already done. Note that the element are also sorted
         * on top of the *data* stack.
         */
        struct material_component * component;
        int iel;
        for (iel = 0, component = s_shared->composition[m->index];
             iel < s_shared->elements_in[m->index]; iel++, component++) {
                /*
                 * Get the requested element's data and put them on top of
                 * the stack for further usage.
                 */
                struct tabulation_element * e =
                    tabulation_element_get(data, component->element);

                if (e == NULL) {
                        /* Create and tabulate the new element. */
                        e = tabulation_element_create(data, component->element);
                        if (e == NULL) return PUMAS_RETURN_MEMORY_ERROR;
                        tabulate_element(e, data->n_kinetics, data->kinetic);
                }

                /* Set the fraction in the current material. */
                e->fraction = component->fraction;
        }

        /* Check and open the output file. */
        int n_d = (data->outdir == NULL) ? 0 : strlen(data->outdir) + 1;
        int n_f = strlen(s_shared->dedx_filename[m->index]) + 1;
        char * path = reallocate(data->path, n_d + n_f);
        if (path == NULL) return PUMAS_RETURN_MEMORY_ERROR;
        data->path = path;
        if (n_d) {
                memcpy(path, data->outdir, n_d);
                path[n_d - 1] = '/';
        }
        memcpy(path + n_d, s_shared->dedx_filename[m->index], n_f);

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
            (s_shared->particle == PUMAS_PARTICLE_MUON) ? "Muon" : "Tau";
        fprintf(stream, " Incident particle is a %s with M = %.5lf MeV\n", type,
            (double)(s_shared->mass * 1E+03));
        fprintf(stream, " Index = %d: %s\n", m->index,
            s_shared->material_name[m->index]);
        fprintf(stream, "      Absorber with <Z/A> = %.5lf\n",
            s_shared->material_ZoA[m->index]);
        fprintf(stream,
            " Sternheimer coef:  a     k=m_s   x_0    x_1    "
            "I[eV]   Cbar  delta0\n");
        fprintf(stream,
            "                %7.4lf %7.4lf %7.4lf %7.4lf %6.1lf "
            "%7.4lf %.2lf\n",
            m->a, m->k, m->x_0, m->x_1, m->I * 1E+09, m->Cbar, m->delta0);
        fprintf(stream, "\n *** Table generated with PUMAS v%d.%d ***\n\n",
            PUMAS_VERSION, PUMAS_SUBVERSION);
        fprintf(stream,
            "      T         p     Ionization  brems     pair     "
            "photonuc  Radloss    dE/dx   CSDA Range  delta   beta\n");
        fprintf(stream,
            "    [MeV]    [MeV/c]  -----------------------"
            "[MeV cm^2/g]------------------------  [g/cm^2]\n");

        /* Loop on the kinetic energy values and print the table. */
        double X = 0., dedx_last = 0.;
        int i;
        for (i = 0; i < data->n_kinetics; i++) {
                /* Compute the electronic energy loss. */
                double elec = electronic_energy_loss(m, data->kinetic[i]);

                double brad[N_DEL_PROCESSES - 1];
                memset(brad, 0x0, sizeof(brad));
                struct tabulation_element * e;
                int iel;
                for (iel = 0, e = data->e_stack;
                     iel < s_shared->elements_in[m->index];
                     iel++, e = e->prev) {
                        int j;
                        for (j = 0; j < N_DEL_PROCESSES - 1; j++)
                                brad[j] += e->fraction *
                                    e->data[(N_DEL_PROCESSES - 1) * i + j];
                }

                /* Update the CDSA range. */
                const double radloss = brad[0] + brad[1] + brad[2];
                const double dedx = radloss + elec;
                if (i == 0)
                        X = 0.5 * data->kinetic[i] / dedx;
                else
                        X += 0.5 * (data->kinetic[i] - data->kinetic[i - 1]) *
                            (1. / dedx + 1. / dedx_last);
                dedx_last = dedx;

                const double p = sqrt(data->kinetic[i] *
                    (data->kinetic[i] + 2. * s_shared->mass));
                const double beta = p / (data->kinetic[i] + s_shared->mass);
                const double MeV = 1E+03;
                const double cmgs = 1E+04;
                fprintf(stream,
                    "  %.3lE %.3lE %.3lE %.3lE %.3lE %.3lE "
                    "%.3lE %.3lE %.3lE %7.4lf %7.5lf\n",
                    data->kinetic[i] * MeV, p * MeV, elec * cmgs,
                    brad[0] * cmgs, brad[1] * cmgs, brad[2] * cmgs,
                    radloss * cmgs, dedx * cmgs, X * MeV / cmgs, m->delta,
                    beta);
        }

        /* Close and return. */
        fclose(stream);
        return PUMAS_RETURN_SUCCESS;
}
#endif
