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

#ifndef pumas_h
#define pumas_h
#ifdef __cplusplus
extern "C" {
#endif

#ifndef PUMAS_API
#define PUMAS_API
#endif

/* For C standard streams. */
#ifndef FILE
#include <stdio.h>
#endif

/**
 * Particle types supported by the PUMAS transport engine.
 */
enum pumas_particle {
        /** The muon or anti-muon lepton. */
        PUMAS_PARTICLE_MUON = 0,
        /** The tau or anti-tau lepton. */
        PUMAS_PARTICLE_TAU
};

/**
 * Keys for some of the tabulated properties used by PUMAS.
 */
enum pumas_property {
        /** The macroscopic inelastic cross-section, in m^(2)/kg. */
        PUMAS_PROPERTY_CROSS_SECTION = 0,
        /** The average energy loss, in GeV/(kg/m^(2)). */
        PUMAS_PROPERTY_ENERGY_LOSS,
        /** The particle range, in kg/m^(2). */
        PUMAS_PROPERTY_GRAMMAGE,
        /** The particle kinetic energy, in GeV. */
        PUMAS_PROPERTY_KINETIC_ENERGY,
        /** The total magnetic rotation angle, in radians/(kg/m^(3)). */
        PUMAS_PROPERTY_MAGNETIC_ROTATION,
        /** The particle proper time, in kg/m^(2). */
        PUMAS_PROPERTY_PROPER_TIME,
        /** The macroscopic elastic scattering 1^(st) path length, in kg/m^(2).
         */
        PUMAS_PROPERTY_SCATTERING_LENGTH
};

/**
 * Schemes for the computation of energy losses.
 */
enum pumas_scheme {
        /** All energy losses are disabled.
         *
         * **Note** : This mode is provided for test purpose only. Running
         * without energy losses requires specifying a range or grammage limit.
         */
        PUMAS_SCHEME_NO_LOSS = -1,
        /** Energy losses are purely determinstic as given by the Continuously
         * Slowing Down Approximation (CSDA).
         */
        PUMAS_SCHEME_CSDA,
        /** Energy losses are simulated using an hybrid Monte-Carlo scheme with
         * hard stochastic interactions and a deterministic Continuous Energy
         * Loss (CEL).
         */
        PUMAS_SCHEME_HYBRID,
        /** In addition to the hybrid scheme small fluctuations of the CEL are
         * also simulated.
         */
        PUMAS_SCHEME_DETAILED
};

/**
 * Modes for handling particle decays.
 */
enum pumas_decay {
        /** Decays are disabled. */
        PUMAS_DECAY_NONE = 0,
        /**
         * Decays are accounted for by a weight factor. This is efficient
         * for muons but irrelevant -numerically instable- for the forward
         * transport of taus since they decay in flight. Hence it is
         * disabled in the latter case.
         */
        PUMAS_DECAY_WEIGHT,
        /** Decays are enabled as a specific Monte-Carlo process. */
        PUMAS_DECAY_PROCESS
};

/** Return codes for the API functions. */
enum pumas_return {
        /** Execution was successful. */
        PUMAS_RETURN_SUCCESS = 0,
        /** End of file was reached. */
        PUMAS_RETURN_END_OF_FILE,
        /** The specified decay mode is not valid. */
        PUMAS_RETURN_DECAY_ERROR,
        /** Some medium has a wrong density value. */
        PUMAS_RETURN_DENSITY_ERROR,
        /** Some data file is not complete. */
        PUMAS_RETURN_INCOMPLETE_FILE,
        /** Some index is out of validity range. */
        PUMAS_RETURN_INDEX_ERROR,
        /** The library is not/already initialised. */
        PUMAS_RETURN_INITIALISATION_ERROR,
        /** An internal library error occured. */
        PUMAS_RETURN_INTERNAL_ERROR,
        /** Some read /write error occured. */
        PUMAS_RETURN_IO_ERROR,
        /** Some file is badly formated. */
        PUMAS_RETURN_FORMAT_ERROR,
        /** Wrong propagation medium. */
        PUMAS_RETURN_MEDIUM_ERROR,
        /** Some memory couldn't be allocated. */
        PUMAS_RETURN_MEMORY_ERROR,
        /** A user supplied limit is required. */
        PUMAS_RETURN_MISSING_LIMIT,
        /** The random callback is not defined. */
        PUMAS_RETURN_MISSING_RANDOM,
        /** Some file couldn't be found. */
        PUMAS_RETURN_PATH_ERROR,
        /** A raise was called without any catch. */
        PUMAS_RETURN_RAISE_ERROR,
        /** Some input string is too long. */
        PUMAS_RETURN_TOO_LONG,
        /** No energy loss path specified. */
        PUMAS_RETURN_UNDEFINED_DEDX,
        /** No MDF file specified. */
        PUMAS_RETURN_UNDEFINED_MDF,
        /** An unkwon element was specified. */
        PUMAS_RETURN_UNKNOWN_ELEMENT,
        /** An unkwon material was specified. */
        PUMAS_RETURN_UNKNOWN_MATERIAL,
        /** The particle type is not known. */
        PUMAS_RETURN_UNKNOWN_PARTICLE,
        /** Some input value is not valid. */
        PUMAS_RETURN_VALUE_ERROR,
        /** The number of PUMAS return codes.  */
        PUMAS_N_RETURNS
};

/**
 * Container for a Monte-Carlo state.
 */
struct pumas_state {
        /** The particle's electric charge. Note that non physical values,
         * i.e. different from 1 or -1, could be set. */
        double charge;
        /** The current kinetic energy, in GeV. */
        double kinetic;
        /** The total travelled distance, in m. */
        double distance;
        /** The total travelled grammage, in kg/m^2. */
        double grammage;
        /** The particle's proper time, in m/c. */
        double time;
        /** The Monte-Carlo weight. */
        double weight;
        /** The absolute location, in m. */
        double position[3];
        /** The momentum's direction. */
        double direction[3];
        /** Status flag telling if the particle has decayed or not.  */
        int decayed;
};

/**
 * The local properties of a propagation medium.
 */
struct pumas_locals {
        /** The material local density, in kg/m^3. */
        double density;
        /** The local magnetic field components, in T. */
        double magnet[3];
};

struct pumas_medium;
/**
 * Callback for setting the local properties of a propagation medium.
 *
 * @param medium    The propagation medium.
 * @param state     The Monte-Carlo state for which the local properties are
 *                  requested.
 * @param locals    A pointer to a `pumas_locals` structure to update.
 * @return A local stepping limit.
 *
 * The callback must return a proposed Monte-Carlo stepping distance, in m,
 * consistent with the size of the propagation medium inhomogeneities,
 * e. g. 1 % of &rho; / |&nabla; &rho;|. Note that returning zero or less
 * signs that the propagation medium is fully uniform.
 *
 * **Warning** : it is an error to return zero or less for any position of the
 * medium if at least one area is not uniform. Instead one should use two
 * different media even though they have the same material base.
 *
 */
typedef double(pumas_locals_cb)(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals);

/**
 * Description of a propagation medium.
 *
 * A propagation medium is fully defined by:
 *
 * - a `material` composition with a uniform relative content.
 * - `pumas_locals` properties set by a user provided `pumas_locals_cb`
 * callback.
 */
struct pumas_medium {
        /**
         * The material index in the Material Description File (MDF). It can
         * be mapped to the corresponding name with the `pumas_material_`
         * functions.
         */
        int material;
        /**
         * The user supplied callback for setting the medium local properties.
         */
        pumas_locals_cb * locals;
};

/** A handle to a recorded Monte-Carlo frame.
 *
 *  This structure exposes data relative to a recorded frame. It is not meant
 * to be modified by the user.
 */
struct pumas_frame {
        /** The recorded state. */
        struct pumas_state state;
        /** The corresponding propagation medium. */
        struct pumas_medium * medium;
        /** Link to the next frame in the record. */
        struct pumas_frame * next;
};

/**
 * A handle for recording Monte-Carlo frames.
 *
 * This structure is a proxy for recording Monte-Carlo states and accessing
 * them. Although it exposes some public data that the user may alter it also
 * encloses other opaque data. Therefore, it **must** be handled with the
 * `pumas_recorder` functions.
 *
 * **Note** : in order to enable or disable a recorder it is enough to
 * link or unlink it from the `recorder` field of any `pumas_context`. Only the
 * corresponding context will be recorded.
 */
struct pumas_recorder {
        /** Link to the 1^(st) recorded frame or `NULL` if none. This field
         * should not be modified.
         */
        struct pumas_frame * first;
        /** The total number of recorded frames. This field should not be
         * modified.
         */
        int length;
        /**
         * The sampling period of the recorder, if strictly positive. If set to
         * zero recording is disabled. Otherwise, setting a negative value
         * enables a test mode where the full detail of the Monte-Carlo
         * stepping is recorded.
         */
        int period;
};

struct pumas_context;
/**
 * Callback for locating the propagation medium of a `pumas_state`.
 *
 * @param context   The Monte-Carlo context requiring a medium.
 * @param state     The Monte-Carlo state for which the medium is requested.
 * @param medium    A pointer to the corresponding medium or `NULL` if the state
 * has exit the simulation area.
 * @return The proposed step size or `0` for an infinite medium.
 *
 * The callback must propose a Monte-Carlo stepping distance, in m,
 * consistent with the geometry. Note that returning zero or less
 * signs that the corresponding medium has no boundaries.
 *
 * **Warning** : it is an error to return zero or less for any state if the
 * extension is finite.
 */
typedef double pumas_medium_cb(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium);

/**
 * Generic function pointer.
 *
 * This is a generic function pointer used to identify the library functions,
 * e.g. for error handling.
 */
typedef void pumas_function_t(void);

/**
 * A container for additional info on errors during initialisation.
 *
 * This structure is a container for getting additional info when an error
 * occurs at initialisation, e.g. when parsing the MDF.
 */
struct pumas_error {
        /** The faulty file if any. */
        const char * file;
        /** The faulty line if a parsing error occurs. */
        int line;
};

/**
 * Callback for error handling.
 *
 * @param rc        The PUMAS return code.
 * @param caller    The caller function where the error occured.
 * @param error     Additional info in case of an initialisation error.
 *
 * The user can provide its own error handler. It will be called at the
 * return of any PUMAS library function providing an error code.
 */
typedef void pumas_handler_cb(enum pumas_return rc, pumas_function_t * caller,
    struct pumas_error * error);

/**
 * Callback providing a stream of pseudo random numbers.
 *
 * @param context The simulation context requiring a random number.
 * @return A uniform pseudo random number in [0;1].
 *
 * **Note** : this is the only random stream used by PUMAS. The user must unsure
 * proper behaviour, i.e. that a flat distribution in [0;1] is indeed returned.
 *
 * **Warning** : if multiple contexts are used the user must ensure that this
 * callback is thread safe, e.g. by using independant streams for each context
 * or a locking mechanism in order to share a single random stream.
 */
typedef double(pumas_random_cb)(struct pumas_context * context);

/**
 * A handle for a simulation stream.
 *
 * This structure is a proxy to thread specific data for a simulation stream.
 * It exposes some public data that the user may configure or alter directly.
 * However, it also encloses other opaque data. Therefore, it **must** be
 * initialised and released with the `pumas_context` functions.
 *
 * + The `medium` field must be set after any initialisation with
 * `pumas_context_create` and prior to any call to `pumas_transport`.
 *
 * + Depending on the level of detail of the simulation a random stream must
 * be provided by the user before any call to `pumas_transport`.
 *
 * + For `kinetic_limit`, `distance_max` or `grammage_max` a strictly positive
 * value activates the corresponding limitation. Setting it to `0` or less
 * disables it however.
 */
struct pumas_context {
        /** A medium callback. */
        pumas_medium_cb * medium;
        /** The pseudo random generator callback. */
        pumas_random_cb * random;
        /** A `pumas_frame` recorder. */
        struct pumas_recorder * recorder;
        /** A pointer to additional memory, if any is requested at
         * initialisation.
         */
        void * user_data;

        /**
         * Flag to enable or disable transverse transport. Default is `0`,
         * i.e. transverse transport is enabled.
         */
        int longitudinal;
        /**
         * Flag to switch between forward and reverse transport. Default is
         * `1`, i.e. forward Monte-Carlo.
         */
        int forward;
        /**
         * The scheme used for the computation of energy losses. Default is
         * `PUMAS_SCHEME_DETAILED`.
         */
        enum pumas_scheme scheme;
        /**
         * The mode for handling decays. Default is `PUMAS_DECAY_WEIGHT` for
         * a muon or `PUMAS_DECAY_PROCESS` for a tau.
         */
        enum pumas_decay decay;

        /** The minimum kinetic energy for forward transport, or maximum one
         * for backward transport.
         */
        double kinetic_limit;
        /** The maximum travelled distance. */
        double distance_max;
        /** The maximum travelled grammage. */
        double grammage_max;
        /** The maximum travelled proper time. */
        double time_max;
};

/**
 * Initialise the PUMAS library.
 *
 * @param particle     The type of the particle to transport.
 * @param mdf_path     The path to a Material Description File (MDF), or `NULL`.
 * @param dedx_path    The path to the energy loss tabulation(s), or `NULL`.
 * @param error        A structure to fill with parsing errors, or `NULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Initialise the library from a MDF and a set of energy loss tabulations. Load
 * the materials data and precompute various properties. *mdf_path* and/or
 * *dedx_path* can be set to `NULL`. If so the corresponding path is read from
 * the `PUMAS_MDF` or `PUMAS_DEDX` environment variable.
 *
 * Call `pumas_finalise` in order to unload the library and release all
 * alocated memory.
 *
 * **Warnings** : this function is not thread safe. Trying to (re-)initialise an
 * already initialised library will generate an error. `pumas_finalise` must
 * be called first.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_END_OF_FILE             And unexpected EOF occured.
 *
 *     PUMAS_RETURN_FORMAT_ERROR            A file has a wrong format.
 *
 *     PUMAS_RETURN_INCOMPLETE_FILE         There are missing entries in
 * the MDF.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library is already initialised.
 *
 *     PUMAS_RETURN_IO_ERROR                A file couldn't be read.
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 *
 *     PUMAS_RETURN_PATH_ERROR              A file couldn't be opened.
 *
 *     PUMAS_RETURN_TOO_LONG                Some XML node in the MDF is
 * too long.
 *
 *     PUMAS_RETURN_UNDEFINED_MDF           No energy loss path was provided.
 *
 *     PUMAS_RETURN_UNDEFINED_MDF           No MDF was provided.
 *
 *     PUMAS_RETURN_UNKNOWN_ELEMENT         An elemnt in the MDF wasn't defined.
 *
 *     PUMAS_RETURN_UNKNOWN_MATERIAL        An material in the MDF wasn't
 * defined.
 *
 *     PUMAS_RETURN_UNKNOWN_PARTICLE        The given type is not supported.
 */
PUMAS_API enum pumas_return pumas_initialise(enum pumas_particle particle,
    const char * mdf_path, const char * dedx_path, struct pumas_error * error);

/**
 * Finalise the PUMAS library.
 *
 * Finalise the library and free the shared memory. Call `pumas_initialise` in
 * order to reload the library.
 *
 * **Warnings** : This function is not thread safe. Finalising the library
 * doesn't release the memory allocated for any `pumas_context`.
 */
PUMAS_API void pumas_finalise(void);

/**
 * Dump the PUMAS library configuration to a stream.
 *
 * @param stream    The stream where to dump.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Dump the current library configuration to a stream as a binary object. Note
 * that only globally shared data are dumped, i.e. material properties and
 * tables as read from a MDF. Simulation contexts, media, recorders, ect ...
 * are not. This binary dump allows for a fast initialisation of the library
 * in subsequent uses.
 *
 * **Warnings** : The binary dump is raw formated, hence *a priori* platform
 * dependent.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 *
 *     PUMAS_RETURN_IO_ERROR                Couldn't write to the stream.
 */
PUMAS_API enum pumas_return pumas_dump(FILE * stream);

/**
 * Load the configuration from a binary dump.
 *
 * @param stream    The stream where to dump.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Load the library configuration from a binary dump and initialise accordingly.
 *
 * **Warnings** : The binary dump is raw formated, hence *a priori* platform
 * dependent. Trying to (re-)initialise an already initialised library will
 * generate an error. `pumas_finalise` must be called first.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_FORMAT_ERROR            The binary dump is not compatible
 * with the current version.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library is already initialised.
 *
 *     PUMAS_RETURN_IO_ERROR                Couldn't read from the stream.
 */
PUMAS_API enum pumas_return pumas_load(FILE * stream);

/**
 * Transport a particle according to the configured `pumas_context`.
 *
 * @param context The simulation context.
 * @param state   The initial or final state.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Depending on the *context* configuration the particle is transported through
 * one or more media, as provided by the *medium* callback. At return, the
 * particle *state* is updated.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_DENSITY_ERROR           A null or negative density was
 * encountered.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initalised.
 *
 *     PUMAS_RETURN_MEDIUM_ERROR            No propagation medium.
 *
 *     PUMAS_RETURN_MISSING_LIMIT           An external limit is needed.
 *
 *     PUMAS_RETURN_MISSING_RANDOM          A *random* callback is needed.
 */
PUMAS_API enum pumas_return pumas_transport(
    struct pumas_context * context, struct pumas_state * state);

/**
 * Print a summary of the current library configuration.
 *
 * @param stream        A stream where the summary will be formated to.
 * @param tabulation    The tabulation separator or `NULL`.
 * @param newline       The newline separator or `NULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The summary provides information on loaded materials as well as some basic
 * statistics. The *tabulation* and *newline* parameters allow to control the
 * output rendering.
 *
 * __Warnings__
 *
 * This function is **not** thread safe. A lock must be set to ensure proper
 * printout in multithreaded applications, if writing concurrently to a same
 * *stream*.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initalised.
 *
 *     PUMAS_RETURN_IO_ERROR                Couldn't write to *stream*.
 */
PUMAS_API enum pumas_return pumas_print(
    FILE * stream, const char * tabulation, const char * newline);

/**
 * Get the version tag of the library.
 *
 * @return The library tag encoded on an `int`.
 *
 * The library tag is encoded on an `int` as tag = 10^(3)*V+S, with V the
 * version index and S the subversion index.
 */
PUMAS_API int pumas_tag();

/**
 * Get info on the transported particle.
 *
 * @param particle      The type of the transported particle or `NULL`.
 * @param lifetime      The type of the transported particle, in m, or `NULL`.
 * @param mass          The mass of the transported particle, in GeV, or `NULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Retrieve info on the transported particle. If not needed, an argument can
 * be set to `NULL`.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initalised.
 */
PUMAS_API enum pumas_return pumas_particle(
    enum pumas_particle * particle, double * lifetime, double * mass);

/**
 * Return a string describing a `pumas_return` code.
 *
 * @param rc    The return code.
 * @return A static string.
 *
 * This function is analog to the standard C `strerror` function but specific
 * to PUMAS return codes. It is thread safe.
 */
PUMAS_API const char * pumas_error_string(enum pumas_return rc);

/**
 * Return a string describing a PUMAS library function.
 *
 * @param function    The library function.
 * @return a static string.
 *
 * This function is meant for verbosing when handling errors. It is thread
 * safe.
 */
PUMAS_API const char * pumas_error_function(pumas_function_t * function);

/**
 * Set or clear the error handler.
 *
 * @param handler    The error handler to set or `NULL`.
 *
 * Set the error handler callback for PUMAS library functions. If *handler* is
 * set to `NULL` error callbacks are disabled.
 *
 * __Warnings__
 *
 * This function is **not** thread safe.
 */
PUMAS_API void pumas_error_handler_set(pumas_handler_cb * handler);

/**
 * Get the current error handler.
 *
 * @return The current error handler or `NULL` if none.
 */
PUMAS_API pumas_handler_cb * pumas_error_handler_get(void);

/**
 * Catch the next error.
 *
 * @param enable   A flag for enabling or disabling error catch.
 *
 * Enable or disable the catch of the next PUMAS library error. While enabled
 * library errors will **not** trigger the error handler. Note however that only
 * the first occuring error will be caught. Call `pumas_error_raise` to enable
 * the error handler again and raise any caught error.
 *
 * __Warnings__
 *
 * This function is not thread safe. Only a single error stream can be handled
 * at a time.
 */
PUMAS_API void pumas_error_catch(int enable);

/**
 * Raise any caught error.
 *
 * @return If no error was caught `PUMAS_RETURN_SUCCESS` is returned otherwise
 * an error code is returned as detailed below.
 *
 * Raise any caught error. Error catching must have been enabled first with
 * `pumas_error_catch` otherwise a specfic `PUMAS_RETURN_RAISE_ERROR` is
 * returned. Note that calling this function disables further error's catching.
 *
 * __Warnings__
 *
 * This function is not thread safe. Only a single error stream can be handled
 * at a time.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_RAISE_ERROR    Error catching hasn't been enabled.
 *
 *     PUMAS_RETURN_*              Any caught error's code.
 */
PUMAS_API enum pumas_return pumas_error_raise(void);

/**
 * Print a formated summary of error data.
 *
 * @param stream      The output stream where to print.
 * @param rc          The error return code.
 * @param function    The faulty function or `NULL`.
 * @param error       The additional error data or `NULL`.
 *
 * The *function* and *error* parameters are optional and can be set to `NULL`.
 * The output summary is formated in JSON.
 *
 * __Warnings__
 *
 * This function is **not** thread safe. A lock must be set to ensure proper
 * printout in multithreaded applications, if writing concurrently to a same
 * *stream*.
 */
PUMAS_API void pumas_error_print(FILE * stream, enum pumas_return rc,
    pumas_function_t * function, struct pumas_error * error);

/**
 * Create a simulation context.
 *
 * @param extra_memory    The size of the user extra memory, if any is claimed.
 * @param context         A handle for the simulation context.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Create a new simulation context with a default configuration. Call
 * `pumas_context_destroy` in order to release all the memory allocated for
 * the context.
 *
 * If `extra_memory` is strictly positive the context will be extended by
 * `extra_memory` bytes for user usage. This memory can then be accessed with
 * the `user_data` field of the returned `pumas_context` structure.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 */
PUMAS_API enum pumas_return pumas_context_create(
    int extra_memory, struct pumas_context ** context);

/**
 * Destroy a simulation context.
 *
 * @param context The simulation context.
 *
 * Call on a previously created context with `pumas_context_create` in order to
 * release the corresponding dynamicaly allocated memory. On return `context`
 * is set to `NULL`.
 */
PUMAS_API void pumas_context_destroy(struct pumas_context ** context);

/**
 * Compute the total grammage that a particle can travel assuming continuous
 * energy loss.
 *
 * @param scheme      The energy loss scheme.
 * @param material    The material index.
 * @param kinetic     The initial kinetic energy, in GeV.
 * @param grammage    The grammage in kg/m^(2).
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The energy loss scheme must be one of `PUMAS_SCHEME_CSDA` or
 * `PUMAS_SCHEME_HYBRID`. For a uniform medium, divide the return value by the
 * density in order to get the corresponding total travelled distance.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The scheme of material index is
 * not valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_property_grammage(
    enum pumas_scheme scheme, int material, double kinetic, double * grammage);

/**
 * Compute the normalised total proper time spent assuming continuous
 * energy loss.
 *
 * @param scheme      The energy loss scheme.
 * @param material    The material index.
 * @param kinetic     The initial kinetic energy, in GeV.
 * @param time        The normalised proper time in kg/m^(2).
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The energy loss scheme must be one of `PUMAS_SCHEME_CSDA` or
 * `PUMAS_SCHEME_HYBRID`. Divide the returned value by the medium density
 * times *c* in order to get the proper time in unit of time.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The scheme of material index is
 * not valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_property_proper_time(
    enum pumas_scheme scheme, int material, double kinetic, double * time);

/**
 * Compute the normalised rotation angle due to a uniform magnetic field for
 * a CSDA particle.
 *
 * @param material    The material index.
 * @param kinetic     The initial kinetic energy, in GeV.
 * @param angle       The normalised rotation angle in kg/m^(2)/T.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Multiply the returned value by the transverse magnetic field amplitude and
 * divide by the density in order to get the rotation angle in radian.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The material index is not valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_property_magnetic_rotation(
    int material, double kinetic, double * angle);

/**
 * Compute the minimum kinetic energy required for travelling over a given
 * `grammage`, assuming continuous energy loss.
 *
 * @param scheme      The energy loss scheme
 * @param material    The material index.
 * @param grammage    The requested grammage, in kg/m^(2).
 * @param kinetic     The required kinetic energy in GeV.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 *  The energy loss scheme must be one of `PUMAS_SCHEME_CSDA` or
 * `PUMAS_SCHEME_HYBRID`.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The scheme of material index is
 * not valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_property_kinetic_energy(
    enum pumas_scheme scheme, int material, double grammage, double * kinetic);

/**
 * Compute the average energy loss per unit weight of material.
 *
 * @param scheme      The energy loss scheme
 * @param material    The material index.
 * @param kinetic     The kinetic energy, in GeV.
 * @param dedx        The computed energy loss in GeV/(kg/m^(2)).
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The energy loss scheme must be one of `PUMAS_SCHEME_CSDA` or
 * `PUMAS_SCHEME_HYBRID`.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The scheme of material index is
 * not valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_property_energy_loss(
    enum pumas_scheme scheme, int material, double kinetic, double * deddx);

/**
 * Compute the Multiple SCattering (MSC) 1^(st) transport path length for a
 * unit weight.
 *
 * @param material    The material index.
 * @param kinetic     The kinetic energy, in GeV.
 * @param length      The computed MSC length in kg/m^(2).
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The MSC 1^(st) transport path length, &lambda;, is related to the standard
 * deviation of the polar scattering angle's as &theta;^(2) = X/(2&lambda;),
 * with X the column depth.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The material index is not valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 *
 *     PUMAS_RETURN_VALUE_ERROR             The MSC path length is infinite.
 */
PUMAS_API enum pumas_return pumas_property_scattering_length(
    int material, double kinetic, double * length);

/**
 * The macroscopic inelastic total cross-section.
 *
 * @param material         The material index.
 * @param kinetic          The kinetic energy, in GeV.
 * @param cross_section    The computed cross-section value.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The returned cross-section value is in unit m^(2)/kg. Multiply by the
 * density in order to get the inverse of the interaction length in unit of
 * distance.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The material index is not valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_property_cross_section(
    int material, double kinetic, double * cross_section);

/**
 * The name of a material.
 *
 * @param index The material index.
 * @param material The corresponding material name.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The material name is defined in the Material Description File (MDF).
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index isn't valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR      The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_material_name(
    int index, const char ** material);

/**
 * The index of a material.
 *
 * @param material    The material name.
 * @param index       The corresponding index.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The material index corresponds to the order of declaration specified in the
 * Material Description File (MDF).
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR      The library isn't initialised.
 *
 *     PUMAS_RETURN_UNKNOWN_MATERIAL          The material isn't defined.
 */
PUMAS_API enum pumas_return pumas_material_index(
    const char * material, int * index);

/**
 * The total number of materials.
 *
 * @return The total number of known materials, basic plus composite.
 */
PUMAS_API int pumas_material_length(void);

/**
 * The number of composite materials.
 *
 * @return The number of composite materials.
 */
PUMAS_API int pumas_composite_length(void);

/**
 * Update the properties of a composite material.
 *
 * @param material  The composite material index.
 * @param fractions The vector of mass fractions of the base materials
 *                  components.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Update the composition and/or the density of a composite material.
 * `fractions` or `densities` can be `NULL` in which case the corresponding
 * property is not updated.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_DENSITY_ERROR             Some density value is null
 * or less.
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index isn't valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR      The library isn't initialised.
 *
 *     PUMAS_RETURN_MEMORY_ERROR              Couldn't allocate memory.
 */
PUMAS_API enum pumas_return pumas_composite_update(
    int material, const double * fractions, const double * densities);

/**
 * Get the properties of a composite material.
 *
 * @param index      The composite material index.
 * @param density    The composite material reference density.
 * @param components The number of base material components of the composite.
 * @param fractions  The vector of mass fractions of the base materials
 *                   components.
 * @param densities  The vector of densities of the base materials components.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Get the properties of a composite material. `density`, `components`,
 * `fractions` or `densities` can be `NULL` in which case the corresponding
 * property is not retrieved.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index isn't valid.
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR      The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_composite_properties(int index,
    double * density, int * components, double * fractions, double * densities);

/**
 * Accessor to the tabulated shared data.
 *
 * @param property    The column index of a property of interest.
 * @param scheme      The energy loss scheme.
 * @param material    The material index.
 * @param row         The kinetic value row index in the table.
 * @param value       The corresponding table value.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * For a given `material` and energy loss `scheme`, this function returns the
 * tabulated data corresponding to the given `property` column and `row` index.
 * Each row of the table corresponds to a different kinetic energy value.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             Some input index is not valid
 * (property, material or scheme).
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 */
PUMAS_API enum pumas_return pumas_table_value(enum pumas_property property,
    enum pumas_scheme scheme, int material, int row, double * value);

/**
 * The depth, i.e. number of kinetic values, of the tabulated data.
 *
 * @return The number of rows in data tables.
 */
PUMAS_API int pumas_table_length(void);

/**
 * Compute the table row index for a given property and its value.
 *
 * @param property    The column index of the property.
 * @param scheme      The energy loss scheme.
 * @param material    The material index.
 * @param value       The property value.
 * @param index       The row index from below for the given value.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * In the case of an out of bounds value the closest index value is provided
 * and `PUMAS_RETURN_VALUE_ERROR` is returned.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             Some input index is not valid
 * (property, material or scheme).
 *
 *     PUMAS_RETURN_INITIALISATION_ERROR    The library isn't initialised.
 *
 *     PUMAS_RETURN_VALUE_ERROR             The provided value is out of the
 * table.
 */
PUMAS_API enum pumas_return pumas_table_index(enum pumas_property property,
    enum pumas_scheme scheme, int material, double value, int * index);

/**
 * Create a new particle recorder.
 *
 * @param context The simulation context or `NULL`.
 * @param recorder A handle for the recorder.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The recorder is configured for a fresh start with default settings. If a
 * non `NULL` context is provided, the recorder will be linked to it.
 *
 * **Note** : though the recorder's creation requests a simulation context, the
 * returned proxy might be moved to any other context after creation. However
 * it should not be linked simultaneously to multiple concurent simulation
 * streams.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_MEMORY_ERROR    Couldn't allocate memory.
 */
PUMAS_API enum pumas_return pumas_recorder_create(
    struct pumas_context * context, struct pumas_recorder ** recorder);

/**
 * Clear all recorded frames.
 *
 * @param recorder The recorder handle.
 *
 * Erase all recorded states from the recorder and reset the frame count.
 */
PUMAS_API void pumas_recorder_clear(struct pumas_recorder * recorder);

/**
 * Destroy a particle recorder releasing all associated memory.
 *
 * @param recorder The recorder handle.
 *
 * **Note** : The recorder is cleared before beeing destroyed. At return
 * `recorder` is set to `NULL`.
 */
PUMAS_API void pumas_recorder_destroy(struct pumas_recorder ** recorder);

/**
 * User supplied callback for memory allocation.
 *
 * @param size    The number of memory bytes to allocate.
 * @return The address of the allocated memory or `NULL` in case of faillure.
 *
 * The provided callback must conform to the `malloc` semantic and behaviour.
 */
typedef void * pumas_allocate_cb(size_t size);

/**
 * Set the memory allocation function for the PUMAS library.
 *
 * @param allocator    The user supplied memory allocator, or `NULL`.
 *
 * This function allows to specify a custom memory allocation function for
 * PUMAS. Passing a `NULL` value results in PUMAS using its default allocator,
 * i.e. `malloc`.
 *
 * __Warnings__
 *
 * This function is **not** thread safe.
 */
PUMAS_API void pumas_memory_allocator(pumas_allocate_cb * allocator);

/**
 * User supplied callback for memory re-allocation.
 *
 * @param ptr     The address of the memory to reallocate.
 * @param size    The number of memory bytes requested for the reallocation.
 * @return The address of the re-allocated memory or `NULL` in case of faillure.
 *
 * The provided callback must conform to the `realloc` semantic and behaviour.
 */
typedef void * pumas_reallocate_cb(void * ptr, size_t size);

/**
 * Set the memory re-allocation function for the PUMAS library.
 *
 * @param reallocator    The user supplied memory reallocator, or `NULL`.
 *
 * This function allows to specify a custom memory re-allocation function for
 * PUMAS. Passing a `NULL` value results in PUMAS using its default
 * reallocator, i.e. `realloc`.
 *
 * __Warnings__
 *
 * This function is **not** thread safe.
 */
PUMAS_API void pumas_memory_reallocator(pumas_reallocate_cb * reallocator);

/**
 * User supplied callback for memory deallocation.
 *
 * @param size    The address of the memory to deallocate.
 *
 * The provided callback must conform to the `free` semantic and behaviour.
 */
typedef void pumas_deallocate_cb(void * ptr);

/**
 * Set the memory deallocation function for the PUMAS library.
 *
 * @param deallocator    The user supplied memory deallocator, or `NULL`.
 *
 * This function allows to specify a custom memory deallocation function for
 * PUMAS. Passing a `NULL` value results in PUMAS using its default
 * deallocator, i.e. `free`.
 *
 * __Warnings__
 *
 * This function is **not** thread safe.
 */
PUMAS_API void pumas_memory_deallocator(pumas_deallocate_cb * deallocator);

#ifdef __cplusplus
}
#endif
#endif
