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
 * Particle types managed by the PUMAS transport engine.
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
        /** The total magnetic rotation angle, in rad kg/m^(3). */
        PUMAS_PROPERTY_MAGNETIC_ROTATION,
        /** The particle proper time, in kg/m^(2). */
        PUMAS_PROPERTY_PROPER_TIME,
        /** The macroscopic elastic scattering 1^(st) path length, in kg/m^(2).
         */
        PUMAS_PROPERTY_SCATTERING_LENGTH
};

/**
 * Modes for the Monte Carlo transport.
 */
enum pumas_mode {
        /** All energy losses are disabled.
         *
         * **Note** : This mode is provided for test purpose only. Running
         * without energy losses requires specifying a range or grammage limit.
         */
        PUMAS_MODE_VIRTUAL = -1,
        /** Energy losses are purely determinstic as given by the Continuously
         * Slowing Down Approximation (CSDA).
         */
        PUMAS_MODE_CSDA = 0,
        /** Energy losses are simulated using an hybrid Monte-Carlo scheme with
         * hard stochastic interactions and a deterministic Continuous Energy
         * Loss (CEL).
         */
        PUMAS_MODE_HYBRID = 1,
        /** In addition to the hybrid scheme small fluctuations of the CEL are
         * also simulated.
         */
        PUMAS_MODE_DETAILED = 2,
        /** Decays are disabled, i.e. muons or taus are stable. */
        PUMAS_MODE_STABLE = 0,
        /**
         * Decays are accounted for by a weight factor. This is efficient
         * for muons but irrelevant -numerically instable- for the forward
         * transport of taus since they decay in flight. Hence it is
         * disabled in the latter case.
         */
        PUMAS_MODE_WEIGHT = 1,
        /** Decays are accounted for with a specific Monte-Carlo process.
         *
         * **Note** : the transported particle stops at the dcay vertex but
         * its decay is not simulated, i.e. no daughter particles are
         * generated. */
        PUMAS_MODE_DECAY = 2,
        /** Do a classical forward Monte Carlo transport. */
        PUMAS_MODE_FORWARD = 0,
        /** Do a reverse Monte Carlo transport. */
        PUMAS_MODE_BACKWARD = 1,
        /** Fully simulate the (multiple)scattering. */
        PUMAS_MODE_FULL_SPACE = 0,
        /** Neglect the tranverse scattering, i.e. only simulate the energy
         * loss.
         *
         * **Note** : charged particles are still deflected by external
         * electromagnetic fields.
         */
        PUMAS_MODE_LONGITUDINAL = 1
};

/** Return codes for the API functions. */
enum pumas_return {
        /** Execution was successful. */
        PUMAS_RETURN_SUCCESS = 0,
        /** The requested accuracy is not valid. */
        PUMAS_RETURN_ACCURACY_ERROR,
        /** The requested cutoff is not valid. */
        PUMAS_RETURN_CUTOFF_ERROR,
        /** End of file was reached. */
        PUMAS_RETURN_END_OF_FILE,
        /** The specified decay mode is not valid. */
        PUMAS_RETURN_DECAY_ERROR,
        /** Some medium has a wrong density value. */
        PUMAS_RETURN_DENSITY_ERROR,
        /** Some medium has a wrong density value. */
        PUMAS_RETURN_DIRECTION_ERROR,
        /** A non unit direction is provided. */
        PUMAS_RETURN_INCOMPLETE_FILE,
        /** Some index is out of validity range. */
        PUMAS_RETURN_INDEX_ERROR,
        /** The Physics is not initialised or a NULL pointer was provided. */
        PUMAS_RETURN_PHYSICS_ERROR,
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
        /** An invalid (unknown) DCS model was requested. */
        PUMAS_RETURN_MODEL_ERROR,
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

/** Flags for transport events. */
enum pumas_event {
        /** No event occured or is foreseen. */
        PUMAS_EVENT_NONE = 0,
        /** A kinetic energy limit was reached or is foreseen. */
        PUMAS_EVENT_LIMIT_ENERGY = 1,
        /** A distance limit was reached or is foreseen. */
        PUMAS_EVENT_LIMIT_DISTANCE = 2,
        /** A grammage limit was reached or is foreseen. */
        PUMAS_EVENT_LIMIT_GRAMMAGE = 4,
        /** A proper time limit was reached or is foreseen. */
        PUMAS_EVENT_LIMIT_TIME = 8,
        /** Shortcut for any external limit. */
        PUMAS_EVENT_LIMIT = 15,
        /** A change of medium occured or is foreseen. */
        PUMAS_EVENT_MEDIUM = 16,
        /** A Bremsstrahlung occured or is foreseen. */
        PUMAS_EVENT_VERTEX_BREMSSTRAHLUNG = 32,
        /** A Pair creation occured or is foreseen. */
        PUMAS_EVENT_VERTEX_PAIR_CREATION = 64,
        /** A Photonuclear interaction occured or is foreseen. */
        PUMAS_EVENT_VERTEX_PHOTONUCLEAR = 128,
        /** A Delta ray occured or is foreseen. */
        PUMAS_EVENT_VERTEX_DELTA_RAY = 256,
        /** Shortcut for any Discrete Energy Loss (DEL). */
        PUMAS_EVENT_VERTEX_DEL = 480,
        /** A hard Coulombian interaction occured or is foreseen. */
        PUMAS_EVENT_VERTEX_COULOMB = 512,
        /** A decay has occured or is foreseen. */
        PUMAS_EVENT_VERTEX_DECAY = 1024,
        /** Shortcut for any interaction vertex. */
        PUMAS_EVENT_VERTEX = 2016,
        /** The particle has a nul or negative weight. */
        PUMAS_EVENT_WEIGHT = 2048,
        /** Extra flag for records tagging the 1st transport step. */
        PUMAS_EVENT_START = 4096,
        /** Extra flag for records tagging the last transport step. */
        PUMAS_EVENT_STOP = 8192
};

/** Indices for customizable Physics processes. */
enum pumas_process {
        /** The Bremstrahlung process */
        PUMAS_PROCESS_BREMSSTRAHLUNG = 0,
        /** The e+e- pair production process */
        PUMAS_PROCESS_PAIR_PRODUCTION,
        /** The photonuclear process */
        PUMAS_PROCESS_PHOTONUCLEAR
};

/** Indices for PUMAS library constants. */
enum pumas_constant {
        /** The Avogadro number in mol. */
        PUMAS_CONSTANT_AVOGADRO_NUMBER = 0,
        /** The electron mass in GeV/c^2. */
        PUMAS_CONSTANT_ELECTRON_MASS,
        /** The muon decay length in m. */
        PUMAS_CONSTANT_MUON_C_TAU,
        /** The muon mass in GeV/c^2. */
        PUMAS_CONSTANT_MUON_MASS,
        /** The tau decay length in m. */
        PUMAS_CONSTANT_TAU_C_TAU,
        /** The tau mass in GeV/c^2. */
        PUMAS_CONSTANT_TAU_MASS,
        /** The number of PUMAS constants.  */
        PUMAS_N_CONSTANTS
};

/**
 * Container for a Monte-Carlo state.
 */
struct pumas_state {
        /** The particle's electric charge. Note that non physical values,
         * i.e. different from 1 or -1, could be set. */
        double charge;
        /** The current kinetic energy, in GeV. */
        double energy;
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
        /** The momentum's unit direction. Must be normalised to one. */
        double direction[3];
        /** Status flag telling if the particle has decayed or not.  */
        int decayed;
};

/**
 * The local properties of a propagation medium.
 */
struct pumas_locals {
        /** The material local density, in kg/m^3. Setting a null or negative
         * value results in the material's default density being used.
         */
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
 * @return The size of local inhomogeneities (see below).
 *
 * The callback must return a length, in m, consistent with the size of the
 * propagation medium inhomogeneities, e. g. &rho; / |&nabla; &rho;| for a
 * density gradient. Returning zero or less signs that the propagation medium is
 * fully uniform.
 *
 * **Note** that inhomogeneities modelled by the `pumas_locals` callback must be
 * **continuous**. If the geometry has a density or magnetic field discontinuity
 * then this must be modelled by using separate media on both sides of the
 * discontinuity.
 *
 * **Warning** : it is an error to return zero or less for any position of the
 * medium if at least one area is not uniform. Instead one should use two
 * different media even though they have the same material base.
 *
 */
typedef double pumas_locals_cb (struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals);

/**
 * Description of a propagation medium.
 *
 * A propagation medium is fully defined by:
 *
 * - a `material` composition with a uniform relative content.
 * - `pumas_locals` properties set by a user provided `pumas_locals_cb`
 * callback or `NULL` e.g. if using the material's default density.
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
         * Setting a `NULL` callback results in the material's default density
         * being used with no magnetic field.
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
        /** The corresponding step event. */
        enum pumas_event event;
        /** Link to the next frame in the record. */
        struct pumas_frame * next;
};

struct pumas_context;
/** A user supplied recorder callback.
 * @param context The recording simulation context.
 * @param state   The recorded particle state.
 * @param medium  The corresponding medium.
 * @param event   The step event.
 *
 * This callback allows to customize the recording of PUMAS Monte-Carlo events.
 *
 * **Note** : by default the recorder uses an in-memory copy with dynamic
 * allocation. Setting a custom recorder disables the default recording.
 */
typedef void pumas_recorder_cb (struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium,
    enum pumas_event event);

/**
 * A handle for recording Monte-Carlo frames.
 *
 * This structure is a proxy for recording Monte-Carlo states and/or accessing
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
         * The sampling period of the recorder, If set to zero or less only
         * medium changes are recorded. Defaults to 1, i.e. all Monte-Carlo
         * steps are recorded.
         */
        int period;
        /**
         * Link to an external (user supplied) recording callback. Note that
         * setting this value disables the in-memory frame recording. Defaults
         * to `NULL`.
         */
        pumas_recorder_cb * record;
        /**
         * A pointer to additional memory, if any is requested at
         * initialisation.
         */
        void * user_data;
};

/** Return codes for the `pumas_medium_cb` callback. */
enum pumas_step {
        /** The proposed step is cross-checked bu PUMAS */
        PUMAS_STEP_CHECK = 0,
        /** The raw proposed step is used by PUMAS */
        PUMAS_STEP_RAW
};

/**
 * Callback for locating the propagation medium of a `pumas_state`.
 *
 * @param context   The Monte-Carlo context requiring a medium.
 * @param state     The Monte-Carlo state for which the medium is requested.
 * @param medium    A pointer to store the medium or `NULL` if not requested.
 * @param step      The proposed step size or zero or less for an infinite
 *                    medium. If not requested this points to `NULL`.
 * @return If the proposed step size should be cross-checked by PUMAS
 * `PUMAS_STEP_CHECK` should be returned otherwise `PUMAS_STEP_RAW`.
 *
 * If *step* is not `NULL`, this callback must propose a Monte-Carlo stepping
 * distance, in m, consistent with the geometry. Note that returning zero or
 * less signs that the corresponding medium has no boundaries. When *medium* is
 * not `NULL` it must be set to the located `pumas_medium`.
 *
 * In addition the user must return a `pumas_step` enum indicating if the
 * proposed *step* needs to be cross-checked by PUMAS or if it should be used
 * raw. Managing steps that end on a geometry boundary can be tricky
 * numerically. Therefore it is recommended to return `PUMAS_STEP_CHECK` if you
 * are unsure of what to do since it is more robust. The raw mode is usefull if
 * your geometry engine already performs those checks in order to avoid double
 * work.
 *
 * **Warning** : it is an error to return zero or less for any state if the
 * extension is finite.
 *
 * **Warning** : in backward Monte Carlo mode the particle is propagated reverse
 * to the state direction. The user must take care to provide a *step* size
 * accordingly, i.e. consistent with the geometry in both forward and backward
 * modes.
 */
typedef enum pumas_step pumas_medium_cb (
    struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium ** medium, double * step);

/**
 * Generic function pointer.
 *
 * This is a generic function pointer used to identify the library functions,
 * e.g. for error handling.
 */
typedef void pumas_function_t (void);

/**
 * Callback for error handling.
 *
 * @param rc          The PUMAS return code.
 * @param caller      The API function where the error occured.
 * @param message     Brief description of the error.
 *
 * The user can override the PUMAS default error handler by providing its own
 * error handler. It will be called at the return of any PUMAS library function
 * providing an error code.
 */
typedef void pumas_handler_cb (enum pumas_return rc, pumas_function_t * caller,
    const char * message);

/**
 * Callback providing a stream of pseudo random numbers.
 *
 * @param context The simulation context requiring a random number.
 * @return A uniform pseudo random number in [0;1].
 *
 * **Note** : this is the only random stream used by PUMAS. If overriding the
 * default `pumas_context` callback then the user must unsure proper behaviour,
 * i.e. that a flat distribution in [0;1] is indeed returned.
 *
 * **Warning** : if multiple contexts are used the user must also ensure that
 * this callback is thread safe, e.g. by using independant streams for each
 * context or a locking mechanism in order to share a single random stream.
 * The default `pumas_context` random callback uses distinct random streams per
 * context which ensures thread safety.
 */
typedef double pumas_random_cb (struct pumas_context * context);

/** Mode flags for the Monte Carlo transport. */
struct pumas_context_mode {
        /**
        * The scheme used for the computation of energy losses. Default
        * is `PUMAS_MOED_DETAILED`.
        */
        enum pumas_mode energy_loss;
        /**
        * The mode for handling decays. Default is `PUMAS_MODE_WEIGHT`
        * for a muon or `PUMAS_MODE_DECAY` for a tau. Set this to
        * `PUMAS_MODE_STABLE` in order to disable decays at all.
        */
        enum pumas_mode decay;
        /**
        * Direction of the Monte Carlo flow. Default is
        * `PUMAS_MODE_FORWARD`. Set this to `PUMAS_MODE_BACKWARD` for a
        * reverse Monte Carlo.
        */
        enum pumas_mode direction;
        /**
        * Algorithm for the simulation of the scattering. Default is
        * `PUMAS_MODE_FULL_SPACE`. Other option is
        * `PUMAS_MODE_LONGITUDNAL` which neglects any transverse
        * scattering.
        */
        enum pumas_mode scattering;
};

/** External limits for the Monte Carlo transport. */
struct pumas_context_limit {
        /**
         * The minimum kinetic energy for forward transport, or the
         * maximum one for backward transport, in GeV.
         */
        double energy;
        /** The maximum travelled distance, in m. */
        double distance;
        /** The maximum travelled grammage, in kg/m^2. */
        double grammage;
        /** The maximum travelled proper time, in m/c. */
        double time;
};

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
 * + The newly created context is configured with a default pseudo random engine
 * using the Mersenne Twister algorithm. The engine can be explictly seeded with
 * the `pumas_context_random_initialise` function. Note that two contexts seeded
 * with the same value are 100% correlated by construction. If no seed is
 * provided one is picked randomly from the OS, e.g. from `/dev/urandom` on
 * UNIX. A custom random engine can be used instead of the default one by
 * overriding the *random* callback.
 *
 * + Note that for `energy`, `distance`, `grammage` or `time` external
 * limits to be taken into account, the corresponding events must be activated
 * as well, with the `event` flag.
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

        /** Monte Carlo transport mode. */
        struct pumas_context_mode mode;
        /**
         * The events that might stop the transport. Default is
         * `PUMAS_EVENT_NONE`, i.e. the transport stops only if the particle
         * exits the simulation media, or if it looses all of its energy.
         */
        enum pumas_event event;

        /** External limits for the Monte Carlo transport. */
        struct pumas_context_limit limit;

        /** Tuning parameter for the accuracy of the Monte Carlo integration.
         *
         * The Monte Carlo transport is discretized in elementary steps. This
         * parameter directly controls the length of these steps. The smaller
         * the *accuracy* value the smaller the step length hence the longer
         * the Monte Carlo simulation.
         */
        double accuracy;
};

/**
 * Opaque handle for Physics tables
 */
struct pumas_physics;

/**
 * Sternheimer coefficients for the parameterisation of the density effect.
 */
struct pumas_physics_density_effect {
        /** Sternheimer *a* Coefficient. */
        double a;
        /** Sternheimer *k* Coefficient. */
        double k;
        /** Sternheimer *x0* Coefficient. */
        double x0;
        /** Sternheimer *x1* Coefficient. */
        double x1;
        /** Sternheimer *Cbar* Coefficient. */
        double Cbar;
        /** Sternheimer *delta0* Coefficient. */
        double delta0;
};

/**
 * Prototype for a Differential Cross-Section (DCS).
 *
 * @param Z       The charge number of the target atom.
 * @param A       The mass number of the target atom.
 * @param m       The projectile rest mass, in GeV
 * @param K       The projectile kinetic energy, in GeV.
 * @param q       The projectile energy loss, in GeV.
 * @return The corresponding value of the atomic DCS, in m^(2)/GeV.
 *
 * The `pumas_dcs_get` function allows to retrieve the DCS for a given process
 * and model. Extra DCSs can be registered with the `pumas_dcs_register`
 * function.
 *
 * **Note** : only the Bremsstrahlung, pair creation and photonuclear processes
 * can be redefined.
 */
typedef double pumas_dcs_t (double Z, double A, double m, double K, double q);

/**
 */
struct pumas_physics_settings {
        /** Relative cutoff between continous and discrete energy.
         *
         * Setting a null or negative value results in the default cutoff value
         * to be used i.e. 5% which is a good compromise between speed and
         * accuracy for transporting a continuous spectrum (see e.g.  Sokalski
         * et al., https://doi.org/10.1103/PhysRevD.64.074015).
         */
        double cutoff;
        /** Physics model for the Bremsstrahlung process.
         *
         *  Available models are:
         *
         *  - "KKP": Kelner, Kokoulin & Petrukhin.
         *  - "ABB": Andreev, Bezrukov & Bugaev.
         *
         * Setting a `NULL` value results in PUMAS default Bremsstrahlung model
         * to be used, i.e. "KKP".
         * */
        const char * bremsstrahlung;
        /** Physics model for e^(+)e^(-) pair production.
         *
         *  Available models are:
         *
         *  - "KKP": Kelner, Kokoulin & Petrukhin.
         *
         * Setting a `NULL` value results in PUMAS default pair production model
         * to be used, i.e. "KKP".
         */
        const char * pair_production;
        /** Physics model for photonuclear interactions. 
         *
         *  Available models are:
         *
         *  - "DRSS": Dutta, Reno, Sarcevic & Seckel.
         *
         * Setting a `NULL` value results in PUMAS default photonuclear model to
         * be used, i.e. "DRSS".
         */
        const char * photonuclear;
};

/**
 * Initialise the Physics.
 *
 * @param physics      Handle for the Physics tables.
 * @param particle     The type of the particle to transport.
 * @param mdf_path     The path to a Material Description File (MDF), or `NULL`.
 * @param dedx_path    The path to the energy loss tabulation(s), or `NULL`.
 * @param settings     Extra physics settings or `NULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Initialise the Physics from a MDF and a set of energy loss tabulations. Load
 * the materials data and precompute various properties. *mdf_path* and/or
 * *dedx_path* can be set to `NULL`. If so the corresponding path is read from
 * the `PUMAS_MDF` or `PUMAS_DEDX` environment variable.
 *
 * Optionaly extra physics settings can be specified by providing a
 * `pumas_physics_settings` structure. If `NULL` is provided then PUMAS default
 * physics settings are used which should perform well for most use cases.
 *
 * Call `pumas_physics_destroy` in order to unload the Physics and release the
 * corresponding alocated memory.
 *
 * **Warnings** : this function is not thread safe.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_CUTOFF_ERROR            A bad cutoff value was provided.
 *
 *     PUMAS_RETURN_END_OF_FILE             And unexpected EOF occured.
 *
 *     PUMAS_RETURN_FORMAT_ERROR            A file has a wrong format.
 *
 *     PUMAS_RETURN_INCOMPLETE_FILE         There are missing entries in
 * the MDF.
 *
 *     PUMAS_RETURN_IO_ERROR                A file couldn't be read.
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 *
 *     PUMAS_RETURN_MODEL_ERROR             A requested DCS model is not valid.
 *
 *     PUMAS_RETURN_PATH_ERROR              A file couldn't be opened.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           A `NULL` physics pointer was
 * provided.
 *
 *     PUMAS_RETURN_TOO_LONG                Some XML node in the MDF is
 * too long.
 *
 *     PUMAS_RETURN_UNDEFINED_DEDX          No energy loss path was provided.
 *
 *     PUMAS_RETURN_UNDEFINED_MDF           No MDF was provided.
 *
 *     PUMAS_RETURN_UNKNOWN_ELEMENT         An element in the MDF wasn't
 * defined.
 *
 *     PUMAS_RETURN_UNKNOWN_MATERIAL        An material in the MDF wasn't
 * defined.
 *
 *     PUMAS_RETURN_UNKNOWN_PARTICLE        The given type is not supported.
 */
PUMAS_API enum pumas_return pumas_physics_create(
    struct pumas_physics ** physics, enum pumas_particle particle,
    const char * mdf_path, const char * dedx_path,
    const struct pumas_physics_settings * settings);

/**
 * Destroy a Physics instance.
 *
 * @param physics      Handle for the Physics tables.
 *
 * Finalise the Physics and free the shared memory. Call
 * `pumas_physics_create` in order to reload the Physics.
 *
 * **Warnings** : This function is not thread safe. Finalising the Physics
 * doesn't release the memory allocated for any `pumas_context`.
 */
PUMAS_API void pumas_physics_destroy(struct pumas_physics ** physics);

/**
 * Dump the Physics configuration to a stream.
 *
 * @param physics   Handle for the Physics tables.
 * @param stream    The stream where to dump.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Dump the Physics tables to a stream as a binary object. Note that only
 * globally shared data are dumped, i.e. material properties and tables as read
 * from a MDF. Simulation contexts, media, recorders, ect. ... are not. This
 * binary dump allows for a fast initialisation of the Physics in subsequent
 * uses.
 *
 * **Warnings** : The binary dump is raw formated, hence *a priori* platform
 * dependent.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 *
 *     PUMAS_RETURN_PATH_ERROR              The output stream in invalid (null).
 *
 *     PUMAS_RETURN_IO_ERROR                Couldn't write to the stream.
 */
PUMAS_API enum pumas_return pumas_physics_dump(
    const struct pumas_physics * physics, FILE * stream);

/**
 * Load the Physics from a binary dump.
 *
 * @param physics   Handle for the Physics tables.
 * @param stream    The stream to load from.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Load the Physics tables from a binary dump and initialise accordingly.
 *
 * **Warnings** : The binary dump is raw formated, hence *a priori* platform
 * dependent. Trying to (re-)initialise an already initialised Physics will
 * generate an error. `pumas_physics_destroy` must be called first.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_FORMAT_ERROR            The binary dump is not compatible
 * with the current version.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 *
 *     PUMAS_RETURN_PATH_ERROR              The input stream in invalid (null).
 *
 *     PUMAS_RETURN_IO_ERROR                Couldn't read from the stream.
 */
PUMAS_API enum pumas_return pumas_physics_load(
    struct pumas_physics ** physics, FILE * stream);

/**
 * Get the cutoff value used by the physics.
 *
 * @param context Handle for the physics tables.
 * @return The cutoff value or -1 if the physics is not properly initialised.
 *
 * The cutoff value between continuous and discrete energy losses is specified
 * during the physics initialisation with `pumas_physics_create`. It cannot be
 * modified afterwards. Instead a new physics object must be created.
 */
PUMAS_API double pumas_physics_cutoff(const struct pumas_physics * physics);

/**
 * Transport a particle according to the configured `pumas_context`.
 *
 * @param context The simulation context.
 * @param state   The initial state or the final state at return.
 * @param event   The `pumas_event` at return, or `ǸULL`.
 * @param media   The initial and final media, or `ǸULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Depending on the *context* configuration the particle is transported through
 * one or more media, as provided by the *medium* callback. At return, the
 * particle *state* is updated. If *event* is not `NULL` it will be filled
 * with the end step event flag. In addition, if *media* is not `ǸULL` it will
 * contain the initial (`media[0]`) and final (`media[1]`) media crossed by the
 * particle.
 *
 * **Note**: the state direction must be a unit vector. Otherwise an error
 * is returned (see below).
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_ACCURACY_ERROR          The requested accuracy is not valid.
 *
 *     PUMAS_RETURN_DENSITY_ERROR           A null or negative density was
 * encountered.
 *
 *     PUMAS_RETURN_DIRECTION_ERROR         A non unit direction was provided.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initalised.
 *
 *     PUMAS_RETURN_MEDIUM_ERROR            No propagation medium.
 *
 *     PUMAS_RETURN_MISSING_LIMIT           An external limit is needed.
 *
 *     PUMAS_RETURN_MISSING_RANDOM          A *random* callback is needed.

 *     PUMAS_RETURN_VALUE_ERROR             State or context is `NULL`.
 */
PUMAS_API enum pumas_return pumas_context_transport(
    struct pumas_context * context, struct pumas_state * state,
    enum pumas_event * event, struct pumas_medium * media[2]);

/**
 * Print a summary of the Physics configuration.
 *
 * @param physics       Handle for the Physics tables.
 * @param stream        A stream where the summary will be formated to.
 * @param tabulation    The tabulation separator or `NULL`.
 * @param newline       The newline separator or `NULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The summary is JSON formated. It provides information on loaded materials as
 * well as some basic statistics. The *tabulation* and *newline* parameters
 * allow to control the output rendering.
 *
 * __Warnings__
 *
 * This function is **not** thread safe. A lock must be set to ensure proper
 * printout in multithreaded applications, if writing concurrently to a same
 * *stream*.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initalised.
 *
 *     PUMAS_RETURN_IO_ERROR                Couldn't write to *stream*.
 */
PUMAS_API enum pumas_return pumas_physics_print(
    const struct pumas_physics * physics, FILE * stream,
    const char * tabulation, const char * newline);

/**
 * Get the version of the PUMAS library.
 *
 * @return The library version encoded on an `int`.
 *
 * The library version is encoded on an `int` as 100*(MAJOR.MINOR). E.g.
 * PUMAS version `0.15` would yield `15`.
 */
PUMAS_API int pumas_version();

/**
 * Get info on the transported particle.
 *
 * @param physics       Handle for the Physics tables.
 * @param particle      The type of the transported particle or `NULL`.
 * @param lifetime      The proper lifetime, in m/c, or `NULL`.
 * @param mass          The mass of the transported particle, in GeV, or `NULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Retrieve info on the transported particle. If not needed, an argument can
 * be set to `NULL`.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_PHYSICS_ERROR    The Physics is not initalised.
 */
PUMAS_API enum pumas_return pumas_physics_particle(
    const struct pumas_physics * physics, enum pumas_particle * particle,
    double * lifetime, double * mass);

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
 * Set the error handler callback for PUMAS library functions.  The user can
 * override the PUMAS default error handler by providing its own error handler.
 * If *handler* is set to `NULL` error callbacks are disabled.
 *
 * __Warnings__
 *
 * This function is **not** thread safe.
 */
PUMAS_API void pumas_error_handler_set(pumas_handler_cb * handler);

/**
 * Get the current error handler. By default PUMAS is configured to printout to
 * stderr whenever an error occurs and to exit back to the OS. See
 * `pumas_error_handler_set` in order to override this behaviour.
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
 * Create a simulation context.
 *
 * @param context         A handle for the simulation context.
 * @param physics         Handle for the physics tables.
 * @param extra_memory    The size of the user extra memory, if any is claimed.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Create a new simulation context with a default configuration. Call
 * `pumas_context_destroy` in order to release all the memory allocated for
 * the context. Note that Physics tables must have been initialised / loaded
 * first.
 *
 * If `extra_memory` is strictly positive the context will be extended by
 * `extra_memory` bytes for user usage. This memory can then be accessed with
 * the `user_data` field of the returned `pumas_context` structure.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_context_create(
    struct pumas_context ** context, const struct pumas_physics * physics,
    int extra_memory);

/**
 * Set the random seed of a simulation context.
 *
 * @param context         A handle for the simulation context.
 * @param seed            The random seed or `NULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Set the seed for the default pseudo random engine. If `NULL` is provided
 * the seed is initialised from the OS, e.g. using `/dev/urandom` on UNIX.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 *
 *     PUMAS_RETURN_PATH_ERROR              The OS random stream could not be
 * read.
 */
PUMAS_API enum pumas_return pumas_context_random_seed_set(
    struct pumas_context * context, const unsigned long * seed);

/**
 * Get the random seed of a simulation context.
 *
 * @param context         A handle for the simulation context.
 * @param seed            The random seed.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Get the seed for the default pseudo random engine. Note that if no seed
 * was explicitly provided, e.g. using `pumas_context_random_seed_set`, then
 * the seed is initialised from the OS, e.g. using `/dev/urandom` on UNIX.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 *
 *     PUMAS_RETURN_PATH_ERROR              The OS random stream could not be
 * read.
 */
PUMAS_API enum pumas_return pumas_context_random_seed_get(
    struct pumas_context * context, unsigned long * seed);

/**
 * Load the random engine state from a stream.
 *
 * @param context         A handle for the simulation context.
 * @param stream          The stream to load from.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Restore the context default random engine state from a stream. See the
 * `pumas_context_random_dump` function for dumping the random engine state
 * to a file.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_FORMAT_ERROR            The binary dump is not compatible
 * with the current version.
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 *
 *     PUMAS_RETURN_IO_ERROR                Couldn't read from the stream.
 *
 *     PUMAS_RETURN_PATH_ERROR              The input stream in invalid (null).
 */
PUMAS_API enum pumas_return pumas_context_random_load(
    struct pumas_context * context, FILE * stream);

/**
 * Dump the random engine state to a stream.
 *
 * @param context         A handle for the simulation context.
 * @param stream          The stream to dump to.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Save the context default random engine state to a stream. See the
 * `pumas_context_random_load` function for loading back the random engine
 * state from a file.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 *
 *     PUMAS_RETURN_IO_ERROR                Couldn't write to the stream.
 *
 *     PUMAS_RETURN_PATH_ERROR              The input stream in invalid (null).*
 */
PUMAS_API enum pumas_return pumas_context_random_dump(
    struct pumas_context * context, FILE * stream);

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
 * Get the Physics used by a simulation context.
 *
 * @param context The simulation context.
 * @return A handle for the Physics tables or `NULL`.
 *
 * The set of Physics tables used by a `pumas_context` cannot be changed.
 * Instead a new context must be created if different Physics is needed.
 */
PUMAS_API const struct pumas_physics * pumas_context_physics_get(
    const struct pumas_context * context);

/**
 * Get the total grammage that a particle can travel assuming continuous
 * energy loss.
 *
 * @param physics     Handle for the Physics tables.
 * @param scheme      The energy loss scheme.
 * @param material    The material index.
 * @param energy      The initial kinetic energy, in GeV.
 * @param grammage    The grammage in kg/m^(2).
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The energy loss scheme must be one of `PUMAS_MODE_CSDA` or
 * `PUMAS_MODE_HYBRID`. For a uniform medium, divide the return value by the
 * density in order to get the corresponding total travelled distance.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The scheme of material index is
 * not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_property_grammage(
    const struct pumas_physics * physics, enum pumas_mode scheme,
    int material, double energy, double * grammage);

/**
 * Get the normalised total proper time spent assuming continuous energy loss.
 *
 * @param physics     Handle for the Physics tables.
 * @param scheme      The energy loss scheme.
 * @param material    The material index.
 * @param energy      The initial kinetic energy, in GeV.
 * @param time        The normalised proper time in kg/m^(2).
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The energy loss scheme must be one of `PUMAS_MODE_CSDA` or
 * `PUMAS_MODE_HYBRID`. Divide the returned value by the medium density
 * times *c* in order to get the proper time in unit of time.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The scheme of material index is
 * not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_property_proper_time(
    const struct pumas_physics * physics, enum pumas_mode scheme,
    int material, double energy, double * time);

/**
 * Get the normalised rotation angle due to a uniform magnetic field for
 * a CSDA particle.
 *
 * @param physics     Handle for the Physics tables.
 * @param material    The material index.
 * @param energy      The initial kinetic energy, in GeV.
 * @param angle       The normalised rotation angle in rad kg/m^(3)/T.
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
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_property_magnetic_rotation(
    const struct pumas_physics * physics, int material, double energy,
    double * angle);

/**
 * Get the minimum kinetic energy required for travelling over a given
 * `grammage`, assuming continuous energy loss.
 *
 * @param physics     Handle for the Physics tables.
 * @param scheme      The energy loss scheme
 * @param material    The material index.
 * @param grammage    The requested grammage, in kg/m^(2).
 * @param energy      The required kinetic energy in GeV.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 *  The energy loss scheme must be one of `PUMAS_MODE_CSDA` or
 * `PUMAS_MODE_HYBRID`.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The scheme of material index is
 * not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_property_kinetic_energy(
    const struct pumas_physics * physics, enum pumas_mode scheme,
    int material, double grammage, double * energy);

/**
 * Get the average energy loss per unit weight of material.
 *
 * @param physics     Handle for the Physics tables.
 * @param scheme      The energy loss scheme
 * @param material    The material index.
 * @param energy      The kinetic energy, in GeV.
 * @param dedx        The computed energy loss in GeV/(kg/m^(2)).
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The energy loss scheme must be one of `PUMAS_MODE_CSDA` or
 * `PUMAS_MODE_HYBRID`.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The scheme of material index is
 * not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_property_energy_loss(
    const struct pumas_physics * physics, enum pumas_mode scheme,
    int material, double energy, double * dedx);

/**
 * Get the Multiple SCattering (MSC) 1^(st) transport path length for a
 * unit weight.
 *
 * @param physics     Handle for the Physics tables.
 * @param material    The material index.
 * @param energy      The kinetic energy, in GeV.
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
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 *
 *     PUMAS_RETURN_VALUE_ERROR             The MSC path length is infinite.
 */
PUMAS_API enum pumas_return pumas_physics_property_scattering_length(
    const struct pumas_physics * physics, int material, double energy,
    double * length);

/**
 * Get the macroscopic total inelastic cross-section.
 *
 * @param physics          Handle for the Physics tables.
 * @param material         The material index.
 * @param energy           The kinetic energy, in GeV.
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
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_property_cross_section(
    const struct pumas_physics * physics, int material, double energy,
    double * cross_section);

/**
 * The total number of atomic elements.
 *
 * @param physics    Handle for the Physics tables.
 * @return The total number of known atomic elements.
 */
PUMAS_API int pumas_physics_element_length(
    const struct pumas_physics * physics);

/**
 * Get the properties of an atomic element.
 *
 * @param physics    Handle for the Physics tables.
 * @param index      The element index.
 * @param Z          The element charge number.
 * @param A          The element mass number in g/mol.
 * @param I          The element mean excitation energy in GeV.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Get the properties of an atomic element. `Z`, `A` or `I` can be `NULL`
 * in which case the corresponding property is not retrieved.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index is not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR             The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_element_properties(
    const struct pumas_physics * physics, int index, double * Z, double * A,
    double * I);

/**
 * The name of an atomic element given its index.
 *
 * @param physics    Handle for the Physics tables.
 * @param index      The atomic element index.
 * @param element    The corresponding element name.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The atomic element name is defined in the Material Description File (MDF).
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index is not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR             The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_element_name(
    const struct pumas_physics * physics, int index, const char ** element);

/**
 * The index of an atomic element given its name.
 *
 * @param physics     Handle for the Physics tables.
 * @param element     The element name.
 * @param index       The corresponding index.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The element index corresponds to the order of declaration specified in the
 * Material Description File (MDF).
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_PHYSICS_ERROR             The Physics is not initialised.
 *
 *     PUMAS_RETURN_UNKNOWN_MATERIAL          The material is not defined.
 */
PUMAS_API enum pumas_return pumas_physics_element_index(
    const struct pumas_physics * physics, const char * element, int * index);

/**
 * The name of a material given its index.
 *
 * @param physics    Handle for the Physics tables.
 * @param index      The material index.
 * @param material   The corresponding material name.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The material name is defined in the Material Description File (MDF).
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index is not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR             The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_material_name(
    const struct pumas_physics * physics, int index, const char ** material);

/**
 * The index of a material given its name.
 *
 * @param physics     Handle for the Physics tables.
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
 *     PUMAS_RETURN_PHYSICS_ERROR             The Physics is not initialised.
 *
 *     PUMAS_RETURN_UNKNOWN_MATERIAL          The material is not defined.
 */
PUMAS_API enum pumas_return pumas_physics_material_index(
    const struct pumas_physics * physics, const char * material, int * index);

/**
 * The total number of materials.
 *
 * @param physics    Handle for the Physics tables.
 * @return The total number of known materials, base plus composite.
 */
PUMAS_API int pumas_physics_material_length(
    const struct pumas_physics * physics);

/**
 * Get the properties of a material.
 *
 * @param physics           Handle for the Physics tables.
 * @param index             The material index.
 * @param length            The number of atomic elements.
 * @param density           The material reference density in kg/m^(3).
 * @param I                 The material mean excitation energy in GeV.
 * @param density_effect    The density effect parameters.
 * @param components        The vector of indices of the atomic elements.
 * @param fractions         The vector of mass fractions of the atomic elements.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Get the properties of a base material. `length`, `density`, `I`,
 * `density_effect`, `components` or `fractions` can be `NULL` in which case
 * the corresponding property is not retrieved.
 *
 * **Note**: the mean excitation energy, *I*, and the *density_effect*
 * parameters are not defined for composite materials resulting in a
 * `PUMAS_RETURN_INDEX_ERROR` if requested nevertheless.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index is not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR             The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_material_properties(
    const struct pumas_physics * physics, int index, int * length,
    double * density, double * I,
    struct pumas_physics_density_effect * density_effect, int * components,
    double * fractions);

/**
 * The number of composite materials.
 *
 * @param physics    Handle for the Physics tables.
 * @return The number of composite materials.
 */
PUMAS_API int pumas_physics_composite_length(
    const struct pumas_physics * physics);

/**
 * Update the properties of a composite material.
 *
 * @param physics    Handle for the Physics tables.
 * @param material   The composite material index.
 * @param fractions  The vector of mass fractions of the base materials
 *                   components.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Update the composition of a composite material. **Note** that negative
 * fraction values are treated as zero.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_VALUE_ERROR               The *fractions* pointer is `NULL`.
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index is not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR             The Physics is not initialised.
 *
 *     PUMAS_RETURN_MEMORY_ERROR              Couldn't allocate memory.
 */
PUMAS_API enum pumas_return pumas_physics_composite_update(
    struct pumas_physics * physics, int material, const double * fractions);

/**
 * Get the properties of a composite material.
 *
 * @param physics    Handle for the Physics tables.
 * @param index      The composite material index.
 * @param length     The number of base material components of the composite.
 * @param components The vector of indices of the base materials.
 * @param fractions  The vector of mass fractions of the base materials
 *                   components.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Get the properties of a composite material. `length`, `components` or
 * `fractions` can be `NULL` in which case the corresponding property is not
 * retrieved.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR               The provided index is not valid.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR             The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_composite_properties(
    const struct pumas_physics * physics, int index, int * length,
    int * components, double * fractions);

/**
 * Accessor to the tabulated Physics data.
 *
 * @param physics     Handle for the Physics tables.
 * @param property    The column index of a property of interest.
 * @param scheme      The energy loss scheme.
 * @param material    The material index.
 * @param row         The kinetic energy value row index in the table.
 * @param value       The corresponding table value.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * For a given `material` and energy loss `scheme`, this function returns the
 * tabulated data corresponding to the given `property` column and `row` index.
 * Each row of the table corresponds to a different kinetic energy value.
 *
 * **Note** that `PUMAS_PROPERTY_SCATTERING_LENGTH` is not supported since it
 * is not tabulated but computed on the fly.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             Some input index is not valid
 * (property, material or scheme).
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 */
PUMAS_API enum pumas_return pumas_physics_table_value(
    const struct pumas_physics * physics, enum pumas_property property,
    enum pumas_mode scheme, int material, int row, double * value);

/**
 * The depth, i.e. number of kinetic energy values, of the tabulated data.
 *
 * @param physics    Handle for the Physics tables.
 * @return The number of rows in data tables.
 */
PUMAS_API int pumas_physics_table_length(const struct pumas_physics * physics);

/**
 * Compute the table row index for a given property and its value.
 *
 * @param physics     Handle for the Physics tables.
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
 * **Note** that only monotone properties are supported, i.e. where there is
 * at most one solution. Those are: `PUMAS_PROPERTY_GRAMMAGE`,
 * `PUMAS_PROPERTY_KINETIC_ENERGY`, `PUMAS_PROPERTY_MAGNETIC_ROTATION` and
 * `PUMAS_PROPERTY_PROPER_TIME`.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             Some input index is not valid
 * (property, material or scheme).
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           The Physics is not initialised.
 *
 *     PUMAS_RETURN_VALUE_ERROR             The provided value is out of the
 * table.
 */
PUMAS_API enum pumas_return pumas_physics_table_index(
    const struct pumas_physics * physics, enum pumas_property property,
    enum pumas_mode scheme, int material, double value, int * index);

/**
 * Create a new particle recorder.
 *
 * @param recorder     A handle for the recorder.
 * @param extra_memory The size of the user extra memory, if any is claimed.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Create a new Monte-Carlo particle recorder. The recorder starts configured
 * with a built-in in-memory frame recorder.
 *
 * If `extra_memory` is strictly positive the recorder will be extended by
 * `extra_memory` bytes for user usage. This memory can then be accessed with
 * the `user_data` field of the returned `pumas_recorder` structure.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_MEMORY_ERROR    Couldn't allocate memory.
 */
PUMAS_API enum pumas_return pumas_recorder_create(
    struct pumas_recorder ** recorder, int extra_memory);

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
typedef void * pumas_allocate_cb (size_t size);

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
typedef void * pumas_reallocate_cb (void * ptr, size_t size);

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
typedef void pumas_deallocate_cb (void * ptr);

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

/** Physical states of materials. */
enum pumas_physics_state {
        /** Undefined physical state. */
        PUMAS_PHYSICS_STATE_UNKNOWN = 0,
        /** Solid physical state. */
        PUMAS_PHYSICS_STATE_SOLID,
        /** Liquid physical state. */
        PUMAS_PHYSICS_STATE_LIQUID,
        /** Gas physical state. */
        PUMAS_PHYSICS_STATE_GAS,
        /** The number of physical states. */
        PUMAS_PHYSICS_N_STATES
};

/**
 * Handle for an atomic element within a material.
 *
 * This structure is a proxy exposing some data of an atomic element within
 * the material last processed by `pumas_tabulation_tabulate`.
 */
struct pumas_physics_element {
        /** Linked list pointer to the previous element. */
        struct pumas_physics_element * prev;
        /** Linked list pointer to the next element. */
        struct pumas_physics_element * next;
        /** The element index. */
        int index;
        /** The mass fraction of the element in the current material. */
        double fraction;
};

/**
 * Description of a base material to tabulate.
 *
 * This structure allows to specify the physical properties of a base material
 * in order to tabulate its energy loss using the `pumas_physics_tabulate`
 * function.
 *
 * **Note** that the material index *must* be set, e.g. using the
 * `pumas_physics_material_index` function. If the Sternheimer coefficients are
 * not explicitly provided there are computed from the density using the
 * Sternheimer and Peierls recipe.
 */
struct pumas_physics_material {
        /** The material index. */
        int index;
        /** The mean excitation energy in GeV. */
        double I;
        /** The material state. */
        enum pumas_physics_state state;
        /** The density effect parameters */
        struct pumas_physics_density_effect density_effect;
};

/**
 * Handle for tabulation data.
 *
 * This structure gathers data related to the tabulation of the energy loss of
 * materials with the `pumas_physics_tabulate` function. **Note** that the two
 * last parameters: *path* and *elements* should not be set directly. They are
 * filled in (updated) by the `pumas_physics_tabulate` function. Note also
 * that if no energy grid is provided then a default one is set.
 *
 * **Warning**: the energy grid should not be changed between successive calls
 * to `pumas_physics_tabulate`. If a new energy grid is needed then a new
 * `pumas_physics_tabulation_data` object must be created.
 */
struct pumas_physics_tabulation_data {
        /** The number of kinetic energy values to tabulate. Providing a value
         * of zero or less results in a default energy grid being set.
         */
        int n_energies;
        /** Array of kinetic energy values to tabulate. Providing a `NULL`
         * value results in a default energy grid being set.
         */
        double * energy;
        /** Flag to enable overwriting an existing energy loss file. */
        int overwrite;
        /** Path to a directory where the tabulation should be written. */
        char * outdir;
        /** Properties of the material to tabulate */
        struct pumas_physics_material material;
        /** Path to the energy loss file of the last tabulated material. */
        char * path;
        /** List of atomic elements contained in the tabulated material(s). */
        struct pumas_physics_element * elements;
};

/**
 * Tabulate the energy loss for the given material and set of energies.
 *
 * @param physics    Handle for the Physics tables.
 * @param data       The tabulation settings.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * This function allows to generate an energy loss file for a given material and
 * a set of kinetic energy values. **Note** that the Physics must have been
 * initialised with the `pumas_physics_create_tabulation` function. The
 * material atomic composition is specified by the MDF provided at
 * initialisation. Additional Physical properties can be specified by filling
 * the input *data* structure.
 *
 * __Warnings__
 *
 * This function is **not** thread safe.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR     The material index is not valid.
 *
 *     PUMAS_RETURN_IO_ERROR        The output file already exists.
 *
 *     PUMAS_RETURN_MEMORY_ERROR    Some memory couldn't be allocated.
 *
 *     PUMAS_RETURN_PATH_ERROR      The output file could not be created.
 *
 */
PUMAS_API enum pumas_return pumas_physics_tabulate(
    struct pumas_physics * physics,
    struct pumas_physics_tabulation_data * data);

/**
 * Clear the temporary memory used for the tabulation of materials.
 *
 * @param data    The tabulation data.
 *
 * This function allows to clear any temporary memory allocated by the
 * `pumas_physics_tabulate` function.
 *
 * __Warnings__
 *
 * This function is **not** thread safe.
 */
PUMAS_API void pumas_physics_tabulation_clear(
    const struct pumas_physics * physics,
    struct pumas_physics_tabulation_data * data);

/**
 * Initialise the Physics in tabulation mode.
 *
 * @param physics      Handle for the Physics tables.
 * @param particle     The type of the particle to transport.
 * @param mdf_path     The path to a Material Description File (MDF) or `NULL`.
 * @param settings     Extra physics settings or `NULL`.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * Initialise the Physics tables in reduced mode using a MDF. The materials data
 * are not loaded. This mode is not suitable for particle's transport. It is
 * meant for pre-computation of the material tables using the
 * `pumas_physics_tabulate` function. **Note** that *mdf_path* can be `NULL`
 * Then it is read from the `PUMAS_MDF` environment variable.
 *
 * Optionaly extra physics settings can be specified by providing a
 * `pumas_physics_settings` structure. If `NULL` is provided then PUMAS default
 * physics settings are used which should perform well for most use cases.
 * **Note** that the *cutoff* parameter is not used in tabulation mode.
 *
 * Call `pumas_physics_destroy` in order to unload the Physics and release
 * allocated memory. **Note** that in addition any temporary memory allocated by
 * `pumas_physics_tabulate` must be explictly freed by calling the
 * `pumas_physics_tabulation_clear` function.
 *
 * **Warnings** : this function is not thread safe.
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
 *     PUMAS_RETURN_IO_ERROR                A file couldn't be read.
 *
 *     PUMAS_RETURN_MEMORY_ERROR            Couldn't allocate memory.
 *
 *     PUMAS_RETURN_MODEL_ERROR             A requested DCS model is not valid.
 *
 *     PUMAS_RETURN_PATH_ERROR              A file couldn't be opened.
 *
 *     PUMAS_RETURN_PHYSICS_ERROR           A `NULL` physics pointer was
 * provided.
 *
 *     PUMAS_RETURN_TOO_LONG                Some XML node in the MDF is
 * too long.
 *
 *     PUMAS_RETURN_UNDEFINED_MDF           No MDF was provided.
 *
 *     PUMAS_RETURN_UNKNOWN_ELEMENT         An element in the MDF wasn't
 * defined.
 *
 *     PUMAS_RETURN_UNKNOWN_MATERIAL        An material in the MDF wasn't
 * defined.
 *
 *     PUMAS_RETURN_UNKNOWN_PARTICLE        The given type is not supported.
 */
PUMAS_API enum pumas_return pumas_physics_create_tabulation(
    struct pumas_physics ** physics, enum pumas_particle particle,
    const char * mdf_path, const struct pumas_physics_settings * settings);

/**
 * Get the physics Differential Cross-Section (DCS) for a given process.
 *
 * @param physics      Handle for the Physics tables or `NULL`.
 * @param process      The Physics process.
 * @param model        The corresponding DCS model.
 * @param dcs          The corresponding DCS function.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * The *dcs* and *model* return values are optionnal. If a `NULL` pointer is
 * provided then the corresponding return value is not filled.
 *
 * **Note** : the DCS model is set at the physics creation and cannot be
 * changed afterwards.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The *process* index is not a valid.
 *
 */
PUMAS_API enum pumas_return pumas_physics_dcs(
    const struct pumas_physics * physics, enum pumas_process process,
    const char ** model, pumas_dcs_t ** dcs);

/**
 * Get a PUMAS library constant.
 *
 * @param index     The constant index.
 * @param value     The corresponding value.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The index is not a valid.
 *
 *     PUMAS_RETURN_VALUE_ERROR             The value pointer is `NULL`.
 */
PUMAS_API enum pumas_return pumas_constant(
    enum pumas_constant index, double * value);

/**
 * Register a Differential Cross Section (DCS) model to PUMAS.
 *
 * @param process   The physics process index.
 * @param model     The model name.
 * @param dcs       The model DCS.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * This function allows to register a DCS model for a physics interaction. Note
 * the model name must not be already used otherwise an error is returned.  In
 * addition, only the following radiative models can be user defined:
 * Bremsstrahlung, e^(+)e^(-) pair production or photonuclear interactions.
 * Ionisation loss & Coulomb multiple scattering are built-in.
 *
 * **Note**: it is not possible to un-register a model.
 *
 * **Warnings** : this function is not thread safe.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The process index is not a valid.
 *
 *     PUMAS_RETURN_MEMORY_ERROR            The maximum number of models was
 * reached.
 *
 *     PUMAS_RETURN_MODEL_ERROR             The model name is already used.
 *
 *     PUMAS_RETURN_VALUE_ERROR             A `NULL` model or dcs was provided.
 */
PUMAS_API enum pumas_return pumas_dcs_register(
    enum pumas_process process, const char * model, pumas_dcs_t * dcs);

/**
 * Get the Differential Cross Section (DCS) for a given model.
 *
 * @param process   The physics process index.
 * @param model     The model name.
 * @param dcs       The corresponding DCS.
 * @return On success `PUMAS_RETURN_SUCCESS` is returned otherwise an error
 * code is returned as detailed below.
 *
 * This function allows to retrieve the DCS for a given physics process and
 * model. See the `struct pumas_physics_settings` for a list of models available
 * by default. Extra models can be registered with the `pumas_dcs_register`
 * function.
 *
 * **Warnings** : this function is not thread safe.
 *
 * __Error codes__
 *
 *     PUMAS_RETURN_INDEX_ERROR             The process index is not a valid.
 *
 *     PUMAS_RETURN_MODEL_ERROR             The model name is not valid.
 */
PUMAS_API enum pumas_return pumas_dcs_get(
    enum pumas_process process, const char * model, pumas_dcs_t ** dcs);

/**
 * Get the default Differential Cross Section (DCS) model for a given process.
 *
 * @param process   The physics process index.
 * @return On success the default model name is returned otherwise `NULL`.
 *
 * Get the name of the default DCS model for a given process.  If the requested
 * process index is not valid then `NULL` is returned.
 */
PUMAS_API const char * pumas_dcs_default(enum pumas_process process);

#ifdef __cplusplus
}
#endif
#endif
