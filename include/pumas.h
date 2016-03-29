/*
 * C API for PUMAS: Semi-Analytical Propagation of MUons with forward or
 * backward Monte-Carlo.
 */
#ifndef pumas_h
#define pumas_h
#ifdef __cplusplus
extern "C" {
#endif

/* For C standard streams. */
#ifndef FILE
#include <stdio.h>
#endif

/**
 * Key tabulated properties used by PUMAS.
 */
enum pumas_property {
	/** The macroscopic inelastic cross-section, in m^(2)/kg. */
	PUMAS_PROPERTY_CROSS_SECTION = 0,
	/** The average energy loss, in GeV/(kg/m^(2)). */
	PUMAS_PROPERTY_ENERGY_LOSS,
	/** The muon range, in kg/m^(2). */
	PUMAS_PROPERTY_GRAMMAGE,
	/** The muon kinetic energy, in GeV. */
	PUMAS_PROPERTY_KINETIC_ENERGY,
	/** The total magnetic rotation angle, in radians/(kg/m^(3)). */
	PUMAS_PROPERTY_MAGNETIC_ROTATION,
	/** The muon proper time, in kg/m^(2). */
	PUMAS_PROPERTY_PROPER_TIME,
	/** The macroscopic elastic scattering 1^(st) path length, in kg/m^(2).
	 */
	PUMAS_PROPERTY_SCATTERING_LENGTH
};

/**
 * Monte-Carlo schemes for the computation of energy losses.
 */
enum pumas_scheme {
	/** All energy losses are disabled.
	 *
	 * **Note** : This mode is provided for test purpose only. Running
	 * without energy losses requires specifying a geometry through a
	 * `pumas_locator_cb` callback.
	 */
	PUMAS_SCHEME_NO_LOSS = -1,
	/** Energy losses are purely continuous, as given by the Continuously
	 * Slowing Down Approximation (CSDA).
	 */
	PUMAS_SCHEME_CSDA,
	/** Energy losses are simulated using an hybrid Monte-Carlo scheme with
	 * hard stochastic interactions and soft continuous energy loss.
	 */
	PUMAS_SCHEME_HYBRID,
	/** In addition to the hybrid scheme small fluctuations of the soft
	 * energy loss are also simulated.
	 */
	PUMAS_SCHEME_DETAILED
};

/** Return codes for the API functions. */
enum pumas_return {
	/*! An error occured. */
	PUMAS_ERROR = -1,
	/*! Execution was successful. */
	PUMAS_SUCCESS = 0
};

/**
 * Container for a Muon Monte-Carlo state.
 */
struct pumas_state {
	/** The particle's electric charge. Note that non physical values,
	 * i.e. different from 1 or -1, can be set if needed. */
	double charge;
	/** The kinetic energy, in GeV. */
	double kinetic;
	/** The total travelled distance, in m. */
	double distance;
	/** The total travelled grammage, in kg/m<sup>2</sup>. */
	double grammage;
	/** The particle's proper time, in m/c. */
	double time;
	/** The particle's Monte-Carlo weight. */
	double weight;
	/** The particle's absolute location, in m. */
	double position[3];
	/** The particle's momentum direction. */
	double direction[3];
};

/**
 * The local properties of a propagation medium.
 */
struct pumas_locals {
	/*! The material's local density, in kg/m^3. */
	double density;
	/*! The local magnetic field components, in T. */
	double magnet[3];
};

/**
 * Callback for setting the local properties of a propagation medium.
 *
 * @param state The muon Monte-Carlo state for which the local properties
 * are requested.
 * @param locals A pointer to a `pumas_locals` structure to update.
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
typedef double (pumas_locals_cb)(const struct pumas_state * state,
	struct pumas_locals * locals);

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
	pumas_locals_cb * setter;
};

/** A handle to a recorded Monte-Carlo frame.
 *
 *  This structure exposes data relative to a recorded frame. It is not meant
 * to be modified by the user.
 */
struct pumas_frame {
	/** The recorded muon state. */
	struct pumas_state state;
	/** The corresponding propagation medium. */
	int medium;
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

/*!
 * Callback for locating the propagation medium of a `pumas_state`.
 *
 * @param state The muon Monte-Carlo state for which the local properties
 * are requested.
 * @return The corresponding medium index or a negative number if the state has
 * exit the simulation area.
 */
typedef int (pumas_locator_cb)(const struct pumas_state * state);

struct pumas_context;
/*!
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
typedef double (pumas_random_cb)(struct pumas_context * context);

/*!
 * A handle for a simulation stream.
 *
 * This structure is a proxy to thread specific data for a simulation stream.
 * It exposes some public data that the user may configure or alter directly.
 * However, it also encloses other opaque data. Therefore, it **must** be
 * initialised and released with the `pumas_context` functions.
 *
 * + The `media` field must be set after any initialisation with
 * `pumas_context_create` and prior to any call to `pumas_propagate`. At least
 * one medium must be registered. There is no soft limitation to the number of
 * media.
 * + The `locator` callback is optional. Setting it to `NULL` implies that the
 * simulation area is composed of a single medium, `*media`, uniform and of
 * infinite extension.
 * + Depending on the level of detail of the simulation a random stream must
 * be provided by the user before any call to `pumas_propagate`.
 * + For `kinetic_limit`, `distance_max` or `grammage_max` a strictly positive
 * value activates the corresponding limitation. Setting it to `0` or less
 * disables it however.
 */
struct pumas_context {
	/*! An array of possible propagation media. */
	const struct pumas_medium * media;
	/*! A medium locator callback. */
	pumas_locator_cb * locator;
	/*! The pseudo random generator callback. */
	pumas_random_cb * random;
	/*! A `pumas_frame` recorder. */
	struct pumas_recorder * recorder;
	/*! A pointer to additional memory, if any is requested at
	 * initialisation.
	 */
	void * user_data;

	/*! Flag to enable or disable transverse transportation. */
	int longitudinal;
	/*! Flag to switch between forward and reverse propagation. */
	int forward;
	/*! The scheme used for the computation of energy losses. */
	enum pumas_scheme scheme;

	/*! The minimum kinetic energy for forward propagation, or maximum one
	 * for backward propagation.
	 */
	double kinetic_limit;
	/*! The maximum propagation distance. */
	double distance_max;
	/*! The maximum propagation grammage. */
	double grammage_max;
};

/*!
 * A handle to a structured error summary.
 *
 * This is a proxy to public error data exposed to the user. It also encloses
 * additional opaque data. Therefore, it must be handled with the `pumas_error`
 * functions.
 */
struct pumas_error {
	/*! The line in the source code where the error occurred. */
	int line;
	/*! Standard [error number](http://en.wikipedia.org/wiki/Errno.h) if
	 * provided.
	 */
	int errnum;
	/*! The name of the function where the error occurred. */
	char * function;
	/*! Any additional error message. It might be `NULL`. */
	char *message;
	/*! A link to the next error in the stack. */
	struct pumas_error* next;
};

/*!
 * Initialise the PUMAS library.
 *
 * @param path The path to a Material Description File (MDF).
 * @param error_size The size of the error buffer.
 * @return `PUMAS_SUCCESS` on success, `PUMAS_ERROR` otherwise.
 *
 * Initialise the library from a MDF. Load the materials data and precompute
 * various properties. If `path` is `NULL` the MDF is read from the `PUMAS_MDF`
 * environment variable. If `error_size` is less or equal to zero no errors
 * details will be logged, though `pumas_error_count` will still report if any
 * error occured. Call `pumas_finalise` in order to unload the library and
 * release all alocated memory.
 *
 * **Warnings** : this function is not thread safe. Trying to (re-)initialise an
 * already initialised library will generate an error. `pumas_finalise` must
 * be called first.
 */
enum pumas_return pumas_initialise(const char * path, int error_size);

/*!
 * Finalise the PUMAS library.
 *
 * @param stream An output stream where to log errors.
 *
 * Finalise the library. Free the shared memory and dump the error stack. Call
 * `pumas_initialise` in order to reload the library. The argument `stream` can
 * be `NULL` in which case no error will be reported, though the error stack
 * will still be cleared.
 *
 * **Warnings** : This function is not thread safe. Finalising the library
 * doesn't release the memory allocated for any `pumas_context`.
 *
 * @see pumas_initialise
 */
void pumas_finalise(FILE * stream);

/*!
 * Propagate a muon according to the configured `pumas_context`.
 *
 * @param context The simulation context.
 * @param state The muon initial and final state.
 * @return `PUMAS_SUCCESS` on success or `PUMAS_ERROR` otherwise.
 *
 * Depending on the `context` configuration the muon is propagated through a
 * single infinite `pumas_medium` or using a `pumas_locator_cb` function. At
 * return, the muon `state` is updated.
 */
enum pumas_return pumas_propagate(struct pumas_context * context,
	struct pumas_state * state);

/*!
 * Write a summary of the current library configuration.
 *
 * @param stream A stream where the summary will be formated to.
 * @return `PUMAS_SUCCESS` on success or `PUMAS_ERROR` otherwise.
 *
 * The summary provides information on loaded materials as well as some basic
 * statistics.
 */
enum pumas_return pumas_inform(FILE * stream);

/*!
 * Get the hash identifier of the library.
 *
 * @return The library hash encoded on an `int`.
 *
 * The library hash is given by the short version of the source code git's hash.
 */
int pumas_hash();

/*!
 * Create a simulation context.
 *
 * @param error_size The size of the error buffer.
 * @param extended_memory The size of the user extended memory, if any is
 * claimed.
 * @return On success, a handle to the simulation context is returned.
 * Otherwise `NULL` is returned.
 *
 * Create a new simulation context with a default configuration. If
 * `error_size` is less or equal to zero no errors will be logged, though
 * `pumas_error_count` will still report a positive number if any error occured.
 * Call `pumas_context_destroy` in order to release all the memory allocated for
 * the context.
 *
 * If `extended_memory` is strictly positive the context will be extended by
 * `extended_memory` bytes for user usage. This memory can then be accessed with
 * the `user_data` field of the returned `pumas_context` structure.
 */
struct pumas_context * pumas_context_create(int error_size,
	int extended_memory);

/*!
 * Destroy a simulation context.
 *
 * @param context The simulation context.
 * @param stream An output stream where to log errors.
 * @return On success `context` is set to `NULL`.
 *
 * Call on a previously created context with `pumas_context_create` in order to
 * release the corresponding dynamicaly allocated memory. The argument `stream`
 * can be `NULL` in which case no error will be reported, though
 * `pumas_error_count` will still report if any error occured.
 */
void pumas_context_destroy(struct pumas_context ** context,
	FILE * stream);

/*!
 * @brief Compute the total grammage that a muon can travel assuming continuous
 * energy loss.
 *
 * @param scheme The energy loss scheme: `PUMAS_SCHEME_CSDA` or
 * `PUMAS_SCHEME_HYBRID`.
 * @param material The material index.
 * @param kinetic The muon initial kinetic energy, in GeV.
 * @return The grammage in kg/m^(2) or a negative value if an error occured.
 *
 * For a uniform medium, divide the return value by the density in order to get
 * the corresponding total travelled distance.
 */
double pumas_property_grammage(enum pumas_scheme scheme, int material,
	double kinetic);

/*!
 * Compute the normalised total proper time spent by a muon assuming continuous
 * energy loss.
 *
 * @param scheme The energy loss scheme: `PUMAS_SCHEME_CSDA` or
 * `PUMAS_SCHEME_HYBRID`.
 * @param material The material index.
 * @param kinetic The muon initial kinetic energy, in GeV.
 * @return The normalised proper time in kg/m^(2) or a negative value if an
 * error occured.
 *
 * Divide the returned value by the medium density times *c* in order to get the
 * proper time in unit of time.
 */
double pumas_property_proper_time(enum pumas_scheme scheme, int material,
	double kinetic);

/*!
 * Compute the normalised rotation angle due to a uniform magnetic field for
 * a CSDA muon.
 *
 * @param material The material index.
 * @param kinetic The muon initial kinetic energy, in GeV.
 * @return The normalised rotation angle in kg/m^(2)/T.
 *
 * Multiply the returned value by the transverse magnetic field amplitude and
 * divide by the density in order to get the rotation angle in radian.
 */
double pumas_property_magnetic_rotation(int material, double kinetic);

/*!
 * Compute the minimum kinetic energy required for travelling over a given
 * `grammage`, assuming continuous energy loss.
 *
 * @param scheme The energy loss scheme: `PUMAS_SCHEME_CSDA` or
 * `PUMAS_SCHEME_HYBRID`.
 * @param material The material index.
 * @param grammage The requested grammage, in kg/m^(2).
 * @return The kinetic energy in GeV or a negative value if an error occured.
 */
double pumas_property_kinetic_energy(enum pumas_scheme scheme, int material,
	double grammage);

/*!
 * Compute the average energy loss per unit weight of material.
 *
 * @param scheme The energy loss scheme
 * @param material The material index.
 * @param kinetic The muon kinetic energy, in GeV.
 * @return The energy loss in GeV/(kg/m^(2)) or a negative value if an error
 * occured.
 */
double pumas_property_energy_loss(enum pumas_scheme scheme, int material,
	double kinetic);

/*!
 * Compute the elastic scattering 1^(st) transport path length for a unit
 * weight.
 *
 * @param material The material index.
 * @param kinetic The muon kinetic energy, in GeV.
 * @return The path length in kg/m^(2).
 *
 * The elastic 1^(st) transport path length, &lambda;, is related to
 * the plan multiple scattering angle's standard deviation as &theta;^(2) =
 * 2&lambda;.
 */
double pumas_property_scattering_length(int material, double kinetic);

/*!
 * The macroscopic inelastic total cross section.
 *
 * @param material The material index.
 * @param kinetic The muon kinetic energy, in GeV.
 * @return the normalised cross section in m^(2)/kg or a negative value if an
 * error occured.
 *
 * Multiply the return value by the density in order to get the inverse of the
 * interaction length in unit of distance.
 */
double pumas_property_cross_section(int material, double kinetic);

/*!
 * The library muon mass.
 *
 * @return The mass value in GeV/c^(2).
 *
 * The library muon mass is a hard coded read only value.
 */
double pumas_muon_mass(void);

/*!
 * The library muon proper lifetime.
 *
 * @return The lifetime, in m.
 *
 * Divide the return value by *c* in order to get the proper lifetime in time
 * unit. The library muon lifetime is a hard coded read only value.
 */
double pumas_muon_lifetime(void);

/*!
 * The name of a material.
 *
 * @param material The material index.
 * @return A `C string` with the material name or `NULL` in case of a bad
 * index value.
 *
 * The material name is defined in the MDF.
 * **Note** : if the material index is out of bounds a library error is
 * registered.
 */
const char* pumas_material_name(int material);

/*!
 * The index of a material.
 *
 * @param material The material name.
 * @return On success the material index is returned, otherwise `PUMAS_ERROR` is
 * returned.
 *
 * The material index corresponds to the order of declaration specified in the
 * MDF.
 * **Note** : if the material is not found an error is registered.
 */
int pumas_material_index(const char * material);

/*!
 * The total number of materials.
 *
 * @return The total number of known materials, basic plus composite.
 */
int pumas_material_length(void);

/*!
 * The number of composite materials.
 *
 * @return The number of composite materials.
 */
int pumas_composite_length(void);

/*!
 * Update the properties of a composite material.
 *
 * @param material The composite material index.
 * @param fractions The vector of mass fractions of the base materials
 * components.
 * @param densities The vector of densities of the base materials components.
 * @return `PUMAS_ERROR` if the update failled, `PUMAS_SUCCESS` otherwise.
 *
 * Update the composition and/or the density of a composite material.
 * `fractions` or `densities` can be `NULL` in which case the corresponding
 * property is not updated.
 */
enum pumas_return pumas_composite_update(int material, const double * fractions,
	const double * densities);

/*!
 * Get the properties of a composite material.
 *
 * @param material The composite material index.
 * @param density The composite material reference density.
 * @param components The number of base material components of the composite.
 * @param fractions The vector of mass fractions of the base materials
 * components.
 * @param densities The vector of densities of the base materials components.
 * @return `PUMAS_ERROR` if the update failled, `PUMAS_SUCCESS` otherwise.
 *
 * Get the properties of a composite material. `density`, `components`,
 * `fractions` or `densities` can be `NULL` in which case the corresponding
 * property is not retrieved.
 */
enum pumas_return pumas_composite_properties(int material, double * density,
	int * components, double * fractions, double * densities);

/*!
 * Accessor to the tabulated shared data.
 *
 * @param property The column index of a property of interest.
 * @param scheme The energy loss scheme.
 * @param material The material index.
 * @param row The kinetic value row index in the table.
 * @return the table value or a negative value if an error occured.
 *
 * For a given `material` and energy loss `scheme`, this function returns the
 * tabulated data corresponding to the given `property` column and `row` index.
 * Each row of the table corresponds to a different kinetic energy value.
 */
double pumas_table_value(enum pumas_property property, enum pumas_scheme scheme,
	int material, int row);

/*!
 * The depth, i.e. number of kinetic values, of the tabulated data.
 *
 * @return The number of rows in data tables.
 */
int pumas_table_length(void);

/*!
 * Compute the the table row index for a given property and its value.
 *
 * @param property The column index of the property.
 * @param scheme The energy loss scheme.
 * @param material The material index.
 * @param value The property value.
 * @return The row index from below for the given value or a negative value if
 * an error occured.
 *
 * In the case of an out of bounds value the closest index value is
 * returned. If any of the other parameters `property`, `scheme` or `material`
 * takes a bad value, a negative index is returned.
 */
int pumas_table_index(enum pumas_property property, enum pumas_scheme scheme,
	int material, double value);

/*!
 * Get the current count of occurred errors.
 *
 * @param context The simulation context.
 * @return The number of occured errors.
 *
 * If `context` is `NULL` the count of library errors is returned. Otherwise,
 * the count of errors related to the given context is returned.
 */
int pumas_error_count(struct pumas_context * context);

/*!
 * Get the first error in the stack.
 *
 * @param context The simulation context.
 * @return If any error occcured a pointer to the first `pumas_error` is
 * returned, otherwise `NULL` is returned.
 *
 * If `context` is `NULL` the error is taken from the libraries' error stack.
 * Otherwise, the error is taken from the given context's error stack.
 */
struct pumas_error * pumas_error_get(struct pumas_context * context);

/*!
 * Clear all errors.
 *
 * @param context The simulation context.
 *
 * Remove all errors from the stack and reset the error counter. If `context` is
 * `NULL` the libraries' error stack is cleared. Otherwise, the given context's
 * error stack is cleared.
 */
void pumas_error_clear(struct pumas_context * context);

/*!
 * Format an error summary to a stream.
 *
 * @param stream The stream where to format the error.
 * @param error The error summary.
 */
void pumas_error_format(FILE* stream, struct pumas_error * error);

/*!
 * Create a new muon recorder.
 *
 * @param context The simulation context.
 * @return On success the new recorder is returned, otherwise `NULL` is
 * returned.
 *
 * The recorder is configured for a fresh start with default settings. It is
 * initially linked to the provided context.
 *
 * **Note** : though the recorder's creation requires a simulation context, the
 * returned proxy might be moved to any other context after creation. However
 * it should not be linked simultaneously to multiple concurent simulation
 * streams.
 */
struct pumas_recorder * pumas_recorder_create(struct pumas_context * context);

/*!
 * Clear all recorded frames.
 *
 * @param recorder The recorder handle.
 *
 * Erase all recorded states from the recorder and reset the frame count.
 */
void pumas_recorder_clear(struct pumas_recorder * recorder);

/*!
 * Destroy a muon recorder releasing all associated memory.
 *
 * @param recorder The recorder handle.
 *
 * **Note** : The recorder is cleared before beeing destroyed. At return
 * `recorder` is set to `NULL`.
 */
void pumas_recorder_destroy(struct pumas_recorder ** recorder);

#ifdef __cplusplus
}
#endif
#endif
