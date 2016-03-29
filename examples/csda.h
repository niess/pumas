#ifndef CSDA_PROPAGATE_H_

/* Configurable settings for the CSDA propagation. */
struct csda_settings {
	/* The initial muon kinetic energy in GeV. */
	double kinetic_energy;

	/* The propagation kinetic threshold in GeV. */
	double kinetic_limit;
	/* The propagation distance in m. */
	double distance_max;

	/* The stepping limitation, if any. */
	double step_max;

	/* The magnetic field amplitude, in T. */
	double magnet_amplitude;
	/* The magnetic angle with respect to the muon initial direction,
	 * in degrees. */
	double magnet_angle;

	/* The medium uniform density in kg/m^3 */
	double density;
	/* The name of the medium material. */
	const char * material;
};

/* Pointer to the propagation settings. */
struct csda_settings * csda_settings(void);

/* Configure a simulation context for CSDA propagation. */
enum pumas_return csda_configure(struct pumas_context * context);

/* Propagate PuMAS a state with CSDA. */
enum pumas_return csda_propagate(struct pumas_context * context,
	struct pumas_state * state);

/* Dump a PumAS state to a stream. */
void csda_dump(FILE * stream, struct pumas_state * state);

#endif /* CSDA_PROPAGATE_H_ */
