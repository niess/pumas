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

#ifndef pumas_tabulate_h
#define pumas_tabulate_h

#ifndef PUMAS_API
#define PUMAS_API
#endif

/* The PUMAS API. */
#include "pumas.h"

/* Matter states, for density effect computation. */
enum tabulation_state {
        TABULATION_STATE_UNKNOWN = 0,
        TABULATION_STATE_SOLID,
        TABULATION_STATE_LIQUID,
        TABULATION_STATE_GAZ,
        TBULATION_N_STATES
};

/* Storage for per element data. */
struct tabulation_element {
        /* Linked list pointers. */
        struct tabulation_element * prev;
        struct tabulation_element * next;
        /* The element index. */
        int index;
        /* Mass fraction in the current material. */
        double fraction;
        /* Placeholder for data. */
        double data[];
};

/* Storage for a base material to process. */
struct tabulation_material {
        /* Linked list pointers. */
        struct tabulation_material * prev;
        struct tabulation_material * next;
        /* The material index. */
        int index;
        /* The material density. */
        double density;
        /* The mean excitation energy. */
        double I;
        /* The material state. */
        enum tabulation_state state;
        /* Sternheimer Coefficients. */
        double a, k, x_0, x_1, Cbar, delta0;
        /* The resulting density effect. */
        double delta;
};

/* Handle for managing tabulation data. */
struct tabulation_data {
        int n_kinetics;
        double * kinetic;
        char * path;
        char * outdir;
        int overwrite;
        struct tabulation_element * e_stack;
        struct tabulation_material * m_stack;
};

/*
 * Dry initialisation for the tabulations. The kinetic energy tables are not
 * loaded.
 */
PUMAS_API enum pumas_return _pumas_initialise_dry(
    enum pumas_particle particle, const char * mdf_path,
    struct pumas_error * error);

/* Clear the temporary data for the tabulations. */
void tabulation_finalise(void);

/* Tabulate the energy loss for the base material on top of the stack. */
PUMAS_API enum pumas_return _pumas_tabulate(struct tabulation_data * data);

/* Manage the materials temporary data for the tabulations. */
struct tabulation_material * tabulation_material_create(
    struct tabulation_data * data, int material);
void tabulation_material_drop(struct tabulation_data * data);
struct tabulation_material * tabulation_material_get(
    struct tabulation_data * data, int material);

#endif
