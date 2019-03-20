/*
 * Copyright (C) 2019 Université Clermont Auvergne, CNRS/IN2P3, LPC
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

/*
 * Unit tests for the PUMAS C library
 */

/* C89 standard library */
#include <stdint.h>
#include <stdlib.h>
/* The Check library */
#include "check.h"
/* The PUMAS library */
#include "pumas.h"


/* Media densities and geometric parameters */
#define TEST_ROCK_DENSITY 2.65E+03
#define TEST_ROCK_DEPTH 1E+02
#define TEST_AIR_DENSITY 1.205E+00
#define TEST_MAX_ALTITUDE 1E+04


/* Global handles for common PUMAS objects */
static struct pumas_context * context = NULL;
static struct pumas_recorder * recorder = NULL;
static struct pumas_state state_data, * state = &state_data;


/* Error handling for the test cases */
struct {
        enum pumas_return rc;
        pumas_function_t * caller;
        char message[2048];
} error_data = { PUMAS_RETURN_SUCCESS, NULL, "" };

static void handle_error(
    enum pumas_return rc, pumas_function_t * caller, const char * message)
{
        error_data.rc = rc;
        error_data.caller = caller;
        strcpy(error_data.message, message);
}

static void reset_error(void)
{
        error_data.rc = PUMAS_RETURN_SUCCESS;
        error_data.caller = NULL;
        error_data.message[0] = 0x0;
}


/* PRNG for test cases */
static double uniform01(struct pumas_context * context)
{
        return (1. + rand()) / (RAND_MAX + 2.);
}


/* Test the error API */
START_TEST (test_api_error)
{
        /* Check the initialisation */
        pumas_handler_cb * handler = pumas_error_handler_get();
        ck_assert_ptr_null(handler);

        /* Check the setters & getters  */
        pumas_error_handler_set(&handle_error);
        handler = pumas_error_handler_get();
        ck_assert_ptr_eq(handler, &handle_error);

        /* Check the error catch & raise */
        reset_error();
        pumas_error_catch(1);
        pumas_material_index(NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(error_data.caller);

        enum pumas_return rc = pumas_error_raise();
        ck_assert_int_eq(rc, PUMAS_RETURN_INITIALISATION_ERROR);
        ck_assert_int_eq(error_data.rc, rc);
        ck_assert_ptr_eq(error_data.caller, &pumas_material_index);

        /* Check that errors are enabled again */
        reset_error();
        pumas_material_index(NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        /* Check that error catching can be disabled */
        reset_error();
        pumas_error_catch(1);
        pumas_material_index(NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        pumas_error_catch(0);
        pumas_material_index(NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        /* Check the stringification of API functions */
#define CHECK_STRING(function)                                                 \
        ck_assert_str_eq(                                                      \
            pumas_error_function((pumas_function_t *)&function), #function)
        /* Generated with:
         * nm -D --defined-only libpumas.so |                                  \
         *     grep "T " |                                                     \
         *     awk '{print "CHECK_STRING("$3");"}'
         */
        CHECK_STRING(pumas_composite_length);
        CHECK_STRING(pumas_composite_properties);
        CHECK_STRING(pumas_composite_update);
        CHECK_STRING(pumas_context_create);
        CHECK_STRING(pumas_context_destroy);
        CHECK_STRING(pumas_dump);
        CHECK_STRING(pumas_error_catch);
        CHECK_STRING(pumas_error_function);
        CHECK_STRING(pumas_error_handler_get);
        CHECK_STRING(pumas_error_handler_set);
        CHECK_STRING(pumas_error_raise);
        CHECK_STRING(pumas_finalise);
        CHECK_STRING(pumas_initialise);
        CHECK_STRING(pumas_load);
        CHECK_STRING(pumas_material_index);
        CHECK_STRING(pumas_material_length);
        CHECK_STRING(pumas_material_name);
        CHECK_STRING(pumas_memory_allocator);
        CHECK_STRING(pumas_memory_deallocator);
        CHECK_STRING(pumas_memory_reallocator);
        CHECK_STRING(pumas_particle);
        CHECK_STRING(pumas_print);
        CHECK_STRING(pumas_property_cross_section);
        CHECK_STRING(pumas_property_energy_loss);
        CHECK_STRING(pumas_property_grammage);
        CHECK_STRING(pumas_property_kinetic_energy);
        CHECK_STRING(pumas_property_magnetic_rotation);
        CHECK_STRING(pumas_property_proper_time);
        CHECK_STRING(pumas_property_scattering_length);
        CHECK_STRING(pumas_recorder_clear);
        CHECK_STRING(pumas_recorder_create);
        CHECK_STRING(pumas_recorder_destroy);
        CHECK_STRING(pumas_table_index);
        CHECK_STRING(pumas_table_length);
        CHECK_STRING(pumas_table_value);
        CHECK_STRING(pumas_tag);
        CHECK_STRING(pumas_transport);
#undef CHECK_STRING

        /* Check the NULL case */
        ck_assert_ptr_null(pumas_error_function(
            (pumas_function_t *)&test_api_error));
}
END_TEST


/* Test the version tag */
START_TEST (test_api_tag)
{
        const int tag = pumas_tag();
        ck_assert_int_eq(tag, 13);
}
END_TEST


static void load_muon(void)
{
#define TEST_MUON_DUMP ".pumas.muon"
        FILE * fid = fopen(TEST_MUON_DUMP, "rb");
        pumas_load(fid);
        fclose(fid);
}

static void dump_muon(void)
{
        FILE * fid = fopen(TEST_MUON_DUMP, "wb+");
        pumas_dump(fid);
        fclose(fid);
}

static void load_tau(void)
{
#define TEST_TAU_DUMP ".pumas.tau"
        FILE * fid = fopen(TEST_TAU_DUMP, "rb");
        pumas_load(fid);
        fclose(fid);
}

static void dump_tau(void)
{
        FILE * fid = fopen(TEST_TAU_DUMP, "wb+");
        pumas_dump(fid);
        fclose(fid);
}


/* Test the library initialisation & finalisation */
START_TEST (test_api_init)
{
        /* Check the particle selection */
        reset_error();
        pumas_initialise(-1, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_UNKNOWN_PARTICLE);

        /* Check the path error */
        reset_error();
        pumas_initialise(PUMAS_PARTICLE_MUON, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_UNDEFINED_DEDX);

        /* Check the load error */
        reset_error();
        pumas_load(NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PATH_ERROR);

        /* Check the initialisation and dump */
        reset_error();
        pumas_initialise(PUMAS_PARTICLE_MUON, "materials/mdf/standard.xml",
            "materials/dedx/muon");
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        reset_error();
        dump_muon();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        pumas_finalise();

        reset_error();
        pumas_initialise(PUMAS_PARTICLE_TAU, "materials/mdf/wet-rock.xml",
            "materials/dedx/tau");
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        reset_error();
        dump_tau();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        /* Check the re-initialisation error  */
        reset_error();
        pumas_initialise(PUMAS_PARTICLE_MUON, "materials/mdf/standard.xml",
            "materials/dedx/muon");
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        load_tau();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);
        pumas_finalise();

        /* Check the load */
        reset_error();
        load_tau();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        enum pumas_particle particle;
        double lifetime, mass;

        reset_error();
        pumas_particle(NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        pumas_particle(&particle, &lifetime, &mass);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(particle, PUMAS_PARTICLE_TAU);
        ck_assert_double_eq(lifetime, 87.03E-06);
        ck_assert_double_eq(mass, 1.77682);

        pumas_finalise();
        pumas_particle(NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        load_muon();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        reset_error();
        pumas_particle(&particle, &lifetime, &mass);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(particle, PUMAS_PARTICLE_MUON);
        ck_assert_double_eq(lifetime, 658.654);
        ck_assert_double_eq(mass, 0.10565839);

        /* Check the dump errors */
        reset_error();
        pumas_dump(NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PATH_ERROR);
        pumas_finalise();

        reset_error();
        pumas_dump(NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        /* Check the transport error */
        reset_error();
        pumas_transport(NULL, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);
}
END_TEST


static int dyn_mem_count = 0;

static void * test_malloc(size_t size)
{
        dyn_mem_count++;
        return malloc(size);
}

static void * test_realloc(void * ptr, size_t size)
{
        if (ptr == NULL)
                dyn_mem_count++;
        return realloc(ptr, size);
}

static void test_free(void * ptr)
{
        if (ptr != NULL)
                dyn_mem_count--;
        free(ptr);
}


/* Test the memory API */
START_TEST (test_api_memory)
{
        /* Set the test memory routines */
        pumas_memory_allocator(&test_malloc);
        pumas_memory_reallocator(&test_realloc);
        pumas_memory_deallocator(&test_free);

        /* Load the muon data */
        load_muon();
        ck_assert_int_gt(dyn_mem_count, 0);

        /* Release the memory */
        pumas_finalise();
        ck_assert_int_eq(dyn_mem_count, 0);

        /* Set the system memory routines */
        pumas_memory_allocator(NULL);
        pumas_memory_reallocator(NULL);
        pumas_memory_deallocator(NULL);
}
END_TEST


/* Test the materials API */
START_TEST (test_api_material)
{
        int i, index;
        const char * name;
        const char * names[] = { "StandardRock", "Water", "Air" };

        /* Check the initialisation error */
        reset_error();
        pumas_material_index("Strawberry", &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_material_name(0, &name);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        index = pumas_material_length();
        ck_assert_int_eq(index, 0);

        /* Load the muon data */
        load_muon();

        /* Check the number of materials and composites */
        index = pumas_material_length();
        ck_assert_int_eq(index, 3);

        index = pumas_composite_length();
        ck_assert_int_eq(index, 0);

        /* Check the invalid material case */
        reset_error();
        pumas_material_index("Strawberry", &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_UNKNOWN_MATERIAL);

        reset_error();
        pumas_material_name(4, &name);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Check the materials mapping */
        for (i = 0; i < sizeof(names) / sizeof(*names); i++) {
                reset_error();
                pumas_material_index(names[i], &index);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(index, i);

                reset_error();
                pumas_material_name(i, &name);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_str_eq(name, names[i]);
        }

        pumas_finalise();
}
END_TEST


/* Test the composites API */
START_TEST (test_api_composite)
{
        int i, index, components;
        const char * name;
        const char * names[] = { "Water", "StandardRock", "WetRock" };
        double density, fractions[2], densities[2];

        /* Check the initialisation error */
        reset_error();
        pumas_composite_properties(0, NULL, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_composite_update(0, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        index = pumas_composite_length();
        ck_assert_int_eq(index, 0);

        /* Load the tau data */
        load_tau();

        /* Check the number of materials and composites */
        index = pumas_material_length();
        ck_assert_int_eq(index, 3);

        index = pumas_composite_length();
        ck_assert_int_eq(index, 1);

        /* Check the materials mapping */
        for (i = 0; i < sizeof(names) / sizeof(*names); i++) {
                reset_error();
                pumas_material_index(names[i], &index);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(index, i);

                reset_error();
                pumas_material_name(i, &name);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_str_eq(name, names[i]);
        }

        /* Check the properties getter */
        reset_error();
        pumas_composite_properties(0, NULL, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_composite_properties(2, NULL, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        double density0 = TEST_ROCK_DENSITY, density1 = 1.00E+03;
        double fraction0 = 0.5, fraction1 = 0.5;
        double rho = 1. / (fraction0 / density0 + fraction1 / density1);

        pumas_composite_properties(
            2, &density, &components, fractions, densities);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(components, 2);
        ck_assert_double_eq(density, rho);
        ck_assert_double_eq(fractions[0], fraction0);
        ck_assert_double_eq(fractions[1], fraction1);
        ck_assert_double_eq(densities[0], density0);
        ck_assert_double_eq(densities[1], density1);

        /* Check the properties setter */
        reset_error();
        pumas_composite_update(2, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        pumas_composite_properties(
            2, &density, &components, fractions, densities);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(components, 2);
        ck_assert_double_eq(density, rho);
        ck_assert_double_eq(fractions[0], fraction0);
        ck_assert_double_eq(fractions[1], fraction1);
        ck_assert_double_eq(densities[0], density0);
        ck_assert_double_eq(densities[1], density1);

        density0 = 2.3E+03, density = 1.2E+03;
        fraction0 = 0.3, fraction1 = 0.7;
        rho = 1. / (fraction0 / density0 + fraction1 / density1);

        {
                double tmp0[] = { fraction0, fraction1 };
                double tmp1[] = { density0, density1 };
                pumas_composite_update(2, tmp0, tmp1);
        }
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        reset_error();
        pumas_composite_properties(
            2, &density, &components, fractions, densities);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(components, 2);
        ck_assert_double_eq(density, rho);
        ck_assert_double_eq(fractions[0], fraction0);
        ck_assert_double_eq(fractions[1], fraction1);
        ck_assert_double_eq(densities[0], density0);
        ck_assert_double_eq(densities[1], density1);

        pumas_finalise();
}
END_TEST


START_TEST (test_api_property)
{
        double value;

        /* Check the initialisation error */
        reset_error();
        pumas_property_cross_section(0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_property_energy_loss(0, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_property_grammage(0, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_property_kinetic_energy(0, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_property_magnetic_rotation(0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_property_proper_time(0, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_property_scattering_length(0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        /* Load the muon data */
        load_muon();

        /* Check the index error */
        reset_error();
        pumas_property_cross_section(3, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_energy_loss(PUMAS_SCHEME_HYBRID, 3, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_energy_loss(PUMAS_SCHEME_DETAILED, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_grammage(PUMAS_SCHEME_HYBRID, 3, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_grammage(PUMAS_SCHEME_DETAILED, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_kinetic_energy(PUMAS_SCHEME_HYBRID, 3, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_kinetic_energy(PUMAS_SCHEME_DETAILED, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_magnetic_rotation(3, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 3, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_proper_time(PUMAS_SCHEME_DETAILED, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_property_scattering_length(3, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Check some values */
        reset_error();
        pumas_property_cross_section(0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_property_cross_section(0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_lt(1. / (value * TEST_ROCK_DENSITY), 1E+03);

        pumas_property_cross_section(0, 1E+06, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_lt(1. / (value * TEST_ROCK_DENSITY), 1E+03);

        pumas_property_energy_loss(PUMAS_SCHEME_CSDA, 0, 1E+00, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 1.808E-04);

        pumas_property_energy_loss(PUMAS_SCHEME_CSDA, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 6.604E-04);

        pumas_property_energy_loss(PUMAS_SCHEME_HYBRID, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 6.604E-04);

        pumas_property_grammage(PUMAS_SCHEME_CSDA, 0, 1E+00, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 5.534E+03, 1.);

        pumas_property_grammage(PUMAS_SCHEME_CSDA, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 2.455E+06, 1E+03);

        pumas_property_grammage(PUMAS_SCHEME_HYBRID, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 2.455E+06);

        pumas_property_kinetic_energy(PUMAS_SCHEME_CSDA, 0, 5.534E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 1E+00, 1E-03);

        pumas_property_kinetic_energy(PUMAS_SCHEME_CSDA, 0, 2.455E+06, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 1E+03, 1.);

        pumas_property_kinetic_energy(PUMAS_SCHEME_HYBRID, 0, 2.455E+06, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 1E+03);

        pumas_property_magnetic_rotation(0, 1E+00, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_gt(value, 5.534E+03 / 3.3374);

        pumas_property_proper_time(PUMAS_SCHEME_CSDA, 0, 1E+00, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 5.534E+03);

        pumas_property_proper_time(PUMAS_SCHEME_CSDA, 0, 1E-02, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 8.498);

        pumas_property_scattering_length(0, 1E-02, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        value *= 1E-03; /* Normalise to electrons in Al */
        ck_assert_double_gt(value, 0.3);
        ck_assert_double_lt(value, 3.);

        /* Check overflows */
        double vmax;
        pumas_property_energy_loss(PUMAS_SCHEME_CSDA, 0, 1E+08, &value);
        pumas_table_value(
            PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_CSDA, 0, 145, &vmax);
        ck_assert_double_gt(value, vmax);

        pumas_property_grammage(PUMAS_SCHEME_CSDA, 0, 1E+08, &value);
        pumas_table_value(
            PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0, 145, &vmax);
        ck_assert_double_gt(value, vmax);

        pumas_property_proper_time(PUMAS_SCHEME_CSDA, 0, 1E+08, &value);
        pumas_table_value(
            PUMAS_PROPERTY_PROPER_TIME, PUMAS_SCHEME_CSDA, 0, 145, &vmax);
        ck_assert_double_gt(value, vmax);

        pumas_property_cross_section(0, 1E+08, &value);
        pumas_table_value(
            PUMAS_PROPERTY_CROSS_SECTION, PUMAS_SCHEME_HYBRID, 0, 145, &vmax);
        ck_assert_double_eq(value, vmax);

        /* Unload the data */
        pumas_finalise();
}
END_TEST


START_TEST (test_api_table)
{
        int index;
        double value;

        /* Check the initialisation error */
        reset_error();
        pumas_table_index(0, 0, 0, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        pumas_table_value(0, 0, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        reset_error();
        index = pumas_table_length();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 0);

        /* Load the muon data */
        load_muon();

        reset_error();
        index = pumas_table_length();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 146);

        /* Check the index error */
        reset_error();
        pumas_table_index(1024, 0, 0, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_table_index(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_DETAILED, 0,
            0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_table_index(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 3,
            0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_table_index(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0,
            1E+18, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        pumas_table_index(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0,
            -1., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        pumas_table_value(1024, 0, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_table_value(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_DETAILED, 0,
            0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_table_value(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 3,
            0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_table_value(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0,
            200, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_table_value(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0,
            -1, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Check some values */
        reset_error();

        pumas_table_index(PUMAS_PROPERTY_KINETIC_ENERGY, PUMAS_SCHEME_CSDA, 0,
            0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 0);

        pumas_table_index(PUMAS_PROPERTY_KINETIC_ENERGY, PUMAS_SCHEME_CSDA, 0,
            1E+00, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_table_index(PUMAS_PROPERTY_KINETIC_ENERGY, PUMAS_SCHEME_CSDA, 0,
            1E+06, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 145);

        pumas_table_value(PUMAS_PROPERTY_KINETIC_ENERGY, PUMAS_SCHEME_CSDA, 0,
            0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_table_value(PUMAS_PROPERTY_KINETIC_ENERGY, PUMAS_SCHEME_CSDA, 0,
            49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 1E+00);

        pumas_table_index(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0,
            0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 0);

        pumas_table_index(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0,
            5.534E+03, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_table_value(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0,
            0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_table_value(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_CSDA, 0,
            49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 5.534E+03, 1.);

        pumas_table_value(PUMAS_PROPERTY_GRAMMAGE, PUMAS_SCHEME_HYBRID, 0,
            49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 5.534E+03);

        pumas_table_value(PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_CSDA, 0,
            0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_table_value(PUMAS_PROPERTY_CROSS_SECTION, PUMAS_SCHEME_HYBRID, 0,
            0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_table_value(PUMAS_PROPERTY_MAGNETIC_ROTATION, 0, 0, 49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_gt(value, 5.534E+03 / 3.3374);

        pumas_table_index(
            PUMAS_PROPERTY_MAGNETIC_ROTATION, 0, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_table_value(PUMAS_PROPERTY_PROPER_TIME, PUMAS_SCHEME_CSDA, 0,
            49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 5.534E+03);

        pumas_table_index(
            PUMAS_PROPERTY_PROPER_TIME, PUMAS_SCHEME_CSDA, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_table_value(PUMAS_PROPERTY_PROPER_TIME, PUMAS_SCHEME_CSDA, 0,
            17, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 8.498);

        pumas_table_index(
            PUMAS_PROPERTY_PROPER_TIME, PUMAS_SCHEME_CSDA, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 17);

        /* Check unsuported properties  */
        pumas_table_value(
            PUMAS_PROPERTY_SCATTERING_LENGTH, 0, 0, 17, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        pumas_table_index(PUMAS_PROPERTY_CROSS_SECTION, 0,
            0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        pumas_table_index(PUMAS_PROPERTY_ENERGY_LOSS, 0,
            0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        pumas_table_index(PUMAS_PROPERTY_SCATTERING_LENGTH, 0,
            0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Unload the data */
        pumas_finalise();
}
END_TEST


static void * fail_malloc(size_t size)
{
        return NULL;
}

static void * fail_realloc(void * ptr, size_t size)
{
        return NULL;
}


/* Test the context API */
START_TEST (test_api_context)
{
        /* Test the initialisation error */
        reset_error();
        context = (void *)0x1;
        pumas_context_create(&context, 0);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);
        ck_assert_ptr_null(context);

        /* Load the muon data */
        load_muon();

        /* Test a memory error */
        pumas_memory_allocator(&fail_malloc);
        reset_error();
        context = (void *)0x1;
        pumas_context_create(&context, 0);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MEMORY_ERROR);
        ck_assert_ptr_null(context);
        pumas_memory_allocator(NULL);

        /* Test the finalisation of a NULL context */
        reset_error();
        pumas_context_destroy(NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        context = NULL;
        pumas_context_destroy(&context);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(context);

        /* Test the allocation of user_data */
        const int n = 1024;
        int i, * data;

        pumas_context_create(&context, n * sizeof *data);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(context->user_data);

        data = context->user_data;
        for (i = 0; i < n; i++) data[i] = i;

        /* Test the initialisation of a muon context */
        ck_assert_ptr_null(context->medium);
        ck_assert_ptr_null(context->random);
        ck_assert_ptr_null(context->recorder);
        ck_assert_int_eq(context->longitudinal, 0);
        ck_assert_int_eq(context->forward, 1);
        ck_assert_int_eq(context->scheme, PUMAS_SCHEME_DETAILED);
        ck_assert_int_eq(context->decay, PUMAS_DECAY_WEIGHT);
        ck_assert_int_eq(context->event, PUMAS_EVENT_NONE);

        ck_assert_double_eq(context->kinetic_limit, 0.);
        ck_assert_double_eq(context->distance_max, 0.);
        ck_assert_double_eq(context->grammage_max, 0.);
        ck_assert_double_eq(context->time_max, 0.);

        /* Test the context destruction */
        pumas_context_destroy(&context);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(context);

        /* Load the tau data */
        pumas_finalise();
        load_tau();

        /* Test the initialisation of a tau context */
        pumas_context_create(&context, -1);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(context->user_data);

        ck_assert_ptr_null(context->medium);
        ck_assert_ptr_null(context->random);
        ck_assert_ptr_null(context->recorder);
        ck_assert_int_eq(context->longitudinal, 0);
        ck_assert_int_eq(context->forward, 1);
        ck_assert_int_eq(context->scheme, PUMAS_SCHEME_DETAILED);
        ck_assert_int_eq(context->decay, PUMAS_DECAY_PROCESS);
        ck_assert_int_eq(context->event, PUMAS_EVENT_NONE);

        ck_assert_double_eq(context->kinetic_limit, 0.);
        ck_assert_double_eq(context->distance_max, 0.);
        ck_assert_double_eq(context->grammage_max, 0.);
        ck_assert_double_eq(context->time_max, 0.);

        /* Free the data */
        pumas_context_destroy(&context);
        pumas_finalise();
}
END_TEST


/* Test the recorder API */
START_TEST (test_api_recorder)
{
        /* Check the memory error */
        pumas_memory_allocator(&fail_malloc);
        reset_error();
        recorder = (void *)0x1;
        pumas_recorder_create(&recorder, 0);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MEMORY_ERROR);
        ck_assert_ptr_null(recorder);
        pumas_memory_allocator(NULL);

        /* Check the null dellocator */
        reset_error();
        recorder = NULL;

        pumas_recorder_destroy(NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        pumas_recorder_destroy(&recorder);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        /* Check the recorder initialisation */
        pumas_recorder_create(&recorder, 0);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(recorder->first);
        ck_assert_int_eq(recorder->length, 0);
        ck_assert_int_eq(recorder->period, 1);
        ck_assert_ptr_null(recorder->record);
        ck_assert_ptr_null(recorder->user_data);

        /* Test the destruction of a valid recorder */
        pumas_recorder_destroy(&recorder);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(recorder);

        /* Test the allocation of user_data */
        const int n = 1024;
        int i, * data;

        pumas_recorder_create(&recorder, n * sizeof *data);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(recorder->user_data);

        data = recorder->user_data;
        for (i = 0; i < n; i++) data[i] = i;
        pumas_recorder_destroy(&recorder);

        /* Test the NULL clear */
        pumas_recorder_clear(NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
}
END_TEST


/* Test the print API */
START_TEST (test_api_print)
{
        /* Test the initialisation_error */
        reset_error();
        pumas_print(NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INITIALISATION_ERROR);

        /* Load the tau data */
        load_tau();

        /* Check the NULL stream case */
        reset_error();
        pumas_print(NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_IO_ERROR);

        /* Check printing to a file */
        FILE * stream = fopen("tau.json", "w+");
        pumas_print(stream, NULL, NULL);
        fputs("\n===\n", stream);
        pumas_print(stream, "    ", "\n");
        fclose(stream);

        /* Free the data */
        pumas_finalise();
}
END_TEST


/* Geometry for test cases  */
static struct {
        int uniform;
        double magnet[3];
        double last_position[3];
} geometry = {1, {0., 0., 0.}, {-1., 1., -1.}};

static double rock_locals(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals)
{
        /* Check the API consistency */
        ck_assert_mem_eq(
            state->position, geometry.last_position, sizeof state->position);

        locals->density = TEST_ROCK_DENSITY;
        memcpy(locals->magnet, geometry.magnet, sizeof locals->magnet);

        return 0.;
}

static double air_locals(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals)
{
        /* Check the API consistency */
        ck_assert_mem_eq(
            state->position, geometry.last_position, sizeof state->position);

        const double lambda = 1E+04;

        locals->density = TEST_AIR_DENSITY * exp(-state->position[2] / lambda);
        memcpy(locals->magnet, geometry.magnet, sizeof locals->magnet);

        return 1E-02 * lambda;
}

static double geometry_medium(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr)
{
        memcpy(geometry.last_position, state->position,
            sizeof geometry.last_position);

        static struct pumas_medium media[2] = {
                {0, &rock_locals}, {1, &air_locals}};

        if (geometry.uniform) {
                if (medium_ptr != NULL) *medium_ptr = media;
                return 0.;
        } else {
                if (medium_ptr != NULL) *medium_ptr = NULL;
                const double z = state->position[2];
                if (z < -0.5 * TEST_ROCK_DEPTH) {
                        return 0.;
                } else if (z < 0.5 * TEST_ROCK_DEPTH) {
                        if (medium_ptr != NULL) *medium_ptr = media;

                        const double dz1 = z + 0.5 * TEST_ROCK_DEPTH;
                        const double dz2 = 0.5 * TEST_ROCK_DEPTH - z;
                        const double dz = dz1 < dz2 ? dz1 : dz2;
                        return dz > 0. ? dz : 1E-03;
                } else if (z < TEST_MAX_ALTITUDE) {
                        if (medium_ptr != NULL) *medium_ptr = media + 1;

                        const double sgn = context->forward ? 1. : -1.;
                        const double uz = state->direction[2] * sgn;
                        if (fabs(uz) < 1E-03)
                                return 1E+03;

                        double s;
                        if (uz > 0.)
                                s = (TEST_MAX_ALTITUDE - z) / uz;
                        else
                                s = (0.5 * TEST_ROCK_DEPTH - z) / uz;
                        return s > 0. ? s : 1E-03;
                } else {
                        return 0.;
                }
        }
}


/* Fixtures for Lossless tests */
static void lossless_setup(void)
{
        /* Load the muon data and create a simulation context */
        load_muon();
        geometry.uniform = 1;
        pumas_context_create(&context, 0);
        context->medium = &geometry_medium;
        context->longitudinal = 1;
        context->scheme = PUMAS_SCHEME_NO_LOSS;
}

static void lossless_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_finalise();
}


/* Helper function for initialising a particle state */
static void initialise_state(void)
{
        state->charge = -1.;
        state->distance = 0.;
        state->grammage = 0.;
        state->time = 0.;
        state->weight = 1.;
        memset(state->position, 0x0, sizeof state->position);
        memset(state->direction, 0x0, sizeof state->direction);
        state->direction[2] = 1.;
        state->decayed = 0;
}


/* Test the Lossless straight case */
START_TEST (test_lossless_straight)
{
        int i, forward;
        enum pumas_event event_data, * event;
        struct pumas_medium * media_data[2], ** media;
        struct pumas_medium * rock;
        geometry_medium(context, state, &rock);

        double ctau, mu;
        pumas_particle(NULL, &ctau, &mu);

        /* Check some basic errors */
        reset_error();
        initialise_state();
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MISSING_LIMIT);


        for (forward = 0; forward < 2; forward++) {
                context->forward = forward;
                const double sgn = forward ? 1. : -1.;

                /* Check the forward distance limit */
                context->event = PUMAS_EVENT_LIMIT_DISTANCE;
                context->distance_max = 1E+03;

                double k = 1.;
                double d = context->distance_max;
                double X = d * TEST_ROCK_DENSITY;
                double gamma = k / mu + 1.;
                double t = d / sqrt(gamma * gamma - 1.);
                for (i = 0; i < 2; i++) {
                        if (i == 0)
                                event = NULL, media = NULL;
                        else
                                event = &event_data, media = media_data;

                        reset_error();
                        initialise_state();
                        state->kinetic = k;

                        pumas_transport(context, state, event, media);
                        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                        ck_assert_double_eq(state->charge, -1.);
                        ck_assert_double_eq(state->kinetic, k);
                        ck_assert_double_eq(state->distance, d);
                        ck_assert_double_eq(state->grammage, X);
                        ck_assert_double_eq_tol(state->time, t, FLT_EPSILON);
                        ck_assert_double_eq_tol(
                            state->weight, exp(-state->time / ctau),
                            FLT_EPSILON);
                        ck_assert_double_eq(state->position[0], 0.);
                        ck_assert_double_eq(state->position[1], 0.);
                        ck_assert_double_eq_tol(
                            state->position[2], sgn * state->distance,
                            FLT_EPSILON);
                        ck_assert_double_eq(state->direction[0], 0.);
                        ck_assert_double_eq(state->direction[1], 0.);
                        ck_assert_double_eq(state->direction[2], 1.);
                        ck_assert_int_eq(state->decayed, 0);

                        if (i == 0) {
                                ck_assert_ptr_null(event);
                                ck_assert_ptr_null(media);
                        } else {
                                ck_assert_int_eq(
                                    *event, PUMAS_EVENT_LIMIT_DISTANCE);
                                ck_assert_ptr_eq(media[0], rock);
                                ck_assert_ptr_eq(media[0], media[1]);
                        }
                }
                context->distance_max = 0.;

                /* Check the forward grammage limit */
                context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;
                context->grammage_max = 1E+03;

                X = context->grammage_max;
                d = X / TEST_ROCK_DENSITY;
                t = d / sqrt(gamma * gamma - 1.);
                for (i = 0; i < 2; i++) {
                        if (i == 0)
                                event = NULL, media = NULL;
                        else
                                event = &event_data, media = media_data;

                        reset_error();
                        initialise_state();
                        state->kinetic = k;

                        pumas_transport(context, state, event, media);
                        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                        ck_assert_double_eq(state->charge, -1.);
                        ck_assert_double_eq(state->kinetic, k);
                        ck_assert_double_eq(state->distance, d);
                        ck_assert_double_eq(state->grammage, X);
                        ck_assert_double_eq_tol(state->time, t, FLT_EPSILON);
                        ck_assert_double_eq_tol(
                            state->weight, exp(-state->time / ctau),
                            FLT_EPSILON);
                        ck_assert_double_eq(state->position[0], 0.);
                        ck_assert_double_eq(state->position[1], 0.);
                        ck_assert_double_eq_tol(
                            state->position[2], sgn * state->distance,
                            FLT_EPSILON);
                        ck_assert_double_eq(state->direction[0], 0.);
                        ck_assert_double_eq(state->direction[1], 0.);
                        ck_assert_double_eq(state->direction[2], 1.);
                        ck_assert_int_eq(state->decayed, 0);

                        if (i == 0) {
                                ck_assert_ptr_null(event);
                                ck_assert_ptr_null(media);
                        } else {
                                ck_assert_int_eq(*event,
                                    PUMAS_EVENT_LIMIT_GRAMMAGE);
                                ck_assert_ptr_eq(media[0], rock);
                                ck_assert_ptr_eq(media[0], media[1]);
                        }
                }

                /* Check the forward time limit */
                context->event = PUMAS_EVENT_LIMIT_TIME;
                context->time_max = 1E+03;

                t = context->time_max;
                d = t * sqrt(gamma * gamma - 1.);
                X = d * TEST_ROCK_DENSITY;
                for (i = 0; i < 2; i++) {
                        if (i == 0)
                                event = NULL, media = NULL;
                        else
                                event = &event_data, media = media_data;

                        reset_error();
                        initialise_state();
                        state->kinetic = k;

                        pumas_transport(context, state, event, media);
                        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                        ck_assert_double_eq(state->charge, -1.);
                        ck_assert_double_eq(state->kinetic, k);
                        ck_assert_double_eq(state->distance, d);
                        ck_assert_double_eq(state->grammage, X);
                        ck_assert_double_eq_tol(state->time, t, FLT_EPSILON);
                        ck_assert_double_eq_tol(
                            state->weight, exp(-state->time / ctau),
                            FLT_EPSILON);
                        ck_assert_double_eq(state->position[0], 0.);
                        ck_assert_double_eq(state->position[1], 0.);
                        ck_assert_double_eq_tol(
                            state->position[2], sgn * state->distance,
                            FLT_EPSILON);
                        ck_assert_double_eq(state->direction[0], 0.);
                        ck_assert_double_eq(state->direction[1], 0.);
                        ck_assert_double_eq(state->direction[2], 1.);
                        ck_assert_int_eq(state->decayed, 0);

                        if (i == 0) {
                                ck_assert_ptr_null(event);
                                ck_assert_ptr_null(media);
                        } else {
                                ck_assert_int_eq(
                                    *event, PUMAS_EVENT_LIMIT_TIME);
                                ck_assert_ptr_eq(media[0], rock);
                                ck_assert_ptr_eq(media[0], media[1]);
                        }
                }
        }
}
END_TEST


START_TEST (test_lossless_magnet)
{
        geometry.magnet[1] = 0.1;
        context->distance_max = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        double ctau, mu;
        pumas_particle(NULL, &ctau, &mu);

        /* Check the MC deflection and the analytical computation when no
         * energy loss
         */
        const double tol = 1E-02;
        const double k = 1E+00;
        const double gamma = k / mu + 1.;
        const double bg = sqrt(gamma * gamma - 1.);
        const double d = context->distance_max;
        const double rL = mu * bg / (geometry.magnet[1] * 0.299792458);
        const double X = d * TEST_ROCK_DENSITY;
        const double t = d / bg;
        const double phi = d / rL;
        double * r = state->position, * u = state->direction;

        reset_error();
        initialise_state();
        state->kinetic = k;

        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->kinetic, k);
        ck_assert_double_eq(state->distance, d);
        ck_assert_double_eq_tol(state->grammage, X, X * 1E-09);
        ck_assert_double_eq_tol(state->time, t, FLT_EPSILON);
        ck_assert_double_eq_tol(
            state->weight, exp(-state->time / ctau), FLT_EPSILON);

        ck_assert_double_eq_tol(r[0], rL * (1. - cos(phi)), rL * tol);
        ck_assert_double_eq(r[1], 0.);
        ck_assert_double_eq_tol(r[2], rL * sin(phi), rL * tol);

        ck_assert_double_eq_tol(u[0], sin(phi), tol);
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_eq_tol(u[2], cos(phi), tol);

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Erase the magnetic field and restore the context */
        geometry.magnet[1] = 0.;
        context->distance_max = 0.;
        context->event = PUMAS_EVENT_NONE;
}
END_TEST


START_TEST (test_lossless_geometry)
{
        int i;
        enum pumas_event event;
        struct pumas_medium * media[2], * rock, * air;
        struct pumas_frame * frame;
        pumas_recorder_create(&recorder, 0);
        context->recorder = recorder;
        recorder->period = 0;

        geometry.uniform = 0;
        initialise_state();
        geometry_medium(context, state, &rock);
        state->position[2] = 0.5 * (0.5 * TEST_ROCK_DEPTH + TEST_MAX_ALTITUDE);
        geometry_medium(context, state, &air);

        /* Test the initially out of world case */
        state->position[2] = 2 * TEST_MAX_ALTITUDE;
        reset_error();
        pumas_transport(context, state, &event, media);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
        ck_assert_ptr_eq(media[0], NULL);
        ck_assert_ptr_eq(media[1], NULL);
        ck_assert_int_eq(recorder->length, 1);
        frame = recorder->first;
        ck_assert_ptr_eq(frame->medium, NULL);
        ck_assert_int_eq(frame->event, PUMAS_EVENT_START | PUMAS_EVENT_MEDIUM |
            PUMAS_EVENT_STOP);
        ck_assert_ptr_null(frame->next);

        for (i = 0; i < 2; i++) {
                /* Test the unconstrained transport */
                context->event = PUMAS_EVENT_NONE;
                initialise_state();
                state->kinetic = 1E+03;
                if (i) {
                        context->forward = 0;
                        state->direction[2] = -1.;
                } else {
                        context->forward = 1;
                }

                reset_error();
                pumas_recorder_clear(recorder);
                pumas_transport(context, state, &event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
                ck_assert_ptr_eq(media[0], rock);
                ck_assert_ptr_eq(media[1], NULL);
                ck_assert_double_eq_tol(
                    state->position[2], TEST_MAX_ALTITUDE, FLT_EPSILON);
                ck_assert_int_eq(recorder->length, 3);
                frame = recorder->first;
                ck_assert_ptr_eq(frame->medium, rock);
                ck_assert_int_eq(frame->event, PUMAS_EVENT_START);
                frame = frame->next;
                ck_assert_ptr_eq(frame->medium, air);
                ck_assert_int_eq(frame->event, PUMAS_EVENT_MEDIUM);
                frame = frame->next;
                ck_assert_ptr_null(frame->medium);
                ck_assert_int_eq(frame->event,
                    PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);

                /* Test the medium stop condition */
                context->event = PUMAS_EVENT_MEDIUM;
                initialise_state();
                state->kinetic = 1E+03;
                if (i) {
                        state->direction[2] = -1.;
                }

                reset_error();
                pumas_recorder_clear(recorder);
                pumas_transport(context, state, &event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
                ck_assert_ptr_eq(media[0], rock);
                ck_assert_ptr_eq(media[1], air);
                ck_assert_int_eq(recorder->length, 2);
                frame = recorder->first;
                ck_assert_ptr_eq(frame->medium, rock);
                ck_assert_int_eq(frame->event, PUMAS_EVENT_START);
                frame = frame->next;
                ck_assert_ptr_eq(frame->medium, air);
                ck_assert_int_eq(frame->event,
                    PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);

                pumas_recorder_clear(recorder);
                pumas_transport(context, state, &event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
                ck_assert_ptr_eq(media[0], air);
                ck_assert_ptr_eq(media[1], NULL);
                ck_assert_int_eq(recorder->length, 2);
                frame = recorder->first;
                ck_assert_ptr_eq(frame->medium, air);
                ck_assert_int_eq(frame->event, PUMAS_EVENT_START);
                frame = frame->next;
                ck_assert_ptr_null(frame->medium);
                ck_assert_int_eq(frame->event,
                    PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);
        }

        pumas_recorder_destroy(&recorder);
        context->recorder = NULL;
        geometry.uniform = 1;
        context->event = PUMAS_EVENT_NONE;
        context->forward = 1;
}
END_TEST


/* Fixtures for CSDA tests */
static void csda_setup(void)
{
        /* Load the muon data and create a simulation context */
        load_muon();
        geometry.uniform = 1;
        pumas_context_create(&context, 0);
        context->medium = &geometry_medium;
        context->longitudinal = 1;
        context->scheme = PUMAS_SCHEME_CSDA;
}

static void csda_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_finalise();
}


/* Test the CSDA straight case */
START_TEST (test_csda_straight)
{
        int i;
        enum pumas_event event_data, * event;
        struct pumas_medium * media_data[2], ** media;
        struct pumas_medium * rock;
        geometry_medium(context, state, &rock);

        double ctau;
        pumas_particle(NULL, &ctau, NULL);

        /* Check some basic errors */
        reset_error();
        pumas_transport(NULL, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        pumas_transport(context, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        initialise_state();
        state->decayed = 1;
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        event = &event_data;
        pumas_transport(context, state, event, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(*event, PUMAS_EVENT_VERTEX_DECAY);

        context->medium = NULL;
        reset_error();
        initialise_state();
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MEDIUM_ERROR);
        context->medium = &geometry_medium;

        context->random = NULL;
        context->decay = PUMAS_DECAY_PROCESS;
        reset_error();
        initialise_state();
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MISSING_RANDOM);
        context->random = &uniform01;
        context->decay = PUMAS_DECAY_WEIGHT;

        /* Check the forward full path transport */
        double X, d, t0;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_property_grammage(
                    PUMAS_SCHEME_CSDA, 0, state->kinetic, &X);
                pumas_property_proper_time(
                    PUMAS_SCHEME_CSDA, 0, state->kinetic, &t0);
                t0 /= TEST_ROCK_DENSITY;
                d = X / TEST_ROCK_DENSITY;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->kinetic, 0.);
                ck_assert_double_eq(state->distance, d);
                ck_assert_double_eq(state->grammage, X);
                ck_assert_double_eq(state->time, t0);
                ck_assert_double_eq_tol(
                    state->weight, exp(-t0 / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], d);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        /* Check the no event case but with limit values set */
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                context->kinetic_limit = 0.5;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->kinetic_limit = 0.;
                context->distance_max = 0.5 * d;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->distance_max = 0.;
                context->grammage_max = 0.5 * X;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->grammage_max = 0.;
                context->time_max = 0.5 * t0;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                context->time_max = 0.;
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        /* Check the forward kinetic limit */
        double X1, d1, t1;
        context->kinetic_limit = 0.5;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;
        pumas_property_grammage(
            PUMAS_SCHEME_CSDA, 0, context->kinetic_limit, &X1);
        pumas_property_proper_time(
            PUMAS_SCHEME_CSDA, 0, context->kinetic_limit, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, context->kinetic_limit);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-(t0 - t1) / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], d - d1);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->kinetic_limit = 0.;

        /* Check the forward grammage limit */
        double k1;
        X1 = 0.5 * X;
        d1 = X1 / TEST_ROCK_DENSITY;
        context->grammage_max = X1;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_proper_time(PUMAS_SCHEME_CSDA, 0, k1, &t1);
                t1 /= TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-(t0 - t1) / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], d - d1);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the forward distance limit */
        d1 = 0.5 * d;
        X1 = d1 * TEST_ROCK_DENSITY;
        context->distance_max = d1;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_proper_time(PUMAS_SCHEME_CSDA, 0, k1, &t1);
                t1 /= TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-(t0 - t1) / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], d - d1);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->distance_max = 0.;

        /* Check the forward time limit */
        t1 = 0.5 * t0;
        context->time_max = t1;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_grammage(PUMAS_SCHEME_CSDA, 0, k1, &X1);
                d1 = X1 / TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-(t0 - t1) / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], d - d1);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->time_max = 0.;

        /* Check the backward transport with a kinetic limit */
        context->forward = 0;
        context->kinetic_limit = 1.;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;
        k1 = 0.5;
        pumas_property_grammage(PUMAS_SCHEME_CSDA, 0, k1, &X1);
        pumas_property_proper_time(PUMAS_SCHEME_CSDA, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        double de0, de1;
        pumas_property_energy_loss(
            PUMAS_SCHEME_CSDA, 0, context->kinetic_limit, &de0);
        pumas_property_energy_loss(PUMAS_SCHEME_CSDA, 0, k1, &de1);
        double w = exp(-(t0 - t1) / ctau) * de0 / de1;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->kinetic, context->kinetic_limit);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->weight, w, FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], -(d - d1));
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->kinetic_limit = 0.;

        /* Check the backward transport with a grammage limit */
        context->grammage_max = X - X1;
        context->kinetic_limit = 0.75; /* Not activated  */
        context->distance_max = 0.5 * (d - d1); /* Not activated  */
        context->time_max = 0.5 * (t0 - t1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->kinetic, 1.);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->weight, w, FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], -(d - d1));
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the backward transport with a distance limit */
        context->distance_max = d - d1;
        context->grammage_max = 0.5 * (X - X1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->kinetic, 1.);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->weight, w, FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], -(d - d1));
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the backward transport with a time limit */
        context->time_max = t0 - t1;
        context->distance_max = 0.5 * (d - d1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq_tol(state->kinetic, 1., FLT_EPSILON);
                ck_assert_double_eq_tol(state->distance, d - d1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->weight, w, FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq_tol(
                    state->position[2], -(d - d1), FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;
        context->time_max = 0.5 * (t0 - t1); /* Not activated */

        /* Check the case of an initial kinetic limit violation */
        context->kinetic_limit = 0.5;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1.);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        context->kinetic_limit = 1;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 0.5;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.5);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->kinetic_limit = 0.;

        /* Check the case of an initial grammage limit violation */
        context->grammage_max = 0.5 * X;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                state->grammage = X;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1.);
                ck_assert_double_eq(state->grammage, X);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                state->grammage = X;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1.);
                ck_assert_double_eq(state->grammage, X);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->grammage_max = 0.;

        /* Check the case of an initial distance limit violation */
        context->distance_max = 0.5 * d;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                state->distance = d;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1.);
                ck_assert_double_eq(state->distance, d);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                state->distance = d;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1.);
                ck_assert_double_eq(state->distance, d);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->distance_max = 0.;

        /* Check the case of an initial time limit violation */
        context->time_max = 0.5 * t0;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                state->time = t0;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1.);
                ck_assert_double_eq(state->time, t0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1.;
                state->time = t0;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1.);
                ck_assert_double_eq(state->time, t0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->time_max = 0.;
        context->event = PUMAS_EVENT_NONE;

        /* Check the overflow case */
        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+08;
                pumas_property_grammage(
                    PUMAS_SCHEME_CSDA, 0, state->kinetic, &X);
                pumas_property_proper_time(
                    PUMAS_SCHEME_CSDA, 0, state->kinetic, &t0);
                t0 /= TEST_ROCK_DENSITY;
                d = X / TEST_ROCK_DENSITY;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->kinetic, 0.);
                ck_assert_double_eq(state->distance, d);
                ck_assert_double_eq(state->grammage, X);
                ck_assert_double_eq(state->time, t0);
                ck_assert_double_eq_tol(
                    state->weight, exp(-t0 / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq(state->position[2], d);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
}
END_TEST


START_TEST (test_csda_record)
{
        double ctau;
        pumas_particle(NULL, &ctau, NULL);

        struct pumas_medium * rock;
        geometry_medium(context, state, &rock);

        pumas_recorder_create(&recorder, 0);
        context->recorder = recorder;

        reset_error();
        initialise_state();
        state->kinetic = 1.;

        double X, d, t0;
        pumas_property_grammage(
            PUMAS_SCHEME_CSDA, 0, state->kinetic, &X);
        pumas_property_proper_time(
            PUMAS_SCHEME_CSDA, 0, state->kinetic, &t0);
        t0 /= TEST_ROCK_DENSITY;
        d = X / TEST_ROCK_DENSITY;
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_int_eq(recorder->length, 2);
        ck_assert_ptr_nonnull(recorder->first);

        struct pumas_frame * frame = recorder->first;
        ck_assert_ptr_eq(frame->medium, rock);
        ck_assert_double_eq(frame->state.charge, -1.);
        ck_assert_double_eq(frame->state.kinetic, 1.);
        ck_assert_double_eq(frame->state.distance, 0.);
        ck_assert_double_eq(frame->state.grammage, 0.);
        ck_assert_double_eq(frame->state.time, 0.);
        ck_assert_double_eq(frame->state.weight, 1.);
        ck_assert_double_eq(frame->state.position[0], 0.);
        ck_assert_double_eq(frame->state.position[1], 0.);
        ck_assert_double_eq(frame->state.position[2], 0.);
        ck_assert_double_eq(frame->state.direction[0], 0.);
        ck_assert_double_eq(frame->state.direction[1], 0.);
        ck_assert_double_eq(frame->state.direction[2], 1.);
        ck_assert_int_eq(frame->state.decayed, 0);
        ck_assert_int_eq(frame->event, PUMAS_EVENT_START);
        ck_assert_ptr_nonnull(frame->next);

        frame = frame->next;
        ck_assert_ptr_eq(frame->medium, rock);
        ck_assert_double_eq(frame->state.charge, -1.);
        ck_assert_double_eq(frame->state.kinetic, 0.);
        ck_assert_double_eq(frame->state.distance, d);
        ck_assert_double_eq(frame->state.grammage, X);
        ck_assert_double_eq(frame->state.time, t0);
        ck_assert_double_eq_tol(
            frame->state.weight, exp(-t0 / ctau), FLT_EPSILON);
        ck_assert_double_eq(frame->state.position[0], 0.);
        ck_assert_double_eq(frame->state.position[1], 0.);
        ck_assert_double_eq(frame->state.position[2], d);
        ck_assert_double_eq(frame->state.direction[0], 0.);
        ck_assert_double_eq(frame->state.direction[1], 0.);
        ck_assert_double_eq(frame->state.direction[2], 1.);
        ck_assert_int_eq(frame->state.decayed, 0);
        ck_assert_int_eq(frame->event,
            PUMAS_EVENT_LIMIT_KINETIC | PUMAS_EVENT_STOP);
        ck_assert_ptr_null(frame->next);

        pumas_recorder_destroy(&recorder);
        context->recorder = NULL;
}
END_TEST


START_TEST (test_csda_magnet)
{
        geometry.magnet[1] = 0.1;

        reset_error();
        initialise_state();
        state->kinetic = 1.;

        double X, d, phi0, phi1, phi, * u = state->direction;
        pumas_property_grammage(
            PUMAS_SCHEME_CSDA, 0, state->kinetic, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_property_magnetic_rotation(0, state->kinetic, &phi0);
        pumas_property_magnetic_rotation(0, 0., &phi1);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_eq(u[0], sin(phi));
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_eq(u[2], cos(phi));

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        geometry.magnet[1] = 0.;
}
END_TEST


START_TEST (test_csda_geometry)
{
        int i;
        enum pumas_event event;
        struct pumas_medium * media[2], * rock, * air;
        struct pumas_frame * frame;
        pumas_recorder_create(&recorder, 0);
        context->recorder = recorder;
        recorder->period = 0;

        geometry.uniform = 0;
        initialise_state();
        geometry_medium(context, state, &rock);
        state->position[2] = 0.5 * (0.5 * TEST_ROCK_DEPTH + TEST_MAX_ALTITUDE);
        geometry_medium(context, state, &air);

        /* Test the initially out of world case */
        state->position[2] = 2 * TEST_MAX_ALTITUDE;
        reset_error();
        pumas_transport(context, state, &event, media);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
        ck_assert_ptr_eq(media[0], NULL);
        ck_assert_ptr_eq(media[1], NULL);
        ck_assert_int_eq(recorder->length, 1);
        frame = recorder->first;
        ck_assert_ptr_eq(frame->medium, NULL);
        ck_assert_int_eq(frame->event, PUMAS_EVENT_START | PUMAS_EVENT_MEDIUM |
            PUMAS_EVENT_STOP);
        ck_assert_ptr_null(frame->next);

        for (i = 0; i < 2; i++) {
                /* Test the unconstrained transport */
                context->event = PUMAS_EVENT_NONE;
                initialise_state();
                state->kinetic = 1E+03;
                if (i) {
                        context->forward = 0;
                        state->direction[2] = -1.;
                } else {
                        context->forward = 1;
                }

                reset_error();
                pumas_recorder_clear(recorder);
                pumas_transport(context, state, &event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
                ck_assert_ptr_eq(media[0], rock);
                ck_assert_ptr_eq(media[1], NULL);
                ck_assert_double_eq_tol(
                    state->position[2], TEST_MAX_ALTITUDE, FLT_EPSILON);
                ck_assert_int_eq(recorder->length, 3);
                frame = recorder->first;
                ck_assert_ptr_eq(frame->medium, rock);
                ck_assert_int_eq(frame->event, PUMAS_EVENT_START);
                frame = frame->next;
                ck_assert_ptr_eq(frame->medium, air);
                ck_assert_int_eq(frame->event, PUMAS_EVENT_MEDIUM);
                frame = frame->next;
                ck_assert_ptr_null(frame->medium);
                ck_assert_int_eq(frame->event,
                    PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);

                /* Test the medium stop condition */
                context->event = PUMAS_EVENT_MEDIUM;
                initialise_state();
                state->kinetic = 1E+03;
                if (i) {
                        state->direction[2] = -1.;
                }

                reset_error();
                pumas_recorder_clear(recorder);
                pumas_transport(context, state, &event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
                ck_assert_ptr_eq(media[0], rock);
                ck_assert_ptr_eq(media[1], air);
                ck_assert_int_eq(recorder->length, 2);
                frame = recorder->first;
                ck_assert_ptr_eq(frame->medium, rock);
                ck_assert_int_eq(frame->event, PUMAS_EVENT_START);
                frame = frame->next;
                ck_assert_ptr_eq(frame->medium, air);
                ck_assert_int_eq(frame->event,
                    PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);

                pumas_recorder_clear(recorder);
                pumas_transport(context, state, &event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
                ck_assert_ptr_eq(media[0], air);
                ck_assert_ptr_eq(media[1], NULL);
                ck_assert_int_eq(recorder->length, 2);
                frame = recorder->first;
                ck_assert_ptr_eq(frame->medium, air);
                ck_assert_int_eq(frame->event, PUMAS_EVENT_START);
                frame = frame->next;
                ck_assert_ptr_null(frame->medium);
                ck_assert_int_eq(frame->event,
                    PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);
        }

        pumas_recorder_destroy(&recorder);
        context->recorder = NULL;
        geometry.uniform = 1;
        context->event = PUMAS_EVENT_NONE;
        context->forward = 1;
}
END_TEST


/* Fixtures for hybrid tests */
static void hybrid_setup(void)
{
        /* Load the tau data and create a simulation context */
        load_muon();
        geometry.uniform = 1;
        pumas_context_create(&context, 0);
        context->medium = &geometry_medium;
        context->longitudinal = 1;
        context->scheme = PUMAS_SCHEME_HYBRID;
        context->random = &uniform01;
}

static void hybrid_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_finalise();
}


START_TEST (test_hybrid_straight)
{
        int i, ik;
        enum pumas_event event_data, * event;
        struct pumas_medium * media_data[2], ** media;
        struct pumas_medium * rock;
        geometry_medium(context, state, &rock);

        double ctau;
        pumas_particle(NULL, &ctau, NULL);

        /* Check the missing random engine case */
        reset_error();
        context->random = NULL;
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MISSING_RANDOM);
        context->random = &uniform01;

        /* Check the forward full path transport */
        double X, d, t0;
        for (ik = 0; ik < 2; ik++) {
                for (i = 0; i < 2; i++) {
                        if (i == 0)
                                event = NULL, media = NULL;
                        else
                                event = &event_data, media = media_data;

                        reset_error();
                        initialise_state();
                        state->kinetic = (ik == 0) ? 1E+08 : 1E+03;
                        pumas_property_grammage(
                            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
                        pumas_property_proper_time(
                            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &t0);
                        t0 /= TEST_ROCK_DENSITY;
                        d = X / TEST_ROCK_DENSITY;
                        pumas_transport(context, state, event, media);
                        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                        ck_assert_double_eq(state->charge, -1.);
                        ck_assert_double_eq(state->kinetic, 0.);
                        ck_assert_double_le(state->distance, d);
                        ck_assert_double_le(state->grammage, X);
                        ck_assert_double_le(state->time, t0);
                        ck_assert_double_eq_tol(
                            state->weight, exp(-state->time / ctau),
                            FLT_EPSILON);
                        ck_assert_double_eq(state->position[0], 0.);
                        ck_assert_double_eq(state->position[1], 0.);
                        ck_assert_double_le(
                            state->position[2], d + FLT_EPSILON);
                        ck_assert_double_eq(state->direction[0], 0.);
                        ck_assert_double_eq(state->direction[1], 0.);
                        ck_assert_double_eq(state->direction[2], 1.);
                        ck_assert_int_eq(state->decayed, 0);

                        if (i == 0) {
                                ck_assert_ptr_null(event);
                                ck_assert_ptr_null(media);
                        } else {
                                ck_assert_int_eq(
                                    *event, PUMAS_EVENT_LIMIT_KINETIC);
                                ck_assert_ptr_eq(media[0], rock);
                                ck_assert_ptr_eq(media[0], media[1]);
                        }
                }
        }

        /* Check the no event case but with limit values set */
        context->event = PUMAS_EVENT_NONE;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                context->kinetic_limit = 0.5E+03;
                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->kinetic_limit = 0.;
                context->distance_max = 0.5 * d;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->distance_max = 0.;
                context->grammage_max = 0.5 * X;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->grammage_max = 0.;
                context->time_max = 0.5 * t0;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                context->time_max = 0.;
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        /* Check the forward kinetic limit */
        double X1, d1, t1;
        context->kinetic_limit = 0.5E+03;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;
        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, context->kinetic_limit, &X1);
        pumas_property_proper_time(
            PUMAS_SCHEME_HYBRID, 0, context->kinetic_limit, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->kinetic, context->kinetic_limit);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], d - d1 + FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->kinetic_limit = 0.;

        /* Check the forward grammage limit */
        double k1;
        X1 = 1E-02 * X;
        d1 = X1 / TEST_ROCK_DENSITY;
        context->grammage_max = X1;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 0, k1, &t1);
                t1 /= TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], d - d1 + FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the forward distance limit */
        d1 = 1E-02 * d;
        X1 = d1 * TEST_ROCK_DENSITY;
        context->distance_max = d1;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 0, k1, &t1);
                t1 /= TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], d - d1 + FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->distance_max = 0.;

        /* Check the forward time limit */
        t1 = 1E-02 * t0;
        context->time_max = t1;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_grammage(PUMAS_SCHEME_HYBRID, 0, k1, &X1);
                d1 = X1 / TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], d - d1 + FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->time_max = 0.;

        /* Check the backward transport with a kinetic limit */
        context->forward = 0;
        context->kinetic_limit = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;
        k1 = 0.5 * context->kinetic_limit;
        pumas_property_grammage(PUMAS_SCHEME_HYBRID, 0, k1, &X1);
        pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_ge(state->kinetic, context->kinetic_limit);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_ge(
                    state->position[2], -(d - d1) - FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->kinetic_limit = 0.;

        /* Check the backward transport with a grammage limit */
        context->grammage_max = X - X1;
        context->kinetic_limit = 0.75E+03; /* Not activated  */
        context->distance_max = 0.5 * (d - d1); /* Not activated  */
        context->time_max = 0.5 * (t0 - t1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->kinetic, k1);
                ck_assert_double_eq_tol(state->distance, d - d1, FLT_EPSILON);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq_tol(
                    state->position[2], -(d - d1), FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the backward transport with a distance limit */
        context->distance_max = d - d1;
        context->grammage_max = 0.5 * (X - X1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->kinetic, k1);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq_tol(
                    state->position[2], -(d - d1), FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the backward transport with a time limit */
        context->time_max = t0 - t1;
        context->distance_max = 0.5 * (d - d1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->kinetic, k1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_lt(state->position[2], 0.);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;
        context->time_max = 0.5 * (t0 - t1); /* Not activated */

        /* Check the case of an initial kinetic limit violation */
        context->kinetic_limit = 0.5E+03;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        context->kinetic_limit = 1E+03;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 0.5E+03;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.5E+03);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->kinetic_limit = 0.;

        /* Check the case of an initial grammage limit violation */
        context->grammage_max = 0.5 * X;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->grammage = X;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->grammage, X);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->grammage = X;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->grammage, X);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->grammage_max = 0.;

        /* Check the case of an initial distance limit violation */
        context->distance_max = 0.5 * d;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->distance = d;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->distance, d);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->distance = d;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->distance, d);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->distance_max = 0.;

        /* Check the case of an initial time limit violation */
        context->time_max = 0.5 * t0;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->time = t0;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->time, t0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E03;
                state->time = t0;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->time, t0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        context->time_max = 0.;
        context->event = PUMAS_EVENT_NONE;
}
END_TEST


START_TEST (test_hybrid_scattering)
{
        context->longitudinal = 0;

        double ctau;
        pumas_particle(NULL, &ctau, NULL);

        /* Check the forward transport with scattering */
        double X, d, t0, * u = state->direction;

        reset_error();
        initialise_state();
        state->kinetic = 1E+03;

        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
        pumas_property_proper_time(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &t0);
        t0 /= TEST_ROCK_DENSITY;
        d = X / TEST_ROCK_DENSITY;

        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->charge, -1.);
        ck_assert_double_eq(state->kinetic, 0.);
        ck_assert_double_le(state->distance, d);
        ck_assert_double_le(state->grammage, X);
        ck_assert_double_le(state->time, t0);
        ck_assert_double_eq_tol(
            state->weight, exp(-state->time / ctau), FLT_EPSILON);
        ck_assert_double_le(state->position[2], d + FLT_EPSILON);
        ck_assert_int_eq(state->decayed, 0);

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Check the backward transport with scattering */
        context->forward = 0;
        context->kinetic_limit = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;

        double k1, X1, d1, t1;
        k1 = 0.5 * context->kinetic_limit;
        pumas_property_grammage(PUMAS_SCHEME_HYBRID, 0, k1, &X1);
        pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        reset_error();
        initialise_state();
        state->kinetic = k1;
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->charge, -1.);
        ck_assert_double_ge(state->kinetic, context->kinetic_limit);
        ck_assert_double_le(state->distance, d - d1 + FLT_EPSILON);
        ck_assert_double_le(state->grammage, X - X1 + FLT_EPSILON);
        ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
        ck_assert_double_ge(
            state->position[2], -(d - d1) - FLT_EPSILON);
        ck_assert_int_eq(state->decayed, 0);

        norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Restore the context status */
        context->event = PUMAS_EVENT_NONE;
        context->kinetic_limit = 0.;
        context->forward = 1;
        context->longitudinal = 1;
}
END_TEST


START_TEST (test_hybrid_record)
{
        double ctau;
        pumas_particle(NULL, &ctau, NULL);

        struct pumas_medium * rock;
        geometry_medium(context, state, &rock);

        pumas_recorder_create(&recorder, 0);
        context->recorder = recorder;

        reset_error();
        initialise_state();
        state->kinetic = 1E+03;

        double X, d, t0;
        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
        pumas_property_proper_time(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &t0);
        t0 /= TEST_ROCK_DENSITY;
        d = X / TEST_ROCK_DENSITY;
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_int_ge(recorder->length, 2);
        ck_assert_ptr_nonnull(recorder->first);

        struct pumas_frame * frame = recorder->first;
        ck_assert_ptr_eq(frame->medium, rock);
        ck_assert_double_eq(frame->state.charge, -1.);
        ck_assert_double_eq(frame->state.kinetic, 1E+03);
        ck_assert_double_eq(frame->state.distance, 0.);
        ck_assert_double_eq(frame->state.grammage, 0.);
        ck_assert_double_eq(frame->state.time, 0.);
        ck_assert_double_eq(frame->state.weight, 1.);
        ck_assert_double_eq(frame->state.position[0], 0.);
        ck_assert_double_eq(frame->state.position[1], 0.);
        ck_assert_double_eq(frame->state.position[2], 0.);
        ck_assert_double_eq(frame->state.direction[0], 0.);
        ck_assert_double_eq(frame->state.direction[1], 0.);
        ck_assert_double_eq(frame->state.direction[2], 1.);
        ck_assert_int_eq(frame->state.decayed, 0);
        ck_assert_int_eq(frame->event, PUMAS_EVENT_START);
        ck_assert_ptr_nonnull(frame->next);

        while (frame->next != NULL) frame = frame->next;
        ck_assert_ptr_eq(frame->medium, rock);
        ck_assert_double_eq(frame->state.charge, -1.);
        ck_assert_double_eq(frame->state.kinetic, 0.);
        ck_assert_double_le(frame->state.distance, d + FLT_EPSILON);
        ck_assert_double_le(frame->state.grammage, X + FLT_EPSILON);
        ck_assert_double_le(frame->state.time, t0 + FLT_EPSILON);
        ck_assert_double_eq(frame->state.position[0], 0.);
        ck_assert_double_eq(frame->state.position[1], 0.);
        ck_assert_double_le(frame->state.position[2], d + FLT_EPSILON);
        ck_assert_double_eq(frame->state.direction[0], 0.);
        ck_assert_double_eq(frame->state.direction[1], 0.);
        ck_assert_double_eq(frame->state.direction[2], 1.);
        ck_assert_int_eq(frame->state.decayed, 0);
        ck_assert_int_eq(frame->event,
            PUMAS_EVENT_LIMIT_KINETIC | PUMAS_EVENT_STOP);
        ck_assert_ptr_null(frame->next);

        pumas_recorder_destroy(&recorder);
        context->recorder = NULL;
}
END_TEST


START_TEST (test_hybrid_magnet)
{
        geometry.magnet[1] = 0.1;

        /* Compare the MC deflection and the analytical computation at low
         * energy, i.e. when no DEL */
        double tol = 5E-02;
        reset_error();
        initialise_state();
        state->kinetic = 1E+00;

        double X, d, phi0, phi1, phi, * u = state->direction;
        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_property_magnetic_rotation(0, state->kinetic, &phi0);
        pumas_property_magnetic_rotation(0, 0., &phi1);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_eq_tol(u[0], sin(phi), tol);
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_eq_tol(u[2], cos(phi), tol);

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Check a high energy case */
        reset_error();
        initialise_state();
        state->kinetic = 1E+03;

        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_property_magnetic_rotation(0, state->kinetic, &phi0);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_le(u[0], sin(phi) + tol);
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_ge(u[2], cos(phi) - tol);

        norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Erase the magnetic field */
        geometry.magnet[1] = 0.;
}
END_TEST


/* Fixtures for the detailed tests */
static void detailed_setup(void)
{
        /* Load the tau data and create a simulation context */
        load_muon();
        geometry.uniform = 1;
        pumas_context_create(&context, 0);
        context->medium = &geometry_medium;
        context->longitudinal = 0;
        context->scheme = PUMAS_SCHEME_DETAILED;
        context->random = &uniform01;
}

static void detailed_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_finalise();
}


START_TEST (test_detailed_straight)
{
        context->longitudinal = 1;

        int i;
        enum pumas_event event_data, * event;
        struct pumas_medium * media_data[2], ** media;
        struct pumas_medium * rock;
        geometry_medium(context, state, &rock);

        double ctau;
        pumas_particle(NULL, &ctau, NULL);

        /* Check the missing random engine case */
        reset_error();
        context->random = NULL;
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MISSING_RANDOM);
        context->random = &uniform01;

        /* Check the forward full path transport */
        double X, d, t0;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_property_grammage(
                    PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
                pumas_property_proper_time(
                    PUMAS_SCHEME_HYBRID, 0, state->kinetic, &t0);
                t0 /= TEST_ROCK_DENSITY;
                d = X / TEST_ROCK_DENSITY;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->kinetic, 0.);
                ck_assert_double_le(state->distance, d);
                ck_assert_double_le(state->grammage, X);
                ck_assert_double_le(state->time, t0);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], d + FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        /* Check the no event case but with limit values set */
        context->event = PUMAS_EVENT_NONE;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                context->kinetic_limit = 0.5E+03;
                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->kinetic_limit = 0.;
                context->distance_max = 0.5 * d;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->distance_max = 0.;
                context->grammage_max = 0.5 * X;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->grammage_max = 0.;
                context->time_max = 0.5 * t0;
                reset_error();
                initialise_state();
                state->kinetic = 1.;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.);
                context->time_max = 0.;
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        /* Check the forward kinetic limit */
        double X1, d1, t1;
        context->kinetic_limit = 0.5E+03;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;
        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, context->kinetic_limit, &X1);
        pumas_property_proper_time(
            PUMAS_SCHEME_HYBRID, 0, context->kinetic_limit, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->kinetic, context->kinetic_limit);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], d - d1 + FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->kinetic_limit = 0.;

        /* Check the forward grammage limit */
        double k1;
        X1 = 1E-02 * X;
        d1 = X1 / TEST_ROCK_DENSITY;
        context->grammage_max = X1;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 0, k1, &t1);
                t1 /= TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], d - d1 + FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the forward distance limit */
        d1 = 1E-02 * d;
        X1 = d1 * TEST_ROCK_DENSITY;
        context->distance_max = d1;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 0, k1, &t1);
                t1 /= TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], d - d1 + FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->distance_max = 0.;

        /* Check the forward time limit */
        t1 = 1E-02 * t0;
        context->time_max = t1;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                k1 = state->kinetic;
                pumas_property_grammage(PUMAS_SCHEME_HYBRID, 0, k1, &X1);
                d1 = X1 / TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq_tol(state->distance / (d - d1), 1., 0.5);
                ck_assert_double_eq_tol(state->grammage / (X - X1), 1., 0.5);
                ck_assert_double_eq_tol(state->time, t1, FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq_tol(
                    state->position[2], state->distance, FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->time_max = 0.;

        /* Check the backward transport with a kinetic limit */
        context->forward = 0;
        context->kinetic_limit = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;
        k1 = 0.5 * context->kinetic_limit;
        pumas_property_grammage(PUMAS_SCHEME_HYBRID, 0, k1, &X1);
        pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_ge(state->kinetic, context->kinetic_limit);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_ge(
                    state->position[2], -(d - d1) - FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->kinetic_limit = 0.;

        /* Check the backward transport with a grammage limit */
        context->grammage_max = X - X1;
        context->kinetic_limit = 0.75E+03; /* Not activated  */
        context->distance_max = 0.5 * (d - d1); /* Not activated  */
        context->time_max = 0.5 * (t0 - t1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->kinetic, k1);
                ck_assert_double_eq_tol(state->distance, d - d1, FLT_EPSILON);
                ck_assert_double_eq(state->grammage, X - X1);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq_tol(
                    state->position[2], -(d - d1), FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the backward transport with a distance limit */
        context->distance_max = d - d1;
        context->grammage_max = 0.5 * (X - X1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->kinetic, k1);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_eq_tol(
                    state->position[2], -(d - d1), FLT_EPSILON);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;

        /* Check the backward transport with a time limit */
        context->time_max = t0 - t1;
        context->distance_max = 0.5 * (d - d1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = k1;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->kinetic, k1);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_lt(state->position[2], 0.);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->grammage_max = 0.;
        context->time_max = 0.5 * (t0 - t1); /* Not activated */

        /* Check the case of an initial kinetic limit violation */
        context->kinetic_limit = 0.5E+03;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        context->kinetic_limit = 1E+03;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 0.5E+03;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 0.5E+03);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_KINETIC);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->kinetic_limit = 0.;

        /* Check the case of an initial grammage limit violation */
        context->grammage_max = 0.5 * X;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->grammage = X;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->grammage, X);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->grammage = X;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->grammage, X);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->grammage_max = 0.;

        /* Check the case of an initial distance limit violation */
        context->distance_max = 0.5 * d;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->distance = d;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->distance, d);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->distance = d;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->distance, d);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_DISTANCE);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->forward = 0;
        context->distance_max = 0.;

        /* Check the case of an initial time limit violation */
        context->time_max = 0.5 * t0;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E+03;
                state->time = t0;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->time, t0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->kinetic = 1E03;
                state->time = t0;
                pumas_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->kinetic, 1E+03);
                ck_assert_double_eq(state->time, t0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_TIME);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->forward = 1;
        context->time_max = 0.;
        context->event = PUMAS_EVENT_NONE;
        context->longitudinal = 0;
}
END_TEST


START_TEST (test_detailed_scattering)
{
        double ctau;
        pumas_particle(NULL, &ctau, NULL);

        /* Check the forward transport with scattering */
        double X, d, t0, * u = state->direction;

        reset_error();
        initialise_state();
        state->kinetic = 1E+03;

        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
        pumas_property_proper_time(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &t0);
        t0 /= TEST_ROCK_DENSITY;
        d = X / TEST_ROCK_DENSITY;

        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->charge, -1.);
        ck_assert_double_eq(state->kinetic, 0.);
        ck_assert_double_le(state->distance, d);
        ck_assert_double_le(state->grammage, X);
        ck_assert_double_le(state->time, t0);
        ck_assert_double_eq_tol(
            state->weight, exp(-state->time / ctau), FLT_EPSILON);
        ck_assert_double_le(state->position[2], d + FLT_EPSILON);
        ck_assert_int_eq(state->decayed, 0);

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Check the backward transport with scattering */
        context->forward = 0;
        context->kinetic_limit = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;

        double k1, X1, d1, t1;
        k1 = 1E-03;
        pumas_property_grammage(PUMAS_SCHEME_HYBRID, 0, k1, &X1);
        pumas_property_proper_time(PUMAS_SCHEME_HYBRID, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        reset_error();
        initialise_state();
        state->kinetic = k1;
        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->charge, -1.);
        ck_assert_double_ge(state->kinetic, context->kinetic_limit);
        ck_assert_double_eq_tol(state->distance / (d - d1), 1., 0.5);
        ck_assert_double_eq_tol(state->grammage / (X - X1), 1., 0.5);
        ck_assert_double_eq_tol(state->time / (t0 - t1), 1., 0.5);
        ck_assert_int_eq(state->decayed, 0);

        norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Restore the context status */
        context->event = PUMAS_EVENT_NONE;
        context->kinetic_limit = 0.;
        context->forward = 1;
}
END_TEST


START_TEST (test_detailed_magnet)
{
        context->longitudinal = 1;
        geometry.magnet[1] = 0.1;

        /* Compare the MC deflection and the analytical computation at low
         * energy, i.e. when no DEL */
        context->forward = 0;
        context->event = PUMAS_EVENT_LIMIT_KINETIC;
        context->kinetic_limit = 1E+01;

        double tol = 5E-02;
        reset_error();
        initialise_state();
        state->kinetic = 1E+00;

        double X, d, phi0, phi1, phi, * u = state->direction;
        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_property_magnetic_rotation(0, state->kinetic, &phi0);
        pumas_property_magnetic_rotation(0, context->kinetic_limit, &phi1);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_eq_tol(u[0], sin(phi), tol);
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_eq_tol(u[2], cos(phi), tol);

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Check a high energy case */
        context->forward = 1;
        context->event = PUMAS_EVENT_NONE;
        context->kinetic_limit = 0.;

        reset_error();
        initialise_state();
        state->kinetic = 1E+03;

        pumas_property_grammage(
            PUMAS_SCHEME_HYBRID, 0, state->kinetic, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_property_magnetic_rotation(0, state->kinetic, &phi0);
        pumas_property_magnetic_rotation(0, 0., &phi1);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_le(u[0], sin(phi) + tol);
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_ge(u[2], cos(phi) - tol);

        norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Restore the context and the geometry */
        geometry.magnet[1] = 0.;
        context->longitudinal = 0;
}
END_TEST


/* Fixtures for tau tests */
static void tau_setup(void)
{
        /* Load the tau data and create a simulation context */
        load_tau();
        geometry.uniform = 1;
        pumas_context_create(&context, 0);
        context->medium = &geometry_medium;
        context->random = &uniform01;
}

static void tau_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_finalise();
}


START_TEST (test_tau_csda)
{
        enum pumas_event event;
        context->longitudinal = 1;
        context->scheme = PUMAS_SCHEME_CSDA;

        /* Test some errors */
        context->decay = PUMAS_DECAY_WEIGHT;
        reset_error();
        initialise_state();
        pumas_transport(context, state, &event, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_DECAY_ERROR);
        context->decay = PUMAS_DECAY_PROCESS;

        /* Test the forward transport */
        reset_error();
        initialise_state();
        state->kinetic = 1.;

        double X, d;
        pumas_property_grammage(
            PUMAS_SCHEME_CSDA, 0, state->kinetic, &X);
        d = X / TEST_ROCK_DENSITY;

        pumas_transport(context, state, &event, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(event, PUMAS_EVENT_VERTEX_DECAY);

        ck_assert_double_eq(state->direction[0], 0.);
        ck_assert_double_eq(state->direction[1], 0.);
        ck_assert_double_eq(state->direction[2], 1.);

        context->longitudinal = 0;
        context->scheme = PUMAS_SCHEME_DETAILED;
}
END_TEST


Suite * create_suite(void)
{
        const int timeout = 60;

        /* The test suite */
        Suite * suite = suite_create("Pumas");

        /* The API test case */
        TCase * tc_api = tcase_create("API");
        suite_add_tcase(suite, tc_api);
        tcase_set_timeout(tc_api, timeout);
        tcase_add_test(tc_api, test_api_error);
        tcase_add_test(tc_api, test_api_tag);
        tcase_add_test(tc_api, test_api_init);
        tcase_add_test(tc_api, test_api_memory);
        tcase_add_test(tc_api, test_api_material);
        tcase_add_test(tc_api, test_api_composite);
        tcase_add_test(tc_api, test_api_property);
        tcase_add_test(tc_api, test_api_table);
        tcase_add_test(tc_api, test_api_context);
        tcase_add_test(tc_api, test_api_recorder);
        tcase_add_test(tc_api, test_api_print);

        /* The no loss test case */
        TCase * tc_lossless = tcase_create("Lossless");
        suite_add_tcase(suite, tc_lossless);
        tcase_set_timeout(tc_lossless, timeout);
        tcase_add_unchecked_fixture(
            tc_lossless, lossless_setup, lossless_teardown);
        tcase_add_test(tc_lossless, test_lossless_straight);
        tcase_add_test(tc_lossless, test_lossless_magnet);
        tcase_add_test(tc_lossless, test_lossless_geometry);

        /* The CSDA test case */
        TCase * tc_csda = tcase_create("CSDA");
        suite_add_tcase(suite, tc_csda);
        tcase_set_timeout(tc_csda, timeout);
        tcase_add_unchecked_fixture(tc_csda, csda_setup, csda_teardown);
        tcase_add_test(tc_csda, test_csda_straight);
        tcase_add_test(tc_csda, test_csda_record);
        tcase_add_test(tc_csda, test_csda_magnet);
        tcase_add_test(tc_csda, test_csda_geometry);

        /* The hybrid test case */
        TCase * tc_hybrid = tcase_create("Hybrid");
        suite_add_tcase(suite, tc_hybrid);
        tcase_set_timeout(tc_hybrid, timeout);
        tcase_add_unchecked_fixture(tc_hybrid, hybrid_setup, hybrid_teardown);
        tcase_add_test(tc_hybrid, test_hybrid_straight);
        tcase_add_test(tc_hybrid, test_hybrid_scattering);
        tcase_add_test(tc_hybrid, test_hybrid_record);
        tcase_add_test(tc_hybrid, test_hybrid_magnet);

        /* The detailed test case */
        TCase * tc_detailed = tcase_create("Detailed");
        suite_add_tcase(suite, tc_detailed);
        tcase_set_timeout(tc_detailed, timeout);
        tcase_add_unchecked_fixture(
            tc_detailed, detailed_setup, detailed_teardown);
        tcase_add_test(tc_detailed, test_detailed_straight);
        tcase_add_test(tc_detailed, test_detailed_scattering);
        tcase_add_test(tc_detailed, test_detailed_magnet);

        /* The tau test case */
        TCase * tc_tau = tcase_create("Tau");
        suite_add_tcase(suite, tc_tau);
        tcase_set_timeout(tc_tau, timeout);
        tcase_add_unchecked_fixture(tc_tau, tau_setup, tau_teardown);
        tcase_add_test(tc_tau, test_tau_csda);

        return suite;
}


int main(void)
{
        /* Initialise the PRNG */
        srand(0);

        /* Configure the tests and the runner */
        Suite * suite = create_suite();
        SRunner * runner = srunner_create(suite);
        srunner_set_fork_status(runner, CK_NOFORK);

        /* Run the tests */
        srunner_run_all(runner, CK_NORMAL);
        const int status = srunner_ntests_failed(runner) ?
            EXIT_FAILURE : EXIT_SUCCESS;
        srunner_free(runner);

        /* Return the test status to the OS */
        exit(status);
}
