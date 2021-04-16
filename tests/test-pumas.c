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
static struct pumas_physics * physics = NULL;
static struct pumas_context * context = NULL;
static struct pumas_recorder * recorder = NULL;
static struct pumas_state state_data, *state = &state_data;

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

/* Test the constant API */
START_TEST(test_api_constant)
{
        double value;

        pumas_constant(PUMAS_CONSTANT_AVOGADRO_NUMBER, &value);
        ck_assert_double_eq(value, 6.02214076E+23);

        pumas_constant(PUMAS_CONSTANT_ELECTRON_MASS, &value);
        ck_assert_double_eq(value, 0.510998910E-03);

        pumas_constant(PUMAS_CONSTANT_MUON_C_TAU, &value);
        ck_assert_double_eq(value, 658.654);

        pumas_constant(PUMAS_CONSTANT_MUON_MASS, &value);
        ck_assert_double_eq(value, 0.10565839);

        pumas_constant(PUMAS_CONSTANT_NEUTRON_MASS, &value);
        ck_assert_double_eq(value, 0.939565);

        pumas_constant(PUMAS_CONSTANT_PROTON_MASS, &value);
        ck_assert_double_eq(value, 0.938272);

        pumas_constant(PUMAS_CONSTANT_TAU_C_TAU, &value);
        ck_assert_double_eq(value, 87.03E-06);

        pumas_constant(PUMAS_CONSTANT_TAU_MASS, &value);
        ck_assert_double_eq(value, 1.77682);

        reset_error();
        pumas_constant(-1, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_constant(8, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_constant(PUMAS_CONSTANT_AVOGADRO_NUMBER, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);
}
END_TEST

/* Test the error API */
START_TEST(test_api_error)
{
        /* Check the setters & getters  */
        pumas_error_handler_set(NULL);
        pumas_handler_cb * handler = pumas_error_handler_get();
        ck_assert_ptr_null(handler);

        pumas_error_handler_set(&handle_error);
        handler = pumas_error_handler_get();
        ck_assert_ptr_eq(handler, &handle_error);

        /* Check the error catch & raise */
        reset_error();
        pumas_error_catch(1);
        pumas_physics_material_index(physics, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(error_data.caller);

        enum pumas_return rc = pumas_error_raise();
        ck_assert_int_eq(rc, PUMAS_RETURN_PHYSICS_ERROR);
        ck_assert_int_eq(error_data.rc, rc);
        ck_assert_ptr_eq(error_data.caller, &pumas_physics_material_index);

        /* Check that errors are enabled again */
        reset_error();
        pumas_physics_material_index(physics, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        /* Check that error catching can be disabled */
        reset_error();
        pumas_error_catch(1);
        pumas_physics_material_index(physics, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        pumas_error_catch(0);
        pumas_physics_material_index(physics, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

/* Check the stringification of API functions */
#define CHECK_STRING(function)                                                 \
        ck_assert_str_eq(                                                      \
            pumas_error_function((pumas_function_t *)&function), #function)
        /* Generated with:
         * nm -D --defined-only libpumas.so |                                  \
         *     grep "T " |                                                     \
         *     awk '{print "CHECK_STRING("$3");"}'
         */
        CHECK_STRING(pumas_constant);
        CHECK_STRING(pumas_context_create);
        CHECK_STRING(pumas_context_random_dump);
        CHECK_STRING(pumas_context_random_load);
        CHECK_STRING(pumas_context_random_seed_get);
        CHECK_STRING(pumas_context_random_seed_set);
        CHECK_STRING(pumas_context_destroy);
        CHECK_STRING(pumas_context_physics_get);
        CHECK_STRING(pumas_context_transport);
        CHECK_STRING(pumas_dcs_default);
        CHECK_STRING(pumas_dcs_get);
        CHECK_STRING(pumas_dcs_register);
        CHECK_STRING(pumas_error_catch);
        CHECK_STRING(pumas_error_function);
        CHECK_STRING(pumas_error_handler_get);
        CHECK_STRING(pumas_error_handler_set);
        CHECK_STRING(pumas_error_raise);
        CHECK_STRING(pumas_memory_allocator);
        CHECK_STRING(pumas_memory_deallocator);
        CHECK_STRING(pumas_memory_reallocator);
        CHECK_STRING(pumas_physics_composite_length);
        CHECK_STRING(pumas_physics_composite_properties);
        CHECK_STRING(pumas_physics_composite_update);
        CHECK_STRING(pumas_physics_create);
        CHECK_STRING(pumas_physics_cutoff);
        CHECK_STRING(pumas_physics_dcs);
        CHECK_STRING(pumas_physics_destroy);
        CHECK_STRING(pumas_physics_dump);
        CHECK_STRING(pumas_physics_load);
        CHECK_STRING(pumas_physics_material_index);
        CHECK_STRING(pumas_physics_material_length);
        CHECK_STRING(pumas_physics_material_name);
        CHECK_STRING(pumas_physics_particle);
        CHECK_STRING(pumas_physics_print);
        CHECK_STRING(pumas_physics_property_cross_section);
        CHECK_STRING(pumas_physics_property_energy_loss);
        CHECK_STRING(pumas_physics_property_grammage);
        CHECK_STRING(pumas_physics_property_kinetic_energy);
        CHECK_STRING(pumas_physics_property_magnetic_rotation);
        CHECK_STRING(pumas_physics_property_proper_time);
        CHECK_STRING(pumas_physics_property_scattering_length);
        CHECK_STRING(pumas_physics_table_index);
        CHECK_STRING(pumas_physics_table_length);
        CHECK_STRING(pumas_physics_table_value);
        CHECK_STRING(pumas_recorder_clear);
        CHECK_STRING(pumas_recorder_create);
        CHECK_STRING(pumas_recorder_destroy);
        CHECK_STRING(pumas_version);
#undef CHECK_STRING

        /* Check the NULL case */
        ck_assert_ptr_null(
            pumas_error_function((pumas_function_t *)&test_api_error));
}
END_TEST

/* Test the version tag */
START_TEST(test_api_version)
{
        const int tag = pumas_version();
        ck_assert_int_eq(tag, 100 * PUMAS_VERSION);
}
END_TEST

static void load_muon(void)
{
#define TEST_MUON_DUMP ".pumas.muon"
        FILE * fid = fopen(TEST_MUON_DUMP, "rb");
        pumas_physics_load(&physics, fid);
        fclose(fid);
}

static void dump_muon(void)
{
        FILE * fid = fopen(TEST_MUON_DUMP, "wb+");
        pumas_physics_dump(physics, fid);
        fclose(fid);
}

static void load_tau(void)
{
#define TEST_TAU_DUMP ".pumas.tau"
        FILE * fid = fopen(TEST_TAU_DUMP, "rb");
        pumas_physics_load(&physics, fid);
        fclose(fid);
}

static void dump_tau(void)
{
        FILE * fid = fopen(TEST_TAU_DUMP, "wb+");
        pumas_physics_dump(physics, fid);
        fclose(fid);
}

/* Test the library initialisation & finalisation */
START_TEST(test_api_init)
{
        /* Check the particle selection */
        reset_error();
        pumas_physics_create(&physics, -1, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_UNKNOWN_PARTICLE);

        /* Check the missing MDF error */
        reset_error();
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_UNDEFINED_MDF);

        /* Check the load error */
        reset_error();
        pumas_physics_load(&physics, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PATH_ERROR);

        /* Check the cutoff error */
        reset_error();
        struct pumas_physics_settings settings = {.cutoff = 2};
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", "materials/dedx/muon",
            &settings);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_CUTOFF_ERROR);

        /* Check the default cutoff */
        reset_error();
        settings.cutoff = 0;
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", "materials/dedx/muon",
            &settings);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(5E-02, pumas_physics_cutoff(physics));
        pumas_physics_destroy(&physics);

        /* Check the cutoff setter */
        reset_error();
        settings.cutoff = 0.1;
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", "materials/dedx/muon",
            &settings);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(settings.cutoff, pumas_physics_cutoff(physics));
        pumas_physics_destroy(&physics);

        /* Check the initialisation and dump */
        reset_error();
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", "materials/dedx/muon", NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        reset_error();
        dump_muon();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        pumas_physics_destroy(&physics);

        reset_error();
        pumas_physics_create(&physics, PUMAS_PARTICLE_TAU,
            "materials/materials.xml", "materials/dedx/tau", NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        reset_error();
        dump_tau();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        pumas_physics_destroy(&physics);

        /* Check the load */
        reset_error();
        load_tau();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        enum pumas_particle particle;
        double lifetime, mass;

        reset_error();
        pumas_physics_particle(physics, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        pumas_physics_particle(physics, &particle, &lifetime, &mass);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(particle, PUMAS_PARTICLE_TAU);
        ck_assert_double_eq(lifetime, 87.03E-06);
        ck_assert_double_eq(mass, 1.77682);

        pumas_physics_destroy(&physics);
        pumas_physics_particle(physics, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        load_muon();
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        reset_error();
        pumas_physics_particle(physics, &particle, &lifetime, &mass);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(particle, PUMAS_PARTICLE_MUON);
        ck_assert_double_eq(lifetime, 658.654);
        ck_assert_double_eq(mass, 0.10565839);

        /* Check the dump errors */
        reset_error();
        pumas_physics_dump(physics, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PATH_ERROR);
        pumas_physics_destroy(&physics);

        reset_error();
        pumas_physics_dump(physics, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        /* Check the transport error */
        reset_error();
        pumas_context_transport(NULL, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);
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
        if (ptr == NULL) dyn_mem_count++;
        return realloc(ptr, size);
}

static void test_free(void * ptr)
{
        if (ptr != NULL) dyn_mem_count--;
        free(ptr);
}

/* Test the memory API */
START_TEST(test_api_memory)
{
        /* Set the test memory routines */
        pumas_memory_allocator(&test_malloc);
        pumas_memory_reallocator(&test_realloc);
        pumas_memory_deallocator(&test_free);

        /* Load the muon data */
        load_muon();
        ck_assert_int_gt(dyn_mem_count, 0);

        /* Release the memory */
        pumas_physics_destroy(&physics);
        ck_assert_int_eq(dyn_mem_count, 0);

        /* Set the system memory routines */
        pumas_memory_allocator(NULL);
        pumas_memory_reallocator(NULL);
        pumas_memory_deallocator(NULL);
}
END_TEST

/* Test the element API */
START_TEST(test_api_element)
{
        int i, index, length;
        const char * name;
        double Z, A, I;

        /* Check the initialisation error */
        reset_error();
        pumas_physics_element_index(physics, "H", &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_element_name(physics, 0, &name);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        index = pumas_physics_element_length(physics);
        ck_assert_int_eq(index, 0);

        reset_error();
        pumas_physics_element_properties(
            physics, 0, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        /* Load the muon data */
        load_muon();

        /* Check the number of elements */
        index = pumas_physics_element_length(physics);
        ck_assert_int_eq(index, 6);

        /* Check the invalid material case */
        reset_error();
        pumas_physics_element_index(physics, "Strawberry", &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_UNKNOWN_ELEMENT);

        reset_error();
        pumas_physics_element_name(physics, 6, &name);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_element_properties(physics, 6, &Z, &A, &I);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Check the elements mapping */
        const char * names[] = {"H", "C", "N", "O", "Rk", "Ar"};
        for (i = 0; i < sizeof(names) / sizeof(*names); i++) {
                reset_error();
                pumas_physics_element_index(physics, names[i], &index);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(index, i);

                reset_error();
                pumas_physics_element_name(physics, i, &name);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_str_eq(name, names[i]);
        }

        /* Check the properties getter */
        reset_error();
        pumas_physics_element_properties(physics, 0, &Z, &A, &I);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(Z, 1);
        ck_assert_double_eq(A, 1.0087);
        ck_assert_double_eq_tol(I * 1E+09, 19.2, 1E-06);

        pumas_physics_destroy(&physics);
}
END_TEST

/* Test the materials API */
START_TEST(test_api_material)
{
        int i, index, length, components[2];
        const char * name;
        const char * names[] = { "StandardRock", "Water", "Air" };
        double density, I, fractions[2];

        /* Check the initialisation error */
        reset_error();
        pumas_physics_material_index(physics, "Strawberry", &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_material_name(physics, 0, &name);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        index = pumas_physics_material_length(physics);
        ck_assert_int_eq(index, 0);

        reset_error();
        pumas_physics_material_properties(
            physics, 0, NULL, NULL, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        /* Load the muon data */
        load_muon();

        /* Check the number of materials and composites */
        index = pumas_physics_material_length(physics);
        ck_assert_int_eq(index, 4);

        index = pumas_physics_composite_length(physics);
        ck_assert_int_eq(index, 1);

        /* Check the invalid material case */
        reset_error();
        pumas_physics_material_index(physics, "Strawberry", &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_UNKNOWN_MATERIAL);

        reset_error();
        pumas_physics_material_name(physics, 4, &name);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_material_properties(physics, 4, &length, &density, &I,
            components, fractions);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Check the materials mapping */
        for (i = 0; i < sizeof(names) / sizeof(*names); i++) {
                reset_error();
                pumas_physics_material_index(physics, names[i], &index);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(index, i);

                reset_error();
                pumas_physics_material_name(physics, i, &name);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_str_eq(name, names[i]);
        }

        /* Check the properties getter */
        double wH = 0.111894, wO = 0.888106;

        reset_error();
        pumas_physics_material_properties(physics, 0, &length, &density, &I,
            components, fractions);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(length, 1);
        ck_assert_double_eq(density, 2.65E+03);
        ck_assert_double_eq_tol(I, 1.364E-07, 1E-11);
        ck_assert_double_eq(components[0], 4);
        ck_assert_double_eq(fractions[0], 1);

        reset_error();
        pumas_physics_material_properties(physics, 1, &length, &density, &I,
            components, fractions);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(length, 2);
        ck_assert_double_eq(density, 1E+03);
        ck_assert_double_eq_tol(I, 7.97E-08, 1E-11);
        ck_assert_double_eq(components[0], 3);
        ck_assert_double_eq(components[1], 0);
        ck_assert_double_eq(fractions[0], wO);
        ck_assert_double_eq(fractions[1], wH);

        pumas_physics_destroy(&physics);
}
END_TEST

/* Test the composites API */
START_TEST(test_api_composite)
{
        int i, index, length, components[3];
        const char * name;
        const char * names[] = { "StandardRock", "Water", "Air", "WetRock" };
        double density, fractions[3];

        /* Check the initialisation error */
        reset_error();
        pumas_physics_composite_properties(physics, 0, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_composite_update(physics, 0, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        index = pumas_physics_composite_length(physics);
        ck_assert_int_eq(index, 0);

        /* Load the tau data */
        load_tau();

        /* Check the number of materials and composites */
        index = pumas_physics_material_length(physics);
        ck_assert_int_eq(index, 4);

        index = pumas_physics_composite_length(physics);
        ck_assert_int_eq(index, 1);

        /* Check the materials mapping */
        for (i = 0; i < sizeof(names) / sizeof(*names); i++) {
                reset_error();
                pumas_physics_material_index(physics, names[i], &index);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_int_eq(index, i);

                reset_error();
                pumas_physics_material_name(physics, i, &name);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_str_eq(name, names[i]);
        }

        /* Check the properties getter */
        reset_error();
        pumas_physics_composite_properties(physics, 0, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_composite_properties(physics, 3, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        double density0 = TEST_ROCK_DENSITY, density1 = 1.00E+03;
        double fraction0 = 0.5, fraction1 = 0.5;
        double rho = 1. / (fraction0 / density0 + fraction1 / density1);
        double wH = 0.111894, wO = 0.888106;

        reset_error();
        pumas_physics_composite_properties(
            physics, 3, &length, components, fractions);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(length, 2);
        ck_assert_int_eq(components[0], 0);
        ck_assert_int_eq(components[1], 1);
        ck_assert_double_eq(fractions[0], fraction0);
        ck_assert_double_eq(fractions[1], fraction1);

        reset_error();
        pumas_physics_material_properties(
            physics, 3, &length, &density, NULL, components, fractions);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(length, 3);
        ck_assert_double_eq(density, rho);
        ck_assert_int_eq(components[0], 4);
        ck_assert_int_eq(components[1], 3);
        ck_assert_int_eq(components[2], 0);
        ck_assert_double_eq(fractions[0], fraction0);
        ck_assert_double_eq(fractions[1], fraction1 * wO);
        ck_assert_double_eq(fractions[2], fraction1 * wH);

        /* Check the properties setter */
        reset_error();
        pumas_physics_composite_update(physics, 2, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        pumas_physics_composite_properties(
            physics, 3, &length, components, fractions);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(length, 2);
        ck_assert_int_eq(components[0], 0);
        ck_assert_int_eq(components[1], 1);
        ck_assert_double_eq(fractions[0], fraction0);
        ck_assert_double_eq(fractions[1], fraction1);

        reset_error();
        pumas_physics_material_properties(
            physics, 3, NULL, &density, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(density, rho);

        fraction0 = 0.3, fraction1 = 0.7;
        rho = 1. / (fraction0 / density0 + fraction1 / density1);

        {
                double tmp[] = { fraction0, fraction1 };
                pumas_physics_composite_update(physics, 3, tmp);
        }
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        reset_error();
        pumas_physics_composite_properties(
            physics, 3, &length, components, fractions);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(length, 2);
        ck_assert_int_eq(components[0], 0);
        ck_assert_int_eq(components[1], 1);
        ck_assert_double_eq(fractions[0], fraction0);
        ck_assert_double_eq(fractions[1], fraction1);

        reset_error();
        pumas_physics_material_properties(
            physics, 3, &length, &density, NULL, components, fractions);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(length, 3);
        ck_assert_double_eq(density, rho);
        ck_assert_int_eq(components[0], 4);
        ck_assert_int_eq(components[1], 3);
        ck_assert_int_eq(components[2], 0);
        ck_assert_double_eq(fractions[0], fraction0);
        ck_assert_double_eq(fractions[1], fraction1 * wO);
        ck_assert_double_eq(fractions[2], fraction1 * wH);

        pumas_physics_destroy(&physics);
}
END_TEST

START_TEST(test_api_property)
{
        double value;

        /* Check the initialisation error */
        reset_error();
        pumas_physics_property_cross_section(physics, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_property_energy_loss(physics, 0, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_property_grammage(physics, 0, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_property_kinetic_energy(physics, 0, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_property_magnetic_rotation(physics, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_property_proper_time(physics, 0, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_property_scattering_length(physics, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        /* Load the muon data */
        load_muon();

        /* Check the index error */
        reset_error();
        pumas_physics_property_cross_section(physics, 4, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_energy_loss(
            physics, PUMAS_MODE_HYBRID, 4, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_energy_loss(
            physics, PUMAS_MODE_DETAILED, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 4, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_DETAILED, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_kinetic_energy(
            physics, PUMAS_MODE_HYBRID, 4, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_kinetic_energy(
            physics, PUMAS_MODE_DETAILED, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_magnetic_rotation(physics, 4, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 4, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_DETAILED, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_property_scattering_length(physics, 4, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Check some values */
        reset_error();
        pumas_physics_property_cross_section(physics, 0, 0., &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_physics_property_cross_section(physics, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_lt(1. / (value * TEST_ROCK_DENSITY), 1E+03);

        pumas_physics_property_cross_section(physics, 0, 1E+06, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_lt(1. / (value * TEST_ROCK_DENSITY), 1E+03);

        pumas_physics_property_energy_loss(
            physics, PUMAS_MODE_CSDA, 0, 1E+00, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 1.823E-04, 1E-07);

        pumas_physics_property_energy_loss(
            physics, PUMAS_MODE_CSDA, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 6.623E-04, 1E-07);

        pumas_physics_property_energy_loss(
            physics, PUMAS_MODE_HYBRID, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 6.623E-04);

        pumas_physics_property_grammage(
            physics, PUMAS_MODE_CSDA, 0, 1E+00, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 5.475E+03, 1.);

        pumas_physics_property_grammage(
            physics, PUMAS_MODE_CSDA, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 2.449E+06, 1E+03);

        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, 1E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 2.449E+06);

        pumas_physics_property_kinetic_energy(
            physics, PUMAS_MODE_CSDA, 0, 5.475E+03, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 1E+00, 1E-03);

        pumas_physics_property_kinetic_energy(
            physics, PUMAS_MODE_CSDA, 0, 2.449E+06, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 1E+03, 1.);

        pumas_physics_property_kinetic_energy(
            physics, PUMAS_MODE_HYBRID, 0, 2.449E+06, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 1E+03);

        pumas_physics_property_magnetic_rotation(physics, 0, 1E+00, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_gt(value, 5.475E+03 / 3.3374);

        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_CSDA, 0, 1E+00, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 5.475E+03);

        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_CSDA, 0, 1E-02, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 8.498);

        pumas_physics_property_scattering_length(physics, 0, 1E-02, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        value *= 1E-03; /* Normalise to electrons in Al */
        ck_assert_double_gt(value, 0.3);
        ck_assert_double_lt(value, 3.);

        /* Check overflows */
        double vmax;
        pumas_physics_property_energy_loss(
            physics, PUMAS_MODE_CSDA, 0, 1E+08, &value);
        pumas_physics_table_value(physics, PUMAS_PROPERTY_ENERGY_LOSS,
            PUMAS_MODE_CSDA, 0, 145, &vmax);
        ck_assert_double_gt(value, vmax);

        pumas_physics_property_grammage(
            physics, PUMAS_MODE_CSDA, 0, 1E+08, &value);
        pumas_physics_table_value(
            physics, PUMAS_PROPERTY_GRAMMAGE, PUMAS_MODE_CSDA, 0, 145, &vmax);
        ck_assert_double_gt(value, vmax);

        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_CSDA, 0, 1E+08, &value);
        pumas_physics_table_value(physics, PUMAS_PROPERTY_PROPER_TIME,
            PUMAS_MODE_CSDA, 0, 145, &vmax);
        ck_assert_double_gt(value, vmax);

        pumas_physics_property_cross_section(physics, 0, 1E+08, &value);
        pumas_physics_table_value(physics, PUMAS_PROPERTY_CROSS_SECTION,
            PUMAS_MODE_HYBRID, 0, 145, &vmax);
        ck_assert_double_eq(value, vmax);

        /* Unload the data */
        pumas_physics_destroy(&physics);
}
END_TEST

START_TEST(test_api_table)
{
        int index;
        double value;

        /* Check the initialisation error */
        reset_error();
        pumas_physics_table_index(physics, 0, 0, 0, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        pumas_physics_table_value(physics, 0, 0, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        reset_error();
        index = pumas_physics_table_length(physics);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 0);

        /* Load the muon data */
        load_muon();

        reset_error();
        index = pumas_physics_table_length(physics);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 146);

        /* Check the index error */
        reset_error();
        pumas_physics_table_index(physics, 1024, 0, 0, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_table_index(physics, PUMAS_PROPERTY_GRAMMAGE,
            PUMAS_MODE_DETAILED, 0, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_table_index(
            physics, PUMAS_PROPERTY_GRAMMAGE, PUMAS_MODE_CSDA, 4, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_table_index(physics, PUMAS_PROPERTY_GRAMMAGE,
            PUMAS_MODE_CSDA, 0, 1E+18, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        pumas_physics_table_index(physics, PUMAS_PROPERTY_GRAMMAGE,
            PUMAS_MODE_CSDA, 0, -1., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        pumas_physics_table_value(physics, 1024, 0, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_table_value(physics, PUMAS_PROPERTY_GRAMMAGE,
            PUMAS_MODE_DETAILED, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_table_value(
            physics, PUMAS_PROPERTY_GRAMMAGE, PUMAS_MODE_CSDA, 4, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_table_value(physics, PUMAS_PROPERTY_GRAMMAGE,
            PUMAS_MODE_CSDA, 0, 200, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_table_value(
            physics, PUMAS_PROPERTY_GRAMMAGE, PUMAS_MODE_CSDA, 0, -1, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Check some values */
        reset_error();

        pumas_physics_table_index(physics, PUMAS_PROPERTY_KINETIC_ENERGY,
            PUMAS_MODE_CSDA, 0, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 0);

        pumas_physics_table_index(physics, PUMAS_PROPERTY_KINETIC_ENERGY,
            PUMAS_MODE_CSDA, 0, 1E+00, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_physics_table_index(physics, PUMAS_PROPERTY_KINETIC_ENERGY,
            PUMAS_MODE_CSDA, 0, 1E+06, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 145);

        pumas_physics_table_value(physics, PUMAS_PROPERTY_KINETIC_ENERGY,
            PUMAS_MODE_CSDA, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_physics_table_value(physics, PUMAS_PROPERTY_KINETIC_ENERGY,
            PUMAS_MODE_CSDA, 0, 49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 1E+00);

        pumas_physics_table_index(
            physics, PUMAS_PROPERTY_GRAMMAGE, PUMAS_MODE_CSDA, 0, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 0);

        pumas_physics_table_index(physics, PUMAS_PROPERTY_GRAMMAGE,
            PUMAS_MODE_CSDA, 0, 5.476E+03, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_physics_table_value(
            physics, PUMAS_PROPERTY_GRAMMAGE, PUMAS_MODE_CSDA, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_physics_table_value(
            physics, PUMAS_PROPERTY_GRAMMAGE, PUMAS_MODE_CSDA, 0, 49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq_tol(value, 5.476E+03, 1.);

        pumas_physics_table_value(physics, PUMAS_PROPERTY_GRAMMAGE,
            PUMAS_MODE_HYBRID, 0, 49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 5.476E+03);

        pumas_physics_table_value(physics, PUMAS_PROPERTY_ENERGY_LOSS,
            PUMAS_MODE_CSDA, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_physics_table_value(physics, PUMAS_PROPERTY_CROSS_SECTION,
            PUMAS_MODE_HYBRID, 0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_physics_table_value(
            physics, PUMAS_PROPERTY_MAGNETIC_ROTATION, 0, 0, 49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_gt(value, 5.476E+03 / 3.3374);

        pumas_physics_table_index(
            physics, PUMAS_PROPERTY_MAGNETIC_ROTATION, 0, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_physics_table_value(physics, PUMAS_PROPERTY_PROPER_TIME,
            PUMAS_MODE_CSDA, 0, 49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 5.476E+03);

        pumas_physics_table_index(physics, PUMAS_PROPERTY_PROPER_TIME,
            PUMAS_MODE_CSDA, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_physics_table_value(physics, PUMAS_PROPERTY_PROPER_TIME,
            PUMAS_MODE_CSDA, 0, 17, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 8.498);

        pumas_physics_table_index(physics, PUMAS_PROPERTY_PROPER_TIME,
            PUMAS_MODE_CSDA, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 17);

        /* Check unsuported properties  */
        pumas_physics_table_value(
            physics, PUMAS_PROPERTY_SCATTERING_LENGTH, 0, 0, 17, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        pumas_physics_table_index(
            physics, PUMAS_PROPERTY_CROSS_SECTION, 0, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        pumas_physics_table_index(
            physics, PUMAS_PROPERTY_ENERGY_LOSS, 0, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        pumas_physics_table_index(
            physics, PUMAS_PROPERTY_SCATTERING_LENGTH, 0, 0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Unload the data */
        pumas_physics_destroy(&physics);
}
END_TEST

static void * fail_malloc(size_t size) { return NULL; }

static void * fail_realloc(void * ptr, size_t size) { return NULL; }

/* Test the context API */
START_TEST(test_api_context)
{
        /* Test the initialisation error */
        reset_error();
        context = (void *)0x1;
        pumas_context_create(&context, physics, 0);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);
        ck_assert_ptr_null(context);

        /* Load the muon data */
        load_muon();

        /* Test the Physics getter in the NULL case */
        ck_assert_ptr_eq(pumas_context_physics_get(context), NULL);

        /* Test a memory error */
        pumas_memory_allocator(&fail_malloc);
        reset_error();
        context = (void *)0x1;
        pumas_context_create(&context, physics, 0);
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
        int i, *data;

        pumas_context_create(&context, physics, n * sizeof *data);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(context->user_data);

        data = context->user_data;
        for (i = 0; i < n; i++) data[i] = i;

        /* Test the initialisation of a muon context */
        ck_assert_ptr_eq(pumas_context_physics_get(context), physics);
        ck_assert_ptr_null(context->medium);
        ck_assert_ptr_nonnull(context->random);
        ck_assert_ptr_null(context->recorder);
        ck_assert_int_eq(context->mode.scattering, PUMAS_MODE_FULL_SPACE);
        ck_assert_int_eq(context->mode.direction, PUMAS_MODE_FORWARD);
        ck_assert_int_eq(context->mode.energy_loss, PUMAS_MODE_DETAILED);
        ck_assert_int_eq(context->mode.decay, PUMAS_MODE_WEIGHT);
        ck_assert_int_eq(context->event, PUMAS_EVENT_NONE);

        ck_assert_double_eq(context->limit.energy, 0.);
        ck_assert_double_eq(context->limit.distance, 0.);
        ck_assert_double_eq(context->limit.grammage, 0.);
        ck_assert_double_eq(context->limit.time, 0.);

        ck_assert_double_eq(context->accuracy, 1E-02);

        /* Test the context destruction */
        pumas_context_destroy(&context);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(context);

        /* Load the tau data */
        pumas_physics_destroy(&physics);
        load_tau();

        /* Test the initialisation of a tau context */
        pumas_context_create(&context, physics, -1);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_null(context->user_data);

        ck_assert_ptr_eq(pumas_context_physics_get(context), physics);
        ck_assert_ptr_null(context->medium);
        ck_assert_ptr_nonnull(context->random);
        ck_assert_ptr_null(context->recorder);
        ck_assert_int_eq(context->mode.scattering, PUMAS_MODE_FULL_SPACE);
        ck_assert_int_eq(context->mode.direction, PUMAS_MODE_FORWARD);
        ck_assert_int_eq(context->mode.energy_loss, PUMAS_MODE_DETAILED);
        ck_assert_int_eq(context->mode.decay, PUMAS_MODE_DECAY);
        ck_assert_int_eq(context->event, PUMAS_EVENT_NONE);

        ck_assert_double_eq(context->limit.energy, 0.);
        ck_assert_double_eq(context->limit.distance, 0.);
        ck_assert_double_eq(context->limit.grammage, 0.);
        ck_assert_double_eq(context->limit.time, 0.);

        /* Free the data */
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);
}
END_TEST

/* Test the random API */
START_TEST(test_api_random)
{
        /* Load the muon data */
        load_muon();

        /* Test the random seed init */
        pumas_context_create(&context, physics, 0);
        unsigned long seed = 1;
        pumas_context_random_seed_set(context, &seed);
        const double u1 = context->random(context);
        pumas_context_random_seed_set(context, &seed);
        ck_assert_double_eq(context->random(context), u1);

        pumas_context_random_seed_get(context, &seed);
        ck_assert_uint_eq(seed, 1);

        pumas_context_destroy(&context);
        pumas_context_create(&context, physics, 0);
        ck_assert_double_le(fabs(context->random(context)), 1.);

        /* Test random dump & load */
#define TEST_PRNG_DUMP "random.state"
        FILE * stream = fopen(TEST_PRNG_DUMP, "w+");
        reset_error();
        pumas_context_random_dump(context, stream);
        fclose(stream);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        const double u2 = context->random(context);

        stream = fopen(TEST_PRNG_DUMP, "r");
        reset_error();
        pumas_context_random_load(context, stream);
        fclose(stream);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(context->random(context), u2);

        pumas_context_destroy(&context);
        pumas_context_create(&context, physics, 0);
        reset_error();
        pumas_context_random_seed_get(context, &seed);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        /* Free the data */
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);
}
END_TEST

/* Test the recorder API */
START_TEST(test_api_recorder)
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
        int i, *data;

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
START_TEST(test_api_print)
{
        /* Test the initialisation_error */
        reset_error();
        pumas_physics_print(physics, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_PHYSICS_ERROR);

        /* Load the tau data */
        load_tau();

        /* Check the NULL stream case */
        reset_error();
        pumas_physics_print(physics, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_IO_ERROR);

        /* Check printing to a file */
        FILE * stream = fopen("tau.json", "w+");
        pumas_physics_print(physics, stream, NULL, NULL);
        fputs("\n===\n", stream);
        pumas_physics_print(physics, stream, "    ", "\n");
        fclose(stream);

        /* Free the data */
        pumas_physics_destroy(&physics);
}
END_TEST


static double dummy_dcs(double Z, double A, double m, double K, double q)
{
        return 0;
}


/* Test the dcs API */
START_TEST(test_api_dcs)
{
        /* Test the process error */
        pumas_dcs_t * dcs = NULL;
        reset_error();
        pumas_dcs_get(-1, "KKP", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_dcs_register(-1, "dummy", &dummy_dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Test the model error  */
        reset_error();
        pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG, "toto", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MODEL_ERROR);

        reset_error();
        pumas_dcs_register(PUMAS_PROCESS_BREMSSTRAHLUNG, "KKP", &dummy_dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MODEL_ERROR);

        /* Test the NULL error */
        reset_error();
        pumas_dcs_register(PUMAS_PROCESS_BREMSSTRAHLUNG, NULL, &dummy_dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        pumas_dcs_register(PUMAS_PROCESS_BREMSSTRAHLUNG, "dummy", NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        /* Test the default name getter */
        const char * model;
        model = pumas_dcs_default(PUMAS_PROCESS_BREMSSTRAHLUNG);
        ck_assert_str_eq("KKP", model);
        model = pumas_dcs_default(PUMAS_PROCESS_PHOTONUCLEAR);
        ck_assert_str_eq("DRSS", model);
        model = pumas_dcs_default(PUMAS_PROCESS_PAIR_PRODUCTION);
        ck_assert_str_eq("KKP", model);
        model = pumas_dcs_default(-1);
        ck_assert_ptr_null(model);

        /* Test the default dcs getter */
        reset_error();
        pumas_dcs_t * dcs_br;
        pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG, NULL, &dcs_br);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(dcs_br);

        pumas_dcs_t * dcs_pp;
        pumas_dcs_get(PUMAS_PROCESS_PAIR_PRODUCTION, NULL, &dcs_pp);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_ne(dcs_br, dcs_pp);
        ck_assert_ptr_nonnull(dcs_pp);

        pumas_dcs_t * dcs_pn;
        pumas_dcs_get(PUMAS_PROCESS_PHOTONUCLEAR, NULL, &dcs_pn);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(dcs_pn);
        ck_assert_ptr_ne(dcs_br, dcs_pn);
        ck_assert_ptr_ne(dcs_pp, dcs_pn);

        pumas_dcs_t * dcs_el;
        pumas_dcs_get(PUMAS_PROCESS_ELASTIC, NULL, &dcs_el);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(dcs_el);
        ck_assert_ptr_ne(dcs_br, dcs_el);
        ck_assert_ptr_ne(dcs_pp, dcs_el);
        ck_assert_ptr_ne(dcs_pn, dcs_el);

        /* Test the getter */
        pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG, "KKP", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(dcs_br, dcs);

        pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG, "ABB", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(dcs);
        ck_assert_ptr_ne(dcs_br, dcs);

        pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG, "SSR", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(dcs);
        ck_assert_ptr_ne(dcs_br, dcs);

        pumas_dcs_get(PUMAS_PROCESS_PAIR_PRODUCTION, "KKP", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(dcs_pp, dcs);

        pumas_dcs_get(PUMAS_PROCESS_PAIR_PRODUCTION, "SSR", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_nonnull(dcs);
        ck_assert_ptr_ne(dcs_pp, dcs);

        pumas_dcs_get(PUMAS_PROCESS_PHOTONUCLEAR, "DRSS", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(dcs_pn, dcs);

        pumas_dcs_get(PUMAS_PROCESS_PHOTONUCLEAR, "BM", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_ne(dcs_pn, dcs);

        pumas_dcs_get(PUMAS_PROCESS_PHOTONUCLEAR, "BBKS", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_ne(dcs_pn, dcs);

        /* Test the setter */
        pumas_dcs_register(PUMAS_PROCESS_BREMSSTRAHLUNG, "dummy0", &dummy_dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG, "dummy0", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(&dummy_dcs, dcs);

        pumas_dcs_register(PUMAS_PROCESS_PAIR_PRODUCTION, "dummy1", &dummy_dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        pumas_dcs_get(PUMAS_PROCESS_PAIR_PRODUCTION, "dummy1", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(&dummy_dcs, dcs);

        pumas_dcs_register(PUMAS_PROCESS_PHOTONUCLEAR, "dummy2", &dummy_dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        pumas_dcs_get(PUMAS_PROCESS_PHOTONUCLEAR, "dummy2", &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(&dummy_dcs, dcs);

        /* Load the muon data */
        load_muon();

        /* Test the physics getter */
        reset_error();
        pumas_physics_dcs(physics, -1, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        reset_error();
        pumas_physics_dcs(
            physics, PUMAS_PROCESS_BREMSSTRAHLUNG, &model, &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(dcs_br, dcs);
        ck_assert_str_eq("KKP", model);

        pumas_physics_dcs(
            physics, PUMAS_PROCESS_PHOTONUCLEAR, &model, &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(dcs_pn, dcs);
        ck_assert_str_eq("DRSS", model);

        pumas_physics_dcs(
            physics, PUMAS_PROCESS_PAIR_PRODUCTION, &model, &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(dcs_pp, dcs);
        ck_assert_str_eq("KKP", model);

        pumas_physics_dcs(
            physics, PUMAS_PROCESS_BREMSSTRAHLUNG, NULL, &dcs);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_ptr_eq(dcs_br, dcs);

        pumas_physics_dcs(
            physics, PUMAS_PROCESS_PHOTONUCLEAR, &model, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_str_eq("DRSS", model);

        /* Test some numerical values by comparing to fig. 4 of
         * D.E. Groom, N.V. Mokhov, and S.I. Striganov,
         * Atomic Data and Nuclear Data Tables 78, Number 2 (July 2001)
         */
        const double Z = 26;
        const double A = 55.8452;
        const double m = 0.10566;
        const double NA = 6.02214076E+23;
        const double k = 1E+04 * NA / A; /* m^2 to cm^2 / g */

        double p = 10;
        double E = sqrt(p * p + m * m);
        double q = 0.1 * E;
        ck_assert_double_eq_tol(
            dcs_br(Z, A, m, E - m, q) * E * k, 2.5E-05, 5E-06);
        ck_assert_double_eq_tol(
            dcs_pp(Z, A, m, E - m, q) * E * k, 1E-05, 5E-06);
        ck_assert_double_eq_tol(
            dcs_pn(Z, A, m, E - m, q) * E * k, 1E-05, 5E-06);

        p = 10E+03;
        E = sqrt(p * p + m * m);
        q = 1E-03 * E;
        ck_assert_double_eq_tol(
            dcs_br(Z, A, m, E - m, q) * E * k, 5E-03, 3E-03);
        ck_assert_double_eq_tol(
            dcs_pp(Z, A, m, E - m, q) * E * k, 0.5, 0.2);
        ck_assert_double_eq_tol(
            dcs_pn(Z, A, m, E - m, q) * E * k, 2E-03, 5E-04);

        pumas_dcs_get(PUMAS_PROCESS_PHOTONUCLEAR, "BBKS", &dcs);
        E = 3E+02;
        ck_assert_double_nonnan(dcs(Z, A, m, E - m, 3E-07 * E));
        E = 1E+10;
        ck_assert_double_nonnan(dcs(Z, A, m, E - m, 3E-07 * E));
        ck_assert_double_nonnan(dcs(Z, A, m, E - m, 0.5 * E));

        /* XXX Check values of elastic DCS ? */

        /* Free the data */
        pumas_physics_destroy(&physics);

        /* Check the dcs setter */
        reset_error();
        struct pumas_physics_settings settings = {
                .bremsstrahlung = "ABB",
                .pair_production = "SSR",
        };
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", "materials/dedx/muon",
            &settings);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        pumas_physics_dcs(physics, PUMAS_PROCESS_BREMSSTRAHLUNG, &model, &dcs);
        ck_assert_str_eq(settings.bremsstrahlung, model);
        pumas_dcs_get(PUMAS_PROCESS_BREMSSTRAHLUNG, "ABB", &dcs_br);
        ck_assert_ptr_eq(dcs_br, dcs);

        pumas_physics_dcs(physics, PUMAS_PROCESS_PAIR_PRODUCTION, &model, &dcs);
        ck_assert_str_eq(settings.pair_production, model);
        pumas_dcs_get(PUMAS_PROCESS_PAIR_PRODUCTION, "SSR", &dcs_pp);
        ck_assert_ptr_eq(dcs_pp, dcs);

        pumas_physics_destroy(&physics);

        reset_error();
        settings.bremsstrahlung = "toto";
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", "materials/dedx/muon",
            &settings);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MODEL_ERROR);

        reset_error();
        settings.bremsstrahlung = NULL;
        settings.pair_production = "toto";
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", "materials/dedx/muon",
            &settings);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MODEL_ERROR);

        reset_error();
        settings.pair_production = NULL;
        settings.photonuclear = "toto";
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", "materials/dedx/muon",
            &settings);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MODEL_ERROR);
}
END_TEST

/* Geometry for test cases  */
static struct {
        int uniform;
        double magnet[3];
        double last_position[3];
        double position_rt[3];
        double direction_rt[3];
} geometry = { 1, { 0., 0., 0. }, { -1., 1., -1. }, { 0., 0., 0. },
        { 0., 0., 0. } };

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

        return lambda;
}

static enum pumas_step geometry_medium(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr,
    double * step_ptr)
{
        /* Cache the position for checking the locals API */
        memcpy(geometry.last_position, state->position,
            sizeof geometry.last_position);

        if ((medium_ptr == NULL) && (step_ptr == NULL)) {
                return PUMAS_STEP_CHECK;
        }

        if (step_ptr == NULL) {
                /* This is a pure location call. Let us check that the
                 * point is along the last ray tracing call
                 */
                const double r[3] = { state->position[0] -
                            geometry.position_rt[0],
                        state->position[1] - geometry.position_rt[1],
                        state->position[2] - geometry.position_rt[2] };
                const double * const u = geometry.direction_rt;
                double det, tmp;
                det = fabs(r[1] * u[2] - r[2] * u[1]);
                tmp = fabs(r[2] * u[0] - r[0] * u[2]);
                if (tmp > det) det = tmp;
                tmp = fabs(r[0] * u[1] - r[1] * u[0]);
                if (tmp > det) det = tmp;
                if (det > FLT_EPSILON) {
                        ck_assert_double_le(det, FLT_EPSILON);
                }
        } else {
                /* Cache the current point for validation of subsequent
                 * location calls.
                 */
                memcpy(geometry.position_rt, state->position,
                    sizeof geometry.position_rt);
                memcpy(geometry.direction_rt, state->direction,
                    sizeof geometry.direction_rt);
        }

        static struct pumas_medium media[2] = { { 0, &rock_locals },
                { 1, &air_locals } };

        if (geometry.uniform) {
                if (medium_ptr != NULL) *medium_ptr = media;
                if (step_ptr != NULL) *step_ptr = 0.;
                return PUMAS_STEP_CHECK;
        } else {
                if (medium_ptr != NULL) *medium_ptr = NULL;
                const double z = state->position[2];
                if (z < -0.5 * TEST_ROCK_DEPTH) {
                        if (step_ptr != NULL) *step_ptr = 0.;
                        return PUMAS_STEP_CHECK;
                } else if (z < 0.5 * TEST_ROCK_DEPTH) {
                        if (medium_ptr != NULL) *medium_ptr = media;

                        const double dz1 = z + 0.5 * TEST_ROCK_DEPTH;
                        const double dz2 = 0.5 * TEST_ROCK_DEPTH - z;
                        const double dz = dz1 < dz2 ? dz1 : dz2;
                        if (step_ptr != NULL) *step_ptr = dz > 0. ? dz : 1E-03;
                        return PUMAS_STEP_CHECK;
                } else if (z < TEST_MAX_ALTITUDE) {
                        if (medium_ptr != NULL) *medium_ptr = media + 1;

                        const double sgn =
                            (context->mode.direction == PUMAS_MODE_FORWARD) ?
                            1. : -1.;
                        const double uz = state->direction[2] * sgn;
                        if (fabs(uz) < 1E-03) {
                                if (step_ptr != NULL) *step_ptr = 1E+03;
                                return PUMAS_STEP_CHECK;
                        }

                        double s;
                        if (uz > 0.)
                                s = (TEST_MAX_ALTITUDE - z) / uz;
                        else
                                s = (0.5 * TEST_ROCK_DEPTH - z) / uz;
                        if (step_ptr != NULL) *step_ptr = s > 0. ? s : 1E-03;
                        return PUMAS_STEP_CHECK;
                } else {
                        if (step_ptr != NULL) *step_ptr = 0.;
                        return PUMAS_STEP_CHECK;
                }
        }
}

/* Fixtures for Lossless tests */
static void lossless_setup(void)
{
        /* Load the muon data and create a simulation context */
        load_muon();
        geometry.uniform = 1;
        pumas_context_create(&context, physics, 0);
        context->medium = &geometry_medium;
        context->mode.scattering = PUMAS_MODE_LONGITUDINAL;
        context->mode.energy_loss = PUMAS_MODE_VIRTUAL;
}

static void lossless_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);
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
START_TEST(test_lossless_straight)
{
        int i, forward;
        enum pumas_event event_data, *event;
        struct pumas_medium *media_data[2], **media;
        struct pumas_medium * rock;
        double tmp;
        geometry_medium(context, state, &rock, &tmp);

        double ctau, mu;
        pumas_physics_particle(physics, NULL, &ctau, &mu);

        /* Check some basic errors */
        reset_error();
        initialise_state();
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MISSING_LIMIT);

        for (forward = 0; forward < 2; forward++) {
                context->mode.direction =
                    forward ? PUMAS_MODE_FORWARD : PUMAS_MODE_BACKWARD;
                const double sgn = forward ? 1. : -1.;

                /* Check the forward distance limit */
                context->event = PUMAS_EVENT_LIMIT_DISTANCE;
                context->limit.distance = 1E+03;

                double k = 1.;
                double d = context->limit.distance;
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
                        state->energy = k;

                        pumas_context_transport(context, state, event, media);
                        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                        ck_assert_double_eq(state->charge, -1.);
                        ck_assert_double_eq(state->energy, k);
                        ck_assert_double_eq(state->distance, d);
                        ck_assert_double_eq(state->grammage, X);
                        ck_assert_double_eq_tol(state->time, t, FLT_EPSILON);
                        ck_assert_double_eq_tol(state->weight,
                            exp(-state->time / ctau), FLT_EPSILON);
                        ck_assert_double_eq(state->position[0], 0.);
                        ck_assert_double_eq(state->position[1], 0.);
                        ck_assert_double_eq_tol(state->position[2],
                            sgn * state->distance, FLT_EPSILON);
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
                context->limit.distance = 0.;

                /* Check the forward grammage limit */
                context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;
                context->limit.grammage = 1E+03;

                X = context->limit.grammage;
                d = X / TEST_ROCK_DENSITY;
                t = d / sqrt(gamma * gamma - 1.);
                for (i = 0; i < 2; i++) {
                        if (i == 0)
                                event = NULL, media = NULL;
                        else
                                event = &event_data, media = media_data;

                        reset_error();
                        initialise_state();
                        state->energy = k;

                        pumas_context_transport(context, state, event, media);
                        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                        ck_assert_double_eq(state->charge, -1.);
                        ck_assert_double_eq(state->energy, k);
                        ck_assert_double_eq(state->distance, d);
                        ck_assert_double_eq(state->grammage, X);
                        ck_assert_double_eq_tol(state->time, t, FLT_EPSILON);
                        ck_assert_double_eq_tol(state->weight,
                            exp(-state->time / ctau), FLT_EPSILON);
                        ck_assert_double_eq(state->position[0], 0.);
                        ck_assert_double_eq(state->position[1], 0.);
                        ck_assert_double_eq_tol(state->position[2],
                            sgn * state->distance, FLT_EPSILON);
                        ck_assert_double_eq(state->direction[0], 0.);
                        ck_assert_double_eq(state->direction[1], 0.);
                        ck_assert_double_eq(state->direction[2], 1.);
                        ck_assert_int_eq(state->decayed, 0);

                        if (i == 0) {
                                ck_assert_ptr_null(event);
                                ck_assert_ptr_null(media);
                        } else {
                                ck_assert_int_eq(
                                    *event, PUMAS_EVENT_LIMIT_GRAMMAGE);
                                ck_assert_ptr_eq(media[0], rock);
                                ck_assert_ptr_eq(media[0], media[1]);
                        }
                }

                /* Check the forward time limit */
                context->event = PUMAS_EVENT_LIMIT_TIME;
                context->limit.time = 1E+03;

                t = context->limit.time;
                d = t * sqrt(gamma * gamma - 1.);
                X = d * TEST_ROCK_DENSITY;
                for (i = 0; i < 2; i++) {
                        if (i == 0)
                                event = NULL, media = NULL;
                        else
                                event = &event_data, media = media_data;

                        reset_error();
                        initialise_state();
                        state->energy = k;

                        pumas_context_transport(context, state, event, media);
                        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                        ck_assert_double_eq(state->charge, -1.);
                        ck_assert_double_eq(state->energy, k);
                        ck_assert_double_eq(state->distance, d);
                        ck_assert_double_eq(state->grammage, X);
                        ck_assert_double_eq_tol(state->time, t, FLT_EPSILON);
                        ck_assert_double_eq_tol(state->weight,
                            exp(-state->time / ctau), FLT_EPSILON);
                        ck_assert_double_eq(state->position[0], 0.);
                        ck_assert_double_eq(state->position[1], 0.);
                        ck_assert_double_eq_tol(state->position[2],
                            sgn * state->distance, FLT_EPSILON);
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

START_TEST(test_lossless_magnet)
{
        geometry.magnet[1] = 0.1;
        context->limit.distance = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        double ctau, mu;
        pumas_physics_particle(physics, NULL, &ctau, &mu);

        /* Check the MC deflection and the analytical computation when no
         * energy loss
         */
        const double tol = 1E-02;
        const double k = 1E+00;
        const double gamma = k / mu + 1.;
        const double bg = sqrt(gamma * gamma - 1.);
        const double d = context->limit.distance;
        const double rL = mu * bg / (geometry.magnet[1] * 0.299792458);
        const double X = d * TEST_ROCK_DENSITY;
        const double t = d / bg;
        const double phi = d / rL;
        double *r = state->position, *u = state->direction;

        reset_error();
        initialise_state();
        state->energy = k;

        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->energy, k);
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
        context->limit.distance = 0.;
        context->event = PUMAS_EVENT_NONE;
}
END_TEST

START_TEST(test_lossless_geometry)
{
        int i;
        enum pumas_event event;
        struct pumas_medium *media[2], *rock, *air;
        struct pumas_frame * frame;
        pumas_recorder_create(&recorder, 0);
        context->recorder = recorder;
        recorder->period = 0;

        geometry.uniform = 0;
        initialise_state();
        double tmp;
        geometry_medium(context, state, &rock, &tmp);
        state->position[2] = 0.5 * (0.5 * TEST_ROCK_DEPTH + TEST_MAX_ALTITUDE);
        geometry_medium(context, state, &air, &tmp);

        /* Test the initially out of world case */
        state->position[2] = 2 * TEST_MAX_ALTITUDE;
        reset_error();
        pumas_context_transport(context, state, &event, media);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
        ck_assert_ptr_eq(media[0], NULL);
        ck_assert_ptr_eq(media[1], NULL);
        ck_assert_int_eq(recorder->length, 1);
        frame = recorder->first;
        ck_assert_ptr_eq(frame->medium, NULL);
        ck_assert_int_eq(frame->event,
            PUMAS_EVENT_START | PUMAS_EVENT_MEDIUM | PUMAS_EVENT_STOP);
        ck_assert_ptr_null(frame->next);

        for (i = 0; i < 2; i++) {
                /* Test the unconstrained transport */
                context->event = PUMAS_EVENT_NONE;
                initialise_state();
                state->energy = 1E+03;
                if (i) {
                        context->mode.direction = PUMAS_MODE_BACKWARD;
                        state->direction[2] = -1.;
                } else {
                        context->mode.direction = PUMAS_MODE_FORWARD;
                }

                reset_error();
                pumas_recorder_clear(recorder);
                pumas_context_transport(context, state, &event, media);
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
                ck_assert_int_eq(
                    frame->event, PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);

                /* Test the medium stop condition */
                context->event = PUMAS_EVENT_MEDIUM;
                initialise_state();
                state->energy = 1E+03;
                if (i) {
                        state->direction[2] = -1.;
                }

                reset_error();
                pumas_recorder_clear(recorder);
                pumas_context_transport(context, state, &event, media);
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
                ck_assert_int_eq(
                    frame->event, PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);

                pumas_recorder_clear(recorder);
                pumas_context_transport(context, state, &event, media);
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
                ck_assert_int_eq(
                    frame->event, PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);
        }

        pumas_recorder_destroy(&recorder);
        context->recorder = NULL;
        geometry.uniform = 1;
        context->event = PUMAS_EVENT_NONE;
        context->mode.direction = PUMAS_MODE_FORWARD;
}
END_TEST

/* Fixtures for CSDA tests */
static void csda_setup(void)
{
        /* Load the muon data and create a simulation context */
        load_muon();
        geometry.uniform = 1;
        pumas_context_create(&context, physics, 0);
        context->medium = &geometry_medium;
        context->mode.scattering = PUMAS_MODE_LONGITUDINAL;
        context->mode.energy_loss = PUMAS_MODE_CSDA;
}

static void csda_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);
}

/* Test the CSDA straight case */
START_TEST(test_csda_straight)
{
        int i;
        enum pumas_event event_data, *event;
        struct pumas_medium *media_data[2], **media;
        struct pumas_medium * rock;
        double tmp;
        geometry_medium(context, state, &rock, &tmp);

        double ctau;
        pumas_physics_particle(physics, NULL, &ctau, NULL);

        /* Check some basic errors */
        reset_error();
        pumas_context_transport(NULL, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        pumas_context_transport(context, NULL, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_VALUE_ERROR);

        reset_error();
        initialise_state();
        state->direction[0] = 10;
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_DIRECTION_ERROR);

        reset_error();
        initialise_state();
        state->decayed = 1;
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        event = &event_data;
        pumas_context_transport(context, state, event, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(*event, PUMAS_EVENT_VERTEX_DECAY);

        context->medium = NULL;
        reset_error();
        initialise_state();
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MEDIUM_ERROR);
        context->medium = &geometry_medium;

        pumas_random_cb * prng = context->random;
        context->random = NULL;
        context->mode.decay = PUMAS_MODE_DECAY;
        reset_error();
        initialise_state();
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MISSING_RANDOM);
        context->random = prng;
        context->mode.decay = PUMAS_MODE_WEIGHT;

        const double accuracy = context->accuracy;
        context->accuracy = 0;
        reset_error();
        initialise_state();
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_ACCURACY_ERROR);
        context->accuracy = 2;
        reset_error();
        initialise_state();
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_ACCURACY_ERROR);
        context->accuracy = accuracy;

        /* Check the forward full path transport */
        double X, d, t0;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_physics_property_grammage(
                    physics, PUMAS_MODE_CSDA, 0, state->energy, &X);
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_CSDA, 0, state->energy, &t0);
                t0 /= TEST_ROCK_DENSITY;
                d = X / TEST_ROCK_DENSITY;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->energy, 0.);
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
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
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

                context->limit.energy = 0.5;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.energy = 0.;
                context->limit.distance = 0.5 * d;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.distance = 0.;
                context->limit.grammage = 0.5 * X;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.grammage = 0.;
                context->limit.time = 0.5 * t0;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                context->limit.time = 0.;
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        /* Check the forward kinetic limit */
        double X1, d1, t1;
        context->limit.energy = 0.5;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_CSDA, 0, context->limit.energy, &X1);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_CSDA, 0, context->limit.energy, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, context->limit.energy);
                ck_assert_double_eq_tol(state->distance, d - d1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-(t0 - t1) / ctau), FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->position[2], d - d1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[2], 1., FLT_EPSILON);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->limit.energy = 0.;

        /* Check the forward grammage limit */
        double k1;
        X1 = 0.5 * X;
        d1 = X1 / TEST_ROCK_DENSITY;
        context->limit.grammage = X1;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_CSDA, 0, k1, &t1);
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
        context->limit.grammage = 0.;

        /* Check the forward distance limit */
        d1 = 0.5 * d;
        X1 = d1 * TEST_ROCK_DENSITY;
        context->limit.distance = d1;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_CSDA, 0, k1, &t1);
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
        context->limit.distance = 0.;

        /* Check the forward time limit */
        t1 = 0.5 * t0;
        context->limit.time = t1;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_grammage(
                    physics, PUMAS_MODE_CSDA, 0, k1, &X1);
                d1 = X1 / TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq_tol(state->distance, d - d1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->weight, exp(-(t0 - t1) / ctau), FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->position[2], d - d1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[2], 1., FLT_EPSILON);
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
        context->limit.time = 0.;

        /* Check the backward transport with a kinetic limit */
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.energy = 1.;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;
        k1 = 0.5;
        pumas_physics_property_grammage(physics, PUMAS_MODE_CSDA, 0, k1, &X1);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_CSDA, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        double de0, de1;
        pumas_physics_property_energy_loss(
            physics, PUMAS_MODE_CSDA, 0, context->limit.energy, &de0);
        pumas_physics_property_energy_loss(
            physics, PUMAS_MODE_CSDA, 0, k1, &de1);
        double w = exp(-(t0 - t1) / ctau) * de0 / de1;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->energy, context->limit.energy);
                ck_assert_double_eq_tol(state->distance, d - d1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->weight, w, FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->position[2], -(d - d1), FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[2], 1., FLT_EPSILON);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->limit.energy = 0.;

        /* Check the backward transport with a grammage limit */
        context->limit.grammage = X - X1;
        context->limit.energy = 0.75;          /* Not activated  */
        context->limit.distance = 0.5 * (d - d1); /* Not activated  */
        context->limit.time = 0.5 * (t0 - t1);    /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq_tol(state->energy, 1., FLT_EPSILON);
                ck_assert_double_eq_tol(state->distance, d - d1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->weight, w, FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->position[2], -(d - d1), FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[2], 1., FLT_EPSILON);
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
        context->limit.grammage = 0.;

        /* Check the backward transport with a distance limit */
        context->limit.distance = d - d1;
        context->limit.grammage = 0.5 * (X - X1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq_tol(state->energy, 1., FLT_EPSILON);
                ck_assert_double_eq_tol(state->distance, d - d1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->time, t0 - t1, FLT_EPSILON);
                ck_assert_double_eq_tol(state->weight, w, FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->position[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(
                    state->position[2], -(d - d1), FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[0], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[1], 0., FLT_EPSILON);
                ck_assert_double_eq_tol(state->direction[2], 1., FLT_EPSILON);
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
        context->limit.grammage = 0.;

        /* Check the backward transport with a time limit */
        context->limit.time = t0 - t1;
        context->limit.distance = 0.5 * (d - d1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq_tol(state->energy, 1., FLT_EPSILON);
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
        context->limit.grammage = 0.;
        context->limit.time = 0.5 * (t0 - t1); /* Not activated */

        /* Check the case of an initial kinetic limit violation */
        context->limit.energy = 0.5;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1.);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->mode.direction = PUMAS_MODE_FORWARD;
        context->limit.energy = 1;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 0.5;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.5);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.energy = 0.;

        /* Check the case of an initial grammage limit violation */
        context->limit.grammage = 0.5 * X;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                state->grammage = X;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1.);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                state->grammage = X;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1.);
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
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.grammage = 0.;

        /* Check the case of an initial distance limit violation */
        context->limit.distance = 0.5 * d;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                state->distance = d;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1.);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                state->distance = d;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1.);
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
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.distance = 0.;

        /* Check the case of an initial time limit violation */
        context->limit.time = 0.5 * t0;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                state->time = t0;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1.);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1.;
                state->time = t0;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1.);
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
        context->limit.time = 0.;
        context->event = PUMAS_EVENT_NONE;

        /* Check the overflow case */
        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+08;
                pumas_physics_property_grammage(
                    physics, PUMAS_MODE_CSDA, 0, state->energy, &X);
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_CSDA, 0, state->energy, &t0);
                t0 /= TEST_ROCK_DENSITY;
                d = X / TEST_ROCK_DENSITY;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->energy, 0.);
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
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
}
END_TEST

START_TEST(test_csda_record)
{
        double ctau;
        pumas_physics_particle(physics, NULL, &ctau, NULL);

        struct pumas_medium * rock;
        double tmp;
        geometry_medium(context, state, &rock, &tmp);

        pumas_recorder_create(&recorder, 0);
        context->recorder = recorder;

        reset_error();
        initialise_state();
        state->energy = 1.;

        double X, d, t0;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_CSDA, 0, state->energy, &X);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_CSDA, 0, state->energy, &t0);
        t0 /= TEST_ROCK_DENSITY;
        d = X / TEST_ROCK_DENSITY;
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_int_eq(recorder->length, 2);
        ck_assert_ptr_nonnull(recorder->first);

        struct pumas_frame * frame = recorder->first;
        ck_assert_ptr_eq(frame->medium, rock);
        ck_assert_double_eq(frame->state.charge, -1.);
        ck_assert_double_eq(frame->state.energy, 1.);
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
        ck_assert_double_eq(frame->state.energy, 0.);
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
        ck_assert_int_eq(
            frame->event, PUMAS_EVENT_LIMIT_ENERGY | PUMAS_EVENT_STOP);
        ck_assert_ptr_null(frame->next);

        pumas_recorder_destroy(&recorder);
        context->recorder = NULL;
}
END_TEST

START_TEST(test_csda_magnet)
{
        geometry.magnet[1] = 0.1;

        reset_error();
        initialise_state();
        state->energy = 1.;

        double X, d, phi0, phi1, phi, *u = state->direction;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_CSDA, 0, state->energy, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_physics_property_magnetic_rotation(
            physics, 0, state->energy, &phi0);
        pumas_physics_property_magnetic_rotation(physics, 0, 0., &phi1);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_eq(u[0], sin(phi));
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_eq(u[2], cos(phi));

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        geometry.magnet[1] = 0.;
}
END_TEST

START_TEST(test_csda_geometry)
{
        int i;
        enum pumas_event event;
        struct pumas_medium *media[2], *rock, *air;
        struct pumas_frame * frame;
        pumas_recorder_create(&recorder, 0);
        context->recorder = recorder;
        recorder->period = 0;

        geometry.uniform = 0;
        initialise_state();
        double tmp;
        geometry_medium(context, state, &rock, &tmp);
        state->position[2] = 0.5 * (0.5 * TEST_ROCK_DEPTH + TEST_MAX_ALTITUDE);
        geometry_medium(context, state, &air, &tmp);

        /* Test the initially out of world case */
        state->position[2] = 2 * TEST_MAX_ALTITUDE;
        reset_error();
        pumas_context_transport(context, state, &event, media);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(event, PUMAS_EVENT_MEDIUM);
        ck_assert_ptr_eq(media[0], NULL);
        ck_assert_ptr_eq(media[1], NULL);
        ck_assert_int_eq(recorder->length, 1);
        frame = recorder->first;
        ck_assert_ptr_eq(frame->medium, NULL);
        ck_assert_int_eq(frame->event,
            PUMAS_EVENT_START | PUMAS_EVENT_MEDIUM | PUMAS_EVENT_STOP);
        ck_assert_ptr_null(frame->next);

        for (i = 0; i < 2; i++) {
                /* Test the unconstrained transport */
                context->event = PUMAS_EVENT_NONE;
                initialise_state();
                state->energy = 1E+03;
                if (i) {
                        context->mode.direction = PUMAS_MODE_BACKWARD;
                        state->direction[2] = -1.;
                } else {
                        context->mode.direction = PUMAS_MODE_FORWARD;
                }

                reset_error();
                pumas_recorder_clear(recorder);
                pumas_context_transport(context, state, &event, media);
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
                ck_assert_int_eq(
                    frame->event, PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);

                /* Test the medium stop condition */
                context->event = PUMAS_EVENT_MEDIUM;
                initialise_state();
                state->energy = 1E+03;
                if (i) {
                        state->direction[2] = -1.;
                }

                reset_error();
                pumas_recorder_clear(recorder);
                pumas_context_transport(context, state, &event, media);
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
                ck_assert_int_eq(
                    frame->event, PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);

                pumas_recorder_clear(recorder);
                pumas_context_transport(context, state, &event, media);
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
                ck_assert_int_eq(
                    frame->event, PUMAS_EVENT_STOP | PUMAS_EVENT_MEDIUM);
        }

        pumas_recorder_destroy(&recorder);
        context->recorder = NULL;
        geometry.uniform = 1;
        context->event = PUMAS_EVENT_NONE;
        context->mode.direction = PUMAS_MODE_FORWARD;
}
END_TEST

/* Fixtures for hybrid tests */
static void hybrid_setup(void)
{
        /* Load the tau data and create a simulation context */
        load_muon();
        geometry.uniform = 1;
        pumas_context_create(&context, physics, 0);
        context->medium = &geometry_medium;
        context->mode.scattering = PUMAS_MODE_LONGITUDINAL;
        context->mode.energy_loss = PUMAS_MODE_HYBRID;
        unsigned long seed = 0;
        pumas_context_random_seed_set(context, &seed);
}

static void hybrid_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);
}

START_TEST(test_hybrid_straight)
{
        int i, ik;
        enum pumas_event event_data, *event;
        struct pumas_medium *media_data[2], **media;
        struct pumas_medium * rock;
        double tmp;
        geometry_medium(context, state, &rock, &tmp);

        double ctau;
        pumas_physics_particle(physics, NULL, &ctau, NULL);

        /* Check the missing random engine case */
        reset_error();
        pumas_random_cb * prng = context->random;
        context->random = NULL;
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MISSING_RANDOM);
        context->random = prng;

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
                        state->energy = (ik == 0) ? 1E+08 : 1E+03;
                        pumas_physics_property_grammage(physics,
                            PUMAS_MODE_HYBRID, 0, state->energy, &X);
                        pumas_physics_property_proper_time(physics,
                            PUMAS_MODE_HYBRID, 0, state->energy, &t0);
                        t0 /= TEST_ROCK_DENSITY;
                        d = X / TEST_ROCK_DENSITY;
                        pumas_context_transport(context, state, event, media);
                        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                        ck_assert_double_eq(state->charge, -1.);
                        ck_assert_double_eq(state->energy, 0.);
                        ck_assert_double_le(state->distance, d);
                        ck_assert_double_le(state->grammage, X);
                        ck_assert_double_le(state->time, t0);
                        ck_assert_double_le(state->weight,
                            exp(-state->time / ctau));
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
                                    *event, PUMAS_EVENT_LIMIT_ENERGY);
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

                context->limit.energy = 0.5E+03;
                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.energy = 0.;
                context->limit.distance = 0.5 * d;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.distance = 0.;
                context->limit.grammage = 0.5 * X;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.grammage = 0.;
                context->limit.time = 0.5 * t0;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                context->limit.time = 0.;
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        /* Check the forward kinetic limit */
        double X1, d1, t1;
        context->limit.energy = 0.5E+03;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, context->limit.energy, &X1);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, context->limit.energy, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->energy, context->limit.energy);
                ck_assert_double_le(state->distance, d - d1 + FLT_EPSILON);
                ck_assert_double_le(state->grammage, X - X1 + FLT_EPSILON);
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
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->limit.energy = 0.;

        /* Check the forward grammage limit */
        double k1;
        X1 = 1E-02 * X;
        d1 = X1 / TEST_ROCK_DENSITY;
        context->limit.grammage = X1;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_HYBRID, 0, k1, &t1);
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
        context->limit.grammage = 0.;

        /* Check the forward distance limit */
        d1 = 1E-02 * d;
        X1 = d1 * TEST_ROCK_DENSITY;
        context->limit.distance = d1;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_HYBRID, 0, k1, &t1);
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
        context->limit.distance = 0.;

        /* Check the forward time limit */
        t1 = 1E-02 * t0;
        context->limit.time = t1;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_grammage(
                    physics, PUMAS_MODE_HYBRID, 0, k1, &X1);
                d1 = X1 / TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->distance, d - d1 + FLT_EPSILON);
                ck_assert_double_le(state->grammage, X - X1 + FLT_EPSILON);
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
        context->limit.time = 0.;

        /* Check the backward transport with a kinetic limit */
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.energy = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;
        k1 = 0.5 * context->limit.energy;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, k1, &X1);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_ge(state->energy, context->limit.energy);
                ck_assert_double_le(state->distance, d - d1 + FLT_EPSILON);
                ck_assert_double_le(state->grammage, X - X1 + FLT_EPSILON);
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
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->limit.energy = 0.;

        /* Check the backward transport with a grammage limit */
        context->limit.grammage = X - X1;
        context->limit.energy = 0.75E+03;      /* Not activated  */
        context->limit.distance = 0.5 * (d - d1); /* Not activated  */
        context->limit.time = 0.5 * (t0 - t1);    /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->energy, k1);
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
        context->limit.grammage = 0.;

        /* Check the backward transport with a distance limit */
        context->limit.distance = d - d1;
        context->limit.grammage = 0.5 * (X - X1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->energy, k1);
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
        context->limit.grammage = 0.;

        /* Check the backward transport with a time limit */
        context->limit.time = t0 - t1;
        context->limit.distance = 0.5 * (d - d1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->energy, k1);
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
        context->limit.grammage = 0.;
        context->limit.time = 0.5 * (t0 - t1); /* Not activated */

        /* Check the case of an initial kinetic limit violation */
        context->limit.energy = 0.5E+03;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->mode.direction = PUMAS_MODE_FORWARD;
        context->limit.energy = 1E+03;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 0.5E+03;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.5E+03);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.energy = 0.;

        /* Check the case of an initial grammage limit violation */
        context->limit.grammage = 0.5 * X;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->grammage = X;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->grammage = X;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.grammage = 0.;

        /* Check the case of an initial distance limit violation */
        context->limit.distance = 0.5 * d;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->distance = d;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->distance = d;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.distance = 0.;

        /* Check the case of an initial time limit violation */
        context->limit.time = 0.5 * t0;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->time = t0;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E03;
                state->time = t0;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        context->limit.time = 0.;
        context->event = PUMAS_EVENT_NONE;
}
END_TEST

START_TEST(test_hybrid_scattering)
{
        context->mode.scattering = PUMAS_MODE_FULL_SPACE;

        double ctau;
        pumas_physics_particle(physics, NULL, &ctau, NULL);

        /* Check the forward transport with scattering */
        double X, d, t0, *u = state->direction;

        reset_error();
        initialise_state();
        state->energy = 1E+03;

        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &X);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &t0);
        t0 /= TEST_ROCK_DENSITY;
        d = X / TEST_ROCK_DENSITY;

        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->charge, -1.);
        ck_assert_double_eq(state->energy, 0.);
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
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.energy = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;

        double k1, X1, d1, t1;
        k1 = 0.5 * context->limit.energy;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, k1, &X1);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        reset_error();
        initialise_state();
        state->energy = k1;
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->charge, -1.);
        ck_assert_double_ge(state->energy, context->limit.energy);
        ck_assert_double_le(state->distance, d - d1 + FLT_EPSILON);
        ck_assert_double_le(state->grammage, X - X1 + FLT_EPSILON);
        ck_assert_double_le(state->time, t0 - t1 + FLT_EPSILON);
        ck_assert_double_ge(state->position[2], -(d - d1) - FLT_EPSILON);
        ck_assert_int_eq(state->decayed, 0);

        norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Restore the context status */
        context->event = PUMAS_EVENT_NONE;
        context->limit.energy = 0.;
        context->mode.direction = PUMAS_MODE_FORWARD;
        context->mode.scattering = PUMAS_MODE_LONGITUDINAL;
}
END_TEST

START_TEST(test_hybrid_record)
{
        double ctau;
        pumas_physics_particle(physics, NULL, &ctau, NULL);

        struct pumas_medium * rock;
        double tmp;
        geometry_medium(context, state, &rock, &tmp);

        pumas_recorder_create(&recorder, 0);
        context->recorder = recorder;

        reset_error();
        initialise_state();
        state->energy = 1E+03;

        double X, d, t0;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &X);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &t0);
        t0 /= TEST_ROCK_DENSITY;
        d = X / TEST_ROCK_DENSITY;
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_int_ge(recorder->length, 2);
        ck_assert_ptr_nonnull(recorder->first);

        struct pumas_frame * frame = recorder->first;
        ck_assert_ptr_eq(frame->medium, rock);
        ck_assert_double_eq(frame->state.charge, -1.);
        ck_assert_double_eq(frame->state.energy, 1E+03);
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
        ck_assert_double_eq(frame->state.energy, 0.);
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
        ck_assert_int_eq(
            frame->event, PUMAS_EVENT_LIMIT_ENERGY | PUMAS_EVENT_STOP);
        ck_assert_ptr_null(frame->next);

        pumas_recorder_destroy(&recorder);
        context->recorder = NULL;
}
END_TEST

START_TEST(test_hybrid_magnet)
{
        geometry.magnet[1] = 0.1;

        /* Compare the MC deflection and the analytical computation at low
         * energy, i.e. when no DEL */
        double tol = 5E-02;
        reset_error();
        initialise_state();
        state->energy = 1E+00;

        double X, d, phi0, phi1, phi, *u = state->direction;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_physics_property_magnetic_rotation(
            physics, 0, state->energy, &phi0);
        pumas_physics_property_magnetic_rotation(physics, 0, 0., &phi1);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_eq_tol(u[0], sin(phi), tol);
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_eq_tol(u[2], cos(phi), tol);

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Check a high energy case */
        reset_error();
        initialise_state();
        state->energy = 1E+03;

        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_physics_property_magnetic_rotation(
            physics, 0, state->energy, &phi0);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_context_transport(context, state, NULL, NULL);
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
        pumas_context_create(&context, physics, 0);
        context->medium = &geometry_medium;
        context->mode.scattering = PUMAS_MODE_FULL_SPACE;
        context->mode.energy_loss = PUMAS_MODE_DETAILED;
        unsigned long seed = 0;
        pumas_context_random_seed_set(context, &seed);
}

static void detailed_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);
}

START_TEST(test_detailed_straight)
{
        context->mode.scattering = PUMAS_MODE_LONGITUDINAL;

        int i;
        enum pumas_event event_data, *event;
        struct pumas_medium *media_data[2], **media;
        struct pumas_medium * rock;
        double tmp;
        geometry_medium(context, state, &rock, &tmp);

        double ctau;
        pumas_physics_particle(physics, NULL, &ctau, NULL);

        /* Check the missing random engine case */
        reset_error();
        pumas_random_cb * prng = context->random;
        context->random = NULL;
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_MISSING_RANDOM);
        context->random = prng;

        /* Check the forward full path transport */
        double X, d, t0;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_physics_property_grammage(
                    physics, PUMAS_MODE_HYBRID, 0, state->energy, &X);
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_HYBRID, 0, state->energy, &t0);
                t0 /= TEST_ROCK_DENSITY;
                d = X / TEST_ROCK_DENSITY;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_eq(state->energy, 0.);
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
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
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

                context->limit.energy = 0.5E+03;
                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.energy = 0.;
                context->limit.distance = 0.5 * d;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.distance = 0.;
                context->limit.grammage = 0.5 * X;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }

                context->limit.grammage = 0.;
                context->limit.time = 0.5 * t0;
                reset_error();
                initialise_state();
                state->energy = 1.;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.);
                context->limit.time = 0.;
                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        /* Check the forward kinetic limit */
        double X1, d1, t1;
        context->limit.energy = 0.5E+03;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, context->limit.energy, &X1);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, context->limit.energy, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->energy, context->limit.energy);
                ck_assert_double_le(state->distance, (d - d1) * 1.25);
                ck_assert_double_le(state->grammage, (X - X1) * 1.25);
                ck_assert_double_le(state->time, (t0 - t1) * 1.25);
                ck_assert_double_eq_tol(
                    state->weight, exp(-state->time / ctau), FLT_EPSILON);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_le(state->position[2], (d - d1) * 1.25);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->limit.energy = 0.;

        /* Check the forward grammage limit */
        double k1;
        X1 = 1E-02 * X;
        d1 = X1 / TEST_ROCK_DENSITY;
        context->limit.grammage = X1;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_HYBRID, 0, k1, &t1);
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
        context->limit.grammage = 0.;

        /* Check the forward distance limit */
        d1 = 1E-02 * d;
        X1 = d1 * TEST_ROCK_DENSITY;
        context->limit.distance = d1;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_proper_time(
                    physics, PUMAS_MODE_HYBRID, 0, k1, &t1);
                t1 /= TEST_ROCK_DENSITY;
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_le(state->distance, d - d1);
                ck_assert_double_le(state->grammage, X - X1);
                ck_assert_double_eq_tol(state->time / (t0 - t1), 1, 0.25);
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
        context->limit.distance = 0.;

        /* Check the forward time limit */
        t1 = 1E-02 * t0;
        context->limit.time = t1;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                k1 = state->energy;
                pumas_physics_property_grammage(
                    physics, PUMAS_MODE_HYBRID, 0, k1, &X1);
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
        context->limit.time = 0.;

        /* Check the backward transport with a kinetic limit */
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.energy = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;
        k1 = 0.5 * context->limit.energy;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, k1, &X1);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_ge(state->energy, context->limit.energy);
                ck_assert_double_le(state->distance, (d - d1) * 1.1);
                ck_assert_double_le(state->grammage, (X - X1) * 1.1);
                ck_assert_double_eq(state->position[0], 0.);
                ck_assert_double_eq(state->position[1], 0.);
                ck_assert_double_ge(state->position[2], -(d - d1) * 1.1);
                ck_assert_double_eq(state->direction[0], 0.);
                ck_assert_double_eq(state->direction[1], 0.);
                ck_assert_double_eq(state->direction[2], 1.);
                ck_assert_int_eq(state->decayed, 0);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->limit.energy = 0.;

        /* Check the backward transport with a grammage limit */
        context->limit.grammage = X - X1;
        context->limit.energy = 0.75E+03;      /* Not activated  */
        context->limit.distance = 0.5 * (d - d1); /* Not activated  */
        context->limit.time = 0.5 * (t0 - t1);    /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->energy, k1);
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
        context->limit.grammage = 0.;

        /* Check the backward transport with a distance limit */
        context->limit.distance = d - d1;
        context->limit.grammage = 0.5 * (X - X1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->energy, k1);
                ck_assert_double_eq(state->distance, d - d1);
                ck_assert_double_eq_tol(state->grammage, X - X1, FLT_EPSILON);
                ck_assert_double_le(state->time, (t0 - t1) * 1.1);
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
        context->limit.grammage = 0.;

        /* Check the backward transport with a time limit */
        context->limit.time = t0 - t1;
        context->limit.distance = 0.5 * (d - d1); /* Not activated  */
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = k1;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->charge, -1.);
                ck_assert_double_gt(state->energy, k1);
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
        context->limit.grammage = 0.;
        context->limit.time = 0.5 * (t0 - t1); /* Not activated */

        /* Check the case of an initial kinetic limit violation */
        context->limit.energy = 0.5E+03;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }

        context->mode.direction = PUMAS_MODE_FORWARD;
        context->limit.energy = 1E+03;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 0.5E+03;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 0.5E+03);

                if (i == 0) {
                        ck_assert_ptr_null(event);
                        ck_assert_ptr_null(media);
                } else {
                        ck_assert_int_eq(*event, PUMAS_EVENT_LIMIT_ENERGY);
                        ck_assert_ptr_eq(media[0], rock);
                        ck_assert_ptr_eq(media[0], media[1]);
                }
        }
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.energy = 0.;

        /* Check the case of an initial grammage limit violation */
        context->limit.grammage = 0.5 * X;
        context->event = PUMAS_EVENT_LIMIT_GRAMMAGE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->grammage = X;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->grammage = X;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.grammage = 0.;

        /* Check the case of an initial distance limit violation */
        context->limit.distance = 0.5 * d;
        context->event = PUMAS_EVENT_LIMIT_DISTANCE;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->distance = d;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->distance = d;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.distance = 0.;

        /* Check the case of an initial time limit violation */
        context->limit.time = 0.5 * t0;
        context->event = PUMAS_EVENT_LIMIT_TIME;

        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E+03;
                state->time = t0;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        for (i = 0; i < 2; i++) {
                if (i == 0)
                        event = NULL, media = NULL;
                else
                        event = &event_data, media = media_data;

                reset_error();
                initialise_state();
                state->energy = 1E03;
                state->time = t0;
                pumas_context_transport(context, state, event, media);
                ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
                ck_assert_double_eq(state->energy, 1E+03);
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

        context->mode.direction = PUMAS_MODE_FORWARD;
        context->limit.time = 0.;
        context->event = PUMAS_EVENT_NONE;
        context->mode.scattering = PUMAS_MODE_FULL_SPACE;
}
END_TEST

START_TEST(test_detailed_scattering)
{
        double ctau;
        pumas_physics_particle(physics, NULL, &ctau, NULL);

        /* Check the forward transport with scattering */
        double X, d, t0, *u = state->direction;

        reset_error();
        initialise_state();
        state->energy = 1E+03;

        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &X);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &t0);
        t0 /= TEST_ROCK_DENSITY;
        d = X / TEST_ROCK_DENSITY;

        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->charge, -1.);
        ck_assert_double_eq(state->energy, 0.);
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
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->limit.energy = 1E+03;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;

        double k1, X1, d1, t1;
        k1 = 1E-03;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, k1, &X1);
        pumas_physics_property_proper_time(
            physics, PUMAS_MODE_HYBRID, 0, k1, &t1);
        t1 /= TEST_ROCK_DENSITY;
        d1 = X1 / TEST_ROCK_DENSITY;

        reset_error();
        initialise_state();
        state->energy = k1;
        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(state->charge, -1.);
        ck_assert_double_ge(state->energy, context->limit.energy);
        ck_assert_double_eq_tol(state->distance / (d - d1), 1., 0.5);
        ck_assert_double_eq_tol(state->grammage / (X - X1), 1., 0.5);
        ck_assert_double_eq_tol(state->time / (t0 - t1), 1., 0.5);
        ck_assert_int_eq(state->decayed, 0);

        norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Restore the context status */
        context->event = PUMAS_EVENT_NONE;
        context->limit.energy = 0.;
        context->mode.direction = PUMAS_MODE_FORWARD;
}
END_TEST

START_TEST(test_detailed_magnet)
{
        context->mode.scattering = PUMAS_MODE_LONGITUDINAL;
        geometry.magnet[1] = 0.1;

        /* Compare the MC deflection and the analytical computation at low
         * energy, i.e. when no DEL */
        context->mode.direction = PUMAS_MODE_BACKWARD;
        context->event = PUMAS_EVENT_LIMIT_ENERGY;
        context->limit.energy = 1E+01;

        double tol = 5E-02;
        reset_error();
        initialise_state();
        state->energy = 1E+00;

        double X, d, phi0, phi1, phi, *u = state->direction;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_physics_property_magnetic_rotation(
            physics, 0, state->energy, &phi0);
        pumas_physics_property_magnetic_rotation(
            physics, 0, context->limit.energy, &phi1);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_eq_tol(u[0], sin(phi), tol);
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_eq_tol(u[2], cos(phi), tol);

        double norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Check a high energy case */
        context->mode.direction = PUMAS_MODE_FORWARD;
        context->event = PUMAS_EVENT_NONE;
        context->limit.energy = 0.;

        reset_error();
        initialise_state();
        state->energy = 1E+03;

        pumas_physics_property_grammage(
            physics, PUMAS_MODE_HYBRID, 0, state->energy, &X);
        d = X / TEST_ROCK_DENSITY;
        pumas_physics_property_magnetic_rotation(
            physics, 0, state->energy, &phi0);
        pumas_physics_property_magnetic_rotation(physics, 0, 0., &phi1);
        phi = -(phi1 - phi0) * geometry.magnet[1] / TEST_ROCK_DENSITY *
            state->charge;

        pumas_context_transport(context, state, NULL, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);

        ck_assert_double_le(u[0], sin(phi) + tol);
        ck_assert_double_eq(u[1], 0.);
        ck_assert_double_ge(u[2], cos(phi) - tol);

        norm2 = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        ck_assert_double_eq_tol(norm2, 1., FLT_EPSILON);

        /* Restore the context and the geometry */
        geometry.magnet[1] = 0.;
        context->mode.scattering = PUMAS_MODE_FULL_SPACE;
}
END_TEST

/* Fixtures for tau tests */
static void tau_setup(void)
{
        /* Load the tau data and create a simulation context */
        load_tau();
        geometry.uniform = 1;
        pumas_context_create(&context, physics, 0);
        context->medium = &geometry_medium;
        unsigned long seed = 0;
        pumas_context_random_seed_set(context, &seed);
}

static void tau_teardown(void)
{
        /* Destroy the simulation context and unload the data */
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);
}

START_TEST(test_tau_csda)
{
        enum pumas_event event;
        context->mode.scattering = PUMAS_MODE_LONGITUDINAL;
        context->mode.energy_loss = PUMAS_MODE_CSDA;

        /* Test some errors */
        context->mode.decay = PUMAS_MODE_WEIGHT;
        reset_error();
        initialise_state();
        pumas_context_transport(context, state, &event, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_DECAY_ERROR);
        context->mode.decay = PUMAS_MODE_DECAY;

        /* Test the forward transport */
        reset_error();
        initialise_state();
        state->energy = 1.;

        double X, d;
        pumas_physics_property_grammage(
            physics, PUMAS_MODE_CSDA, 0, state->energy, &X);
        d = X / TEST_ROCK_DENSITY;

        pumas_context_transport(context, state, &event, NULL);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(event, PUMAS_EVENT_VERTEX_DECAY);

        ck_assert_double_eq(state->direction[0], 0.);
        ck_assert_double_eq(state->direction[1], 0.);
        ck_assert_double_eq(state->direction[2], 1.);

        context->mode.scattering = PUMAS_MODE_FULL_SPACE;
        context->mode.energy_loss = PUMAS_MODE_DETAILED;
}
END_TEST

START_TEST(test_tabulation)
{
        double kinetic[2] = { 1E-01, 1E+01 };
        struct pumas_physics_settings settings = {
            .n_energies = 2, .energy = kinetic, .dry = 1, .update = 1};
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON,
            "materials/materials.xml", ".", &settings);

        FILE * stream;
        stream = fopen("./air.txt", "r");
        ck_assert_ptr_nonnull(stream);
        fclose(stream);

        stream = fopen("./standard_rock.txt", "r");
        ck_assert_ptr_nonnull(stream);
        fclose(stream);

        stream = fopen("./water.txt", "r");
        ck_assert_ptr_nonnull(stream);
        fclose(stream);
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
        tcase_add_test(tc_api, test_api_version);
        tcase_add_test(tc_api, test_api_constant);
        tcase_add_test(tc_api, test_api_init);
        tcase_add_test(tc_api, test_api_memory);
        tcase_add_test(tc_api, test_api_element);
        tcase_add_test(tc_api, test_api_material);
        tcase_add_test(tc_api, test_api_composite);
        tcase_add_test(tc_api, test_api_property);
        tcase_add_test(tc_api, test_api_table);
        tcase_add_test(tc_api, test_api_context);
        tcase_add_test(tc_api, test_api_random);
        tcase_add_test(tc_api, test_api_recorder);
        tcase_add_test(tc_api, test_api_print);
        tcase_add_test(tc_api, test_api_dcs);

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

        /* The tabulation test case */
        TCase * tc_tabulation = tcase_create("Tabulation");
        suite_add_tcase(suite, tc_tabulation);
        tcase_set_timeout(tc_tabulation, timeout);
        tcase_add_test(tc_tabulation, test_tabulation);

        return suite;
}

int main(void)
{
        /* Initialise the PRNG */
        srand(0);

        /* Override the default error handler */
        pumas_error_handler_set(&handle_error);

        /* Configure the tests and the runner */
        Suite * suite = create_suite();
        SRunner * runner = srunner_create(suite);
        srunner_set_fork_status(runner, CK_NOFORK);

        /* Run the tests */
        srunner_run_all(runner, CK_NORMAL);
        const int status =
            srunner_ntests_failed(runner) ? EXIT_FAILURE : EXIT_SUCCESS;
        srunner_free(runner);

        /* Return the test status to the OS */
        exit(status);
}
