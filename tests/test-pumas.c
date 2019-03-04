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


/* Test the error API */
START_TEST (test_error)
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
            (pumas_function_t *)&test_error));
}
END_TEST


/* Test the version tag */
START_TEST (test_tag)
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
START_TEST (test_init)
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
}
END_TEST


/* Test the materials API */
START_TEST (test_material)
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
START_TEST (test_composite)
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

        double density0 = 2.65E+03, density1 = 1.00E+03;
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


START_TEST (test_property)
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

        pumas_property_cross_section(0, 1E+06, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_lt(1. / (value * 2.65E+03), 1E+03);

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

        /* Unload the data */
        pumas_finalise();
}
END_TEST


START_TEST (test_table)
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

        pumas_table_value(PUMAS_PROPERTY_CROSS_SECTION, PUMAS_SCHEME_HYBRID,
            0, 0, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 0.);

        pumas_table_value(PUMAS_PROPERTY_CROSS_SECTION, PUMAS_SCHEME_HYBRID,
            0, 145, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_gt(value, 0.);
        ck_assert_double_lt(1. / (value * 2.65E+03), 1E+03);

        pumas_table_index(PUMAS_PROPERTY_CROSS_SECTION, PUMAS_SCHEME_HYBRID,
            0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 145);

        pumas_table_index(PUMAS_PROPERTY_CROSS_SECTION, PUMAS_SCHEME_HYBRID,
            0, 0., &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 0);

        pumas_table_index(PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_CSDA, 0,
            value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 0);

        pumas_table_value(PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_CSDA, 0,
            49, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 1.808E-04);

        pumas_table_index(PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_CSDA, 0,
            value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 49);

        pumas_table_value(PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_CSDA, 0,
            97, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_eq(value, 6.604E-04);

        pumas_table_index(PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_CSDA, 0,
            value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 97);

        pumas_table_value(PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_HYBRID, 0,
            97, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_double_lt(value, 6.604E-04);

        pumas_table_index(PUMAS_PROPERTY_ENERGY_LOSS, PUMAS_SCHEME_HYBRID, 0,
            value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_SUCCESS);
        ck_assert_int_eq(index, 97);

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

        /* Check the case of non tabulated properties  */
        pumas_table_value(
            PUMAS_PROPERTY_SCATTERING_LENGTH, 0, 0, 17, &value);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        pumas_table_index(PUMAS_PROPERTY_SCATTERING_LENGTH, 0,
            0, value, &index);
        ck_assert_int_eq(error_data.rc, PUMAS_RETURN_INDEX_ERROR);

        /* Unload the data */
        pumas_finalise();
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
        tcase_add_test(tc_api, test_error);
        tcase_add_test(tc_api, test_tag);
        tcase_add_test(tc_api, test_init);
        tcase_add_test(tc_api, test_material);
        tcase_add_test(tc_api, test_composite);
        tcase_add_test(tc_api, test_property);
        tcase_add_test(tc_api, test_table);

        return suite;
}


int main(void)
{
        /* Configure the tests amd the runner */
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
