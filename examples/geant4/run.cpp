/*
 * This is free and unencumbered software released into the public domain.
 *
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, for any purpose, commercial or non-commercial, and by any
 * means.
 *
 * In jurisdictions that recognize copyright laws, the author or authors
 * of this software dedicate any and all copyright interest in the
 * software to the public domain. We make this dedication for the benefit
 * of the public at large and to the detriment of our heirs and
 * successors. We intend this dedication to be an overt act of
 * relinquishment in perpetuity of all present and future rights to this
 * software under copyright law.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 *
 * For more information, please refer to <http://unlicense.org>
 */

/* This example illustrates how PUMAS can be used together with Geant4. The
 * geometry is loaded from a GDML file generated with the generate.cpp program.
 * A backward muon is transported through the geometry using a G4Navigator and
 * Monte Carlo steps are printed out.
 *
 * Note that the material names and properties defined in the GDML file must
 * be consistent with PUMAS ones (loaded from the physics dump).
 */

/* PUMAS API */
#include "pumas.h"
/* Geant4 wrappers for PUMAS */
#include "g4pumas.h"
/* Vanilla Geant4 includes */
#include "G4Box.hh"
#include "G4GDMLParser.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"


/* Handles for Geant4 and PUMAS */
static G4RunManager * manager = NULL;
static struct pumas_physics * physics = NULL;
static struct pumas_context * context = NULL;


/* Gracefully exit to the OS */
static int ExitGracefully(int rc)
{
        /* Clear wrapper data */
        g4pumas::ClearMedia();

        /* Clear Geant4 data */
        delete manager;

        /* Clear PUMAS data */
        if (context != nullptr)
                pumas_recorder_destroy(&context->recorder);
        pumas_context_destroy(&context);
        pumas_physics_destroy(&physics);

        exit(rc);
}


/* Error handler for PUMAS with a graceful exit */
static void HandleError(
    enum pumas_return, pumas_function_t, const char * message)
{
        /* Dump the error summary */
        std::fputs("pumas: library error. See details below\n", stderr);
        std::fprintf(stderr, "error: %s\n", message);

        /* Exit to the OS */
        ExitGracefully(EXIT_FAILURE);
}


/* Geant4 GDML geometry loader */
struct DetectorConstruction : public G4VUserDetectorConstruction
{
        G4VPhysicalVolume * Construct()
        {
                /* Load the GDML geometry */
                G4GDMLParser gdml;
                gdml.Read("examples/gdml/geometry.gdml");
                return gdml.GetWorldVolume();
        }
};


/* Empty physics. We need some in order to initialise Geant4 despite it
 * is actually not used in this example
 */
struct PhysicsList: public G4VModularPhysicsList {};


/* Printing callback to be called during PUMAS stepping */
static void PrintStep(struct pumas_context *,
    struct pumas_state * state, struct pumas_medium * medium,
    enum pumas_event event)
{
        if (event == PUMAS_EVENT_NONE) {
                return;
        } else if (event & PUMAS_EVENT_START) {
                /* Print a header */
                std::puts("");
                std::printf("  Energy          X          Y          Z     "
                    "Event  Medium\n");
                std::printf("   (GeV)         (m)        (m)        (m)\n\n");
        }

        std::string medium_name;
        if (medium != nullptr) {
                auto g4medium = (g4pumas::Geant4Medium *)medium;
                medium_name = g4medium->physical->GetLogicalVolume()->GetName();
        } else {
                medium_name = "(void)";
        }

        std::printf("%.3E   %10.3lf %10.3lf %10.3lf  %5d  %s\n", state->energy,
            state->position[0], state->position[1], state->position[2],
            event, medium_name.c_str());

        if (event & PUMAS_EVENT_STOP) {
                /* Print a blank tailer */
                std::puts("");
        }
}


int main ()
{
        /* Set the error handler callback. Whenever an error occurs during a
         * PUMAS function call, the supplied error handler will be evaluated,
         * resulting in an exit to the OS
         */
        pumas_error_handler_set(&HandleError);

        /* Initialise PUMAS physics from a binary dump, e.g. generated by the
         * `load` example
         */
        const char * dump_file = "materials/dump";
        auto fid = std::fopen(dump_file, "rb");
        if (fid == NULL) {
                std::perror(dump_file);
                ExitGracefully(EXIT_FAILURE);
        }
        pumas_physics_load(&physics, fid);
        std::fclose(fid);

        /* Initialize the G4 kernel */
        manager = new G4RunManager;
        manager->SetUserInitialization(new DetectorConstruction);
        manager->SetUserInitialization(new PhysicsList);
        manager->Initialize();

        /* Create a new PUMAS simulation context and configure it for backward
         * transport
         */
        pumas_context_create(&context, physics, 0);
        context->mode.direction = PUMAS_MODE_BACKWARD;

        /* Set the medium callback and the PRNG to the Geant4 wrappers */
        context->medium = &g4pumas::Medium;
        context->random = &g4pumas::Random;

        /* Attach a recorder for printing out Monte Carlo steps */
        pumas_recorder_create(&context->recorder, 0);
        context->recorder->record = &PrintStep;

        /* Run the Monte-Carlo */
        struct pumas_state state = {
            -1, 1, 0, 0, 0, 1, {0, 0, -1000.5}, {0, 0, -1}, 0};

        g4pumas::ResetTransport();
        pumas_context_transport(context, &state, NULL, NULL);

        ExitGracefully(EXIT_SUCCESS);
}
