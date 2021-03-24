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

/* This is an example of Geant4 wrapper for PUMAS. It allows to use a
 * `G4Navigator` from PUMAS. Note that for the sake of clarity this is not a
 * thread safe implementation.  It could be adapted for multithreaded usage by
 * migrating the static globals data (e.g. volumes2media) as `user_data` of a
 * `pumas_context`.
 */

/* Geant4 includes */
#include "G4Material.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "Randomize.hh"
/* Geant4 wrapper for PUMAS */
#include "g4pumas.h"


/* Mapping between Geant4 volumes and PUMAS media */
static std::map<G4VPhysicalVolume *, struct pumas_medium *> volumes2media;


void g4pumas::ClearMedia()
{
        for (auto& kv : volumes2media) {
                delete kv.second;
        }
        volumes2media.clear();
}

/* Wrapper for the Geant4 random engine */
double g4pumas::Random(struct pumas_context *)
{
        return G4UniformRand();
}


/* Wrapper for local properties, i.e. the medium density  */
static double Locals(struct pumas_medium * medium,
    struct pumas_state *, struct pumas_locals * locals)
{
        auto g4Medium = (g4pumas::Geant4Medium *)medium;
        if (g4Medium->physical != nullptr) {
                auto material =
                    g4Medium->physical->GetLogicalVolume()->GetMaterial();
                locals->density = material->GetDensity() /
                    CLHEP::kg * CLHEP::m3;
        }

        return 0;
}


/* Wrapper for the the Geant4 geometry */
enum pumas_step g4pumas::Medium(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** mediumPtr,
    double * stepPtr)
{
        G4ThreeVector r(state->position[0] * CLHEP::m,
            state->position[1] * CLHEP::m, state->position[2] * CLHEP::m);
        const double sgn =
            (context->mode.direction == PUMAS_MODE_FORWARD) ? 1 : -1;
        G4ThreeVector u(sgn * state->direction[0], sgn * state->direction[1],
            state->direction[2]);

        auto navigator = G4TransportationManager::GetTransportationManager()
            ->GetNavigatorForTracking();
        auto physical = navigator->LocateGlobalPointAndSetup(r, &u);

        if (physical == nullptr) {
                /* Whenever this could happen let us return a void
                 * medium
                 */
                if (mediumPtr != nullptr) *mediumPtr = nullptr;
                if (stepPtr != nullptr) *stepPtr = 0;
                return PUMAS_STEP_CHECK;
        }

        /* Fetch the physical volume */
        if (mediumPtr != nullptr) {
                /* Fetch or register the volume if new */
                auto medium = volumes2media[physical];
                if (medium == nullptr) {
                        /* Get the PUMAS material index */
                        int index;
                        auto material = physical->GetLogicalVolume()
                            ->GetMaterial();
                        auto physics = pumas_context_physics_get(context);
                        auto rc = pumas_physics_material_index(physics,
                            material->GetName().c_str(), &index);
                        if (rc != PUMAS_RETURN_SUCCESS) {
                                /* If errors are silenced then return a void
                                 * medium.
                                 */
                                *mediumPtr = nullptr;
                                if (stepPtr != nullptr) *stepPtr = 0;
                                return PUMAS_STEP_CHECK;
                        }

                        /* Allocate and register the new medium */
                        auto tmp = new g4pumas::Geant4Medium;
                        tmp->medium.material = index;
                        tmp->medium.locals = &Locals;
                        tmp->physical = physical;
                        medium = (struct pumas_medium *)tmp;
                        volumes2media[physical] = medium;
                }
                *mediumPtr = medium;
        }

        /* Compute the step length */
        if (stepPtr != nullptr) {
                G4double safety = 0.;
                *stepPtr = navigator->ComputeStep(
                    r, u, kInfinity, safety) / CLHEP::m;
        }

        /* Note that Geant4 treatment of boundaries is not consistent with
         * PUMAS algorithm. Therefore we cannot use PUMAS_STEP_RAW.
         */
        return PUMAS_STEP_CHECK;
}


void g4pumas::ResetTransport()
{
        G4TransportationManager::GetTransportationManager()
            ->GetNavigatorForTracking()->ResetStackAndState();
}
