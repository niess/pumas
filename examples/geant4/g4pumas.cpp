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
                return PUMAS_STEP_APPROXIMATE;
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
                                return PUMAS_STEP_APPROXIMATE;
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

        return PUMAS_STEP_APPROXIMATE;
}


void g4pumas::ResetTransport()
{
        G4TransportationManager::GetTransportationManager()
            ->GetNavigatorForTracking()->ResetStackAndState();
}
