#pragma once

/* PUMAS API */
#include "pumas.h"

/* Forward declaration of Geant4 class(es) */
class G4VPhysicalVolume;


/* Geant4 wrappers for PUMAS */
namespace g4pumas {
        /* Extended medium structure linking back to Geant4 data */
        struct Geant4Medium {
                struct pumas_medium medium;
                G4VPhysicalVolume * physical;
        };

        /* Random callback using the Geant PRNG */
        double Random(struct pumas_context * context);

        /* Medium callback using a G4Navigator */
        enum pumas_step Medium(struct pumas_context * context,
            struct pumas_state * state, struct pumas_medium ** mediumPtr,
            double * stepPtr);

        /* Reset the navigator whenever a particle is relocated, e.g. for
         * a new transport
         */
        void ResetTransport();

        /* Clear the media mapping (release allocated memory) */
        void ClearMedia();
}
