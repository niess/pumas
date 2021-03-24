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
 * `G4Navigator` from PUMAS. Note that this is not a thread safe implementation.
 */
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

        /* Random callback using the Geant4 PRNG */
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
