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

/* This example initialise PUMAS physics from a Material Description File (MDF)
 * and generates a binary dump for subsequent usage.
 *
 * Note that this example needs to be run before some other examples since the
 * generated physics dump is needed by those. Note also that for this example
 * to work you need the corresponding MDF and energy loss tables. Those can be
 * downloaded with git, as following:
 *
 * git clone https://gitub.com/niess/pumas-materials materials
 */

/* Standard library includes */
#include <stdlib.h>
/* The PUMAS API */
#include "pumas.h"

/* The executable main entry point */
int main(int argc, char * argv[])
{
        /* Parse any arguments */
        const char * mdf = (argc >= 2) ?
            argv[1] : "materials/mdf/examples/standard.xml";
        const char * dedx = (argc >= 3) ? argv[2] : "materials/dedx/muon";
        const char * dump = (argc >= 4) ? argv[3] : "materials/examples.pumas";

        /* Compute the material's physics data */
        struct pumas_physics * physics;
        pumas_physics_create(&physics, PUMAS_PARTICLE_MUON, mdf, dedx);

        /* Dump the physics data for subsquent usage */
        FILE * stream = fopen(dump, "wb");
        if (stream != NULL) {
                pumas_physics_dump(physics, stream);
                fclose(stream);
        } else {
                fprintf(stderr, "error: could not open %s\n", mdf);
                pumas_physics_destroy(&physics);

                exit(EXIT_FAILURE);
        }

        /* Clean and exit to the OS */
        pumas_physics_destroy(&physics);

        exit(EXIT_SUCCESS);
}
