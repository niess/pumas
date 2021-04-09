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

/* This example implements a smart materials loader for PUMAS. Whenever a
 * binary dump is found, PUMAS is initialised from it. Otherwise, PUMAS is
 * initialised from the MDF. In the later case, a binary dump is generated in
 * order to speed up further initialisations.
 *
 * This example also illustrates how to catch and raise errors issued by PUMAS.
 *
 * Note that for this example to work you need the corresponding MDF and energy
 * loss tables. Those can be downloaded with git, as following:
 *
 * git clone https://gitub.com/niess/pumas-materials materials
 */

/* Standard library includes */
#include <stdlib.h>
/* The PUMAS API */
#include "pumas.h"

/* Load PUMAS materials from a binary dump, if available, or from a MDF file
 * otherwise. On success, dump the shared data as a binary object
 */
static enum pumas_return load_pumas_materials(struct pumas_physics ** physics,
    const char * mdf, const char * dedx, const char * dump)
{
        const enum pumas_particle particle = PUMAS_PARTICLE_MUON;

        /* Check for a binary dump */
        pumas_error_catch(1);
        FILE * stream = fopen(dump, "rb");
        if (stream != NULL) {
                pumas_physics_load(physics, stream);
                fclose(stream);
                return pumas_error_raise();
        }
        pumas_error_catch(0);

        /* If no binary dump, initialise from the MDF and dump */
        enum pumas_return rc;
        if ((rc = pumas_physics_create(physics, particle, mdf, dedx, NULL)) !=
            PUMAS_RETURN_SUCCESS)
                return rc;

        /* Dump the library configuration */
        stream = fopen(dump, "wb+");
        if (stream == NULL) {
                /* Check for an error handler and call it whenever */
                rc = PUMAS_RETURN_PATH_ERROR;
                pumas_handler_cb * handler = pumas_error_handler_get();
                if (handler != NULL) handler(rc, NULL, NULL);
                return rc;
        }

        pumas_error_catch(1);
        pumas_physics_dump(*physics, stream);
        fclose(stream);
        return pumas_error_raise();
}

/* The executable main entry point */
int main(int argc, char * argv[])
{
        /* Parse any arguments */
        const char * mdf = (argc >= 2) ?
            argv[1] : "examples/data/materials.xml";
        const char * dedx = (argc >= 3) ? argv[2] : "@";
        const char * dump = (argc >= 4) ?
            argv[3] : "examples/data/materials.pumas";

        /* Load or compute the material's physics data */
        struct pumas_physics * physics;
        load_pumas_materials(&physics, mdf, dedx, dump);

        /* Clean and exit to the OS */
        pumas_physics_destroy(&physics);

        exit(EXIT_SUCCESS);
}
