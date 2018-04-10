/* This example implements a smart materials loader for PUMAS. Whenever a
 * binary dump is found, PUMAS is initialised from it. Otherwise, PUMAS is
 * initialised from the MDF. In the later case, a binary dump is generated in
 * order to speed up further initialisations.
 *
 * This example is in the public domain. Feel free to copy and modify it.
 */

/* Standard library includes */
#include <stdlib.h>
/* The PUMAS API */
#include "pumas.h"

/* Load PUMAS materials from a binary dump, if available, or from a MDF file
 * otherwise. On success, dump the shared data as a binary object
 */
static enum pumas_return load_pumas_materials(
    const char * mdf, const char * dedx, const char * dump)
{
        const enum pumas_particle particle = PUMAS_PARTICLE_MUON;

        /* Check for a binary dump */
        pumas_error_catch(1);
        FILE * stream = fopen(dump, "rb");
        if (stream != NULL) {
                pumas_load(stream);
                fclose(stream);
                return pumas_error_raise();
        }

        /* If no binary dump, initialise from the MDF and dump */
        enum pumas_return rc;
        if ((rc = pumas_initialise(particle, mdf, dedx, NULL)) !=
            PUMAS_RETURN_SUCCESS)
                return pumas_error_raise();

        /* Dump the library configuration */
        stream = fopen(dump, "wb+");
        if (stream == NULL) {
                /* Disable the error catching */
                pumas_error_catch(0);

                /* Check for an error handler and call it whenever */
                rc = PUMAS_RETURN_PATH_ERROR;
                pumas_handler_cb * handler = pumas_error_handler_get();
                if (handler != NULL) handler(rc, NULL, NULL);

                return rc;
        }
        pumas_dump(stream);
        fclose(stream);
        return pumas_error_raise();
}

/* The executable main entry point */
int main(int narg, char * argv[])
{
        /* Check the number of arguments */
        if (narg < 3) {
                fprintf(stderr,
                    "Usage: %s [path/to/mdf] [path/to/dedx] [path/to/dump]\n",
                    argv[0]);
                exit(EXIT_FAILURE);
        }

        /* Load and pre-compute the given material data */
        enum pumas_return rc = load_pumas_materials(argv[1], argv[2], argv[3]);

        /* Finalise PUMAS and return to the OS */
        pumas_finalise();
        exit((rc == PUMAS_RETURN_SUCCESS) ? EXIT_SUCCESS : EXIT_FAILURE);
}
