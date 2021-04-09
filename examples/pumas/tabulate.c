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

/* This example generates PUMAS energy loss tabulations from a Material
 * Description File (MDF).
 */

/* Standard library includes */
#include <stdlib.h>
/* The PUMAS API */
#include "pumas.h"

/* The executable main entry point */
int main(int argc, char * argv[])
{
        /* Parse any arguments */
        char * mdf = (argc >= 2) ? argv[1] : "examples/data/materials.xml";
        char * dedx = (argc >= 3) ? argv[2] : "examples/data";

        /* Create a physics object in dry mode since only energy loss
         * tabulations are needed
         */
        struct pumas_physics * physics;
        struct pumas_physics_settings settings = { .dry = 1, .update = 1};
        pumas_physics_create(
            &physics, PUMAS_PARTICLE_MUON, mdf, dedx, &settings);

        /* Exit to the OS */
        exit(EXIT_SUCCESS);
}
