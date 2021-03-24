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

/* This program generates a GDML file for the examples. The GDML file
 * contains a set of materials consistent with the one used for PUMAS examples,
 * i.e. standard rock, water and air. A simple geometry is defined containing
 * a 1 km thick water layer on top of 1 km of rocks within a world bow filled
 * with air.
 *
 * Note that this example requires a prior installation of Geant4 with GDML
 * enabled (`-DGEANT4_USE_GDML=ON` CMake flag as well a `libxerces-c`).
 */

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"


int main ()
{
        /* Create the materials */
        auto H  = new G4Element(
             "H",  "H",  1.,  1.008700 * CLHEP::g / CLHEP::mole);
        auto C  = new G4Element(
             "C",  "C",  6., 12.010800 * CLHEP::g / CLHEP::mole);
        auto N  = new G4Element(
             "N",  "N",  7., 14.007200 * CLHEP::g / CLHEP::mole);
        auto O  = new G4Element(
             "O",  "O",  8., 15.999300 * CLHEP::g / CLHEP::mole);
        auto Rk = new G4Element(
            "Rk", "Rk", 11., 22.000000 * CLHEP::g / CLHEP::mole);
        auto Ar = new G4Element(
            "Ar", "Ar", 18., 39.948100 * CLHEP::g / CLHEP::mole);

        auto rock = new G4Material("StandardRock",
            2.65 * CLHEP::g / CLHEP::cm3, 1, kStateSolid);
        rock->AddElement(Rk, 1);
        rock->GetIonisation()->SetMeanExcitationEnergy(136.4 * CLHEP::eV);

        auto water = new G4Material("Water",
            1.00 * CLHEP::g / CLHEP::cm3, 2, kStateLiquid);
        water->AddElement(O, 0.888106);
        water->AddElement(H, 0.111894);
        water->GetIonisation()->SetMeanExcitationEnergy(79.7 * CLHEP::eV);

        auto air = new G4Material("Air",
            0.001205 * CLHEP::g / CLHEP::cm3, 4, kStateGas);
        air->AddElement( N, 0.755267);
        air->AddElement( O, 0.231781);
        air->AddElement(Ar, 0.012827);
        air->AddElement( C, 0.000124);
        air->GetIonisation()->SetMeanExcitationEnergy(85.7 * CLHEP::eV);

        /* Create the geometry */
        G4double h = 2 * CLHEP::km;
        auto worldGeometry = new G4Box("World", h, h, h);
        auto worldLogical = new G4LogicalVolume(
            worldGeometry, air, "World", 0, 0, 0);
        auto worldPhysical = new G4PVPlacement(
            0, G4ThreeVector(0., 0., 0.), worldLogical, "World", 0, false, 0);

        auto seabedGeometry = new G4Box("Seabed", h, h, 0.25 * h);
        auto seabedLogical = new G4LogicalVolume(
            seabedGeometry, rock, "Seabed", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 0., -0.75 * h), seabedLogical,
            "Seabed", worldLogical, false, 0);

        auto seaGeometry = new G4Box("Sea", h, h, 0.25 * h);
        auto seaLogical = new G4LogicalVolume(
            seaGeometry, water, "Sea", 0, 0, 0);
        new G4PVPlacement(0, G4ThreeVector(0., 0., -0.25 * h), seaLogical,
            "Sea", worldLogical, false, 0);

        G4GDMLParser gdml;
        gdml.Write("examples/data/geometry.gdml", worldPhysical);

        exit(EXIT_SUCCESS);
}
