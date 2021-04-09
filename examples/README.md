# PUMAS examples
*Examples of usage of the PUMAS library.*

## Description

This folder contains several examples of usage of the PUMAS library organized
as follow:

-   The [pumas](pumas) folder contains examples requiring only the PUMAS
    library.

    -   [pumas/dump.c](pumas/dump.c) shows how to generate a PUMAS physics dump
        from a Material Decsription File (MDF) and energy loss tables. This dump
        can be used subquently for a fast initialisation of the physics.

    -   [pumas/straight.c](pumas/straight.c) shows how to compute a transmitted
        flux of muons through a constant thickness of a uniform material.

    -   [pumas/geometry.c](pumas/geometry.c) provides an example of geometry
        implementation composed of two layers: standard rock and air.

    -   [pumas/loader.c](pumas/straight.c) is an example of smart loader for
        physics data using a dump whenever available or a MDF.

-   The [turtle](turtle) folder contains an example of Earth geometry
    using the [TURTLE](https://github.com/niess/turtle) library.

-   The [geant4](geant4) folder contains examples of usage of PUMAS with
    [Geant4](https://geant4.web.cern.ch/node/1).

    - [geant4/g4pumas.h](geant4/g4pumas.h) and
      [geant4/g4pumas.cpp](geant4/g4pumas.cpp) are an example of Geant4 wrapper
      for PUMAS. It allows to use a `G4Navigator` from PUMAS. **Note** that this
      is **not** a thread safe implementation.

    - [geant4/run.cpp](geant4/run.cpp) illustrates how to use the `g4turtle`
      wrapper for propagating muons with PUMAS through a Geant4 geometry loaded
      from a GDML file ([data/geometry.gdml](data/geometry.gdml)).

    - [geant4/generate.cpp](geant4/generate.cpp) is an example of generation
       of a GDML geometry using Geant4.

## Installation

On UNIX the native PUMAS examples can be compiled be using the provided
[Makefile](../Makefile) e.g. as:
```bash
make examples
```
The compiled examples are located under the `bin/` folder e.g. as
`example-dump`.  Alternativelly CMake can be used as well by setting
`-DPUMAS_BUILD_EXAMPLES=ON`.

The other examples require a prior installation of
[Geant4](https://geant4.web.cern.ch/node/1) or of
[TURTLE](https://github.com/niess/turtle). Then the corresponding examples can
be compiled using the PUMAS Makefile as:
```bash
make examples-turtle
make examples-geant4
```
Building the [Geant4](https://geant4.web.cern.ch/node/1) or
[TURTLE](https://github.com/niess/turtle) examples with CMake is not currently
implemented.

## License
The examples are provided independently of the PUMAS library under a separate
public domain license allowing them to be copied without any restriction.
