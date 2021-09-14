# PUMAS examples
*Examples of usage of the PUMAS library.*

## Description

This folder contains several examples of usage of the PUMAS library organized
as follow:

-   The [pumas](pumas) folder contains examples requiring only the PUMAS
    library.

    -   [pumas/tabulate.c](pumas/tabulate.c) shows how to generate stopping
        power tables in the PDG format using PUMAS and a Materials Description
        File (MDF).

    -   [pumas/dump.c](pumas/dump.c) shows how to generate a binary dump of
        PUMAS physics tables from a MDF. This dump can be used subquently for a
        fast initialisation of the physics, e.g. as in the
        [pumas/geometry.c](pumas/geometry.c) example below.

    -   [pumas/straight.c](pumas/straight.c) shows how to compute a transmitted
        flux of muons through a constant thickness of a uniform material.

    -   [pumas/geometry.c](pumas/geometry.c) provides an example of geometry
        implementation composed of two layers: standard rock and air.

    -   [pumas/loader.c](pumas/straight.c) is an example of smart loader for
        physics data using a binary dump when available or a MDF otherwise.

-   The [turtle](turtle) folder contains an example of Earth geometry
    using the [TURTLE][TURTLE] library.

-   The [geant4](geant4) folder contains examples of usage of PUMAS with
    [Geant4][Geant4].

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

On UNIX the examples under [pumas/](pumas/) can be compiled be using the
provided [Makefile](../Makefile) e.g. as:
```bash
make examples
```
The compiled examples are located under the `bin/` folder e.g. as
`example-dump`.  Alternativelly CMake can be used as well by setting
`-DPUMAS_BUILD_EXAMPLES=ON`.

The other examples require a prior installation of [Geant4][Geant4] or of
[TURTLE][TURTLE]. Then the corresponding examples can be compiled using the
[Makefile](../Makefile) as:
```bash
make examples-turtle
make examples-geant4
```
Building the [Geant4][Geant4] or [TURTLE][TURTLE] examples with CMake is not
currently implemented.

## License
The examples are provided independently of the PUMAS library under a separate
public domain license allowing them to be copied without any restriction.


[Geant4]: https://geant4.web.cern.ch
[TURTLE]: https://niess.github.io/turtle-pages
