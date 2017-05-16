# PUMAS
( **S**emi **A**nalytical **MU**ons -or taus- **P**ropagation, *backwards* )

## Description

PUMAS performs the transport of relativistic μ's or &tau;'s' in both forward
and backward Monte-Carlo. The library is written in C99 with the Standard
Library as sole dependency. PUMAS is thread safe by design. The library is also
shipped with an executable, `pumas-tabulate`, allowing to generate muons or taus
energy loss tables in the
[PDG](http://pdglive.lbl.gov/2016/AtomicNuclearProperties/index.html)
format. These tables, or those provided by the PDG, are needed as input for the
PUMAS library initialisation.

## Installation

Building the library requires only the files [pumas.h](include/pumas.h) and
[pumas.c](src/pumas.c). In order to build the `pumas-tabulate` executable one
also needs [optparse](https://github.com/skeeto/optparse) (packaged with the
PUMAS source) and the files [pumas-tabulate.c](src/pumas-tabulate.c) and
[pumas-tabulate.h](src/pumas-tabulate.h). In the later case the PUMAS library
needs to be compiled with `_BUILD_TABULATE` being defined. On Linux you might
try and adapt the provided [Makefile](Makefile). Alternatively, CMake can also
be used on other systems with the provided [CMakeLists.txt](CMakeLists.txt).

## Documentation

The API documentation can be found
[here](https://niess.github.io/pumas/docs/index.html#HEAD).

## Materials

A compilation of materials for PUMAS is hosted on GitHub as separate project,
[pumas-materials](https://github.com/niess/pumas-materials). These tables have
been generated with `pumas-tabulate`. For muons, one can also use the whole set
of tabulations provided online by the
[PDG](http://pdg.lbl.gov/2016/AtomicNuclearProperties/index.html).
Note however that those are slightly less accurate than PUMAS ones, above
100 TeV.

## License
The PUMAS library is  under the **GNU LGPLv3** license. See the provided
`LICENSE` and `COPYING.LESSER` files.
