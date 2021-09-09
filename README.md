# PUMAS ![Build](https://github.com/niess/pumas/workflows/Build/badge.svg) [![codecov](https://codecov.io/gh/niess/pumas/branch/master/graph/badge.svg)](https://codecov.io/gh/niess/pumas) [![Documentation](https://readthedocs.org/projects/pumas/badge/?version=latest)](https://pumas.readthedocs.io/en/latest/?badge=latest)
( **S**emi **A**nalytical **MU**ons -or taus- **P**ropagation, *backwards* )

## Description

The PUMAS library allows to transport muon or tau leptons by forward or backward
Monte Carlo. The library is written in C99 with the standard library as sole
dependency. PUMAS is thread safe by design. The library also allows to generate
stopping power tables in the Particle Data Group
([PDG](https://pdg.lbl.gov/2020/AtomicNuclearProperties/index.html)) format. A
compilation of tables generated with PUMAS is available from the
[pumas-materials](https://github.com/niess/pumas-materials) repository.

## Installation

Building the library requires only the files [pumas.h](include/pumas.h) and
[pumas.c](src/pumas.c). On UNIX you might directly use the provided
[Makefile](Makefile). Alternatively, or on other platforms, CMake can be used
with the provided [CMakeLists.txt](CMakeLists.txt). More detailed instructions
are provided in the
[documentation](https://pumas.readthedocs.io/en/latest/installation).

## Documentation

The documentation is available online from [Read the
Docs](https://pumas.readthedocs.io/en/latest). It contains a description of the
library [API](https://pumas.readthedocs.io/en/latest/api) as well as
[tutorials](https://pumas.readthedocs.io/en/latest/tutorials). You might also
directly browse the [examples](examples) provided with the source.

## License
The PUMAS library is  under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files. The
[examples](examples) however have a separate public domain license allowing them
to be copied without any restriction.
