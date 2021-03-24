# PUMAS ![Build](https://github.com/niess/pumas/workflows/Build/badge.svg) [![codecov](https://codecov.io/gh/niess/pumas/branch/master/graph/badge.svg)](https://codecov.io/gh/niess/pumas) [![Documentation](https://readthedocs.org/projects/pumas/badge/?version=latest)](https://pumas.readthedocs.io/en/latest/?badge=latest)
( **S**emi **A**nalytical **MU**ons -or taus- **P**ropagation, *backwards* )

## Description

PUMAS performs the transport of relativistic &mu; or &tau; in both forward
and backward Monte Carlo. The library is written in C99 with the Standard
Library as sole dependency. PUMAS is thread safe by design. The library also
allows to generate muons or taus energy loss tables in the
[PDG](http://pdg.lbl.gov/2016/AtomicNuclearProperties/index.html)
format. These tables, or those provided by the PDG, are needed as input for the
Physics.

## Installation

Building the library requires only the files [pumas.h](include/pumas.h) and
[pumas.c](src/pumas.c). On UNIX you might directly use the provided
[Makefile](Makefile). Alternatively, or on other platforms, CMake can be used
with the provided [CMakeLists.txt](CMakeLists.txt). More detailed instructions
can be found in the
[documentation](https://pumas.readthedocs.io/en/latest/installation).

## Materials

A compilation of materials for PUMAS is hosted on GitHub as a separate project,
[pumas-materials](https://github.com/niess/pumas-materials). These tables have
been generated with the PUMAS library. For muons one can also use the whole set
of tabulations provided online by the
[PDG](https://pdg.lbl.gov/2020/AtomicNuclearProperties/index.html). Note however
that those are less accurate than PUMAS ones above 100 TeV due to improvements
on the photonuclear cross-section (see e.g.
[Sokalski _et al._](https://arxiv.org/abs/hep-ph/0201122)).

## Documentation

The documentation is available online from [Read the
Docs](https://pumas.readthedocs.io/en/latest). It contains a description of the
library [API](https://pumas.readthedocs.io/en/latest/api) as well as
[tutorials](https://pumas.readthedocs.io/en/latest/tutorials). You might also
directly browse the [examples](examples) provided with the source.

_**Note** that for the examples to work you need their corresponding MDF and
energy loss tables distributed separately. Those can be downloaded with git, as
following:_
```bash
git clone https://gitub.com/niess/pumas-materials materials
```

## License
The PUMAS library is  under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files. The
[examples](examples) however have a separate public domain license allowing them
to be copied without any restriction.
