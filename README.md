[![Build Status](https://travis-ci.com/niess/pumas.svg?branch=master)](https://travis-ci.com/niess/pumas)
[![codecov](https://codecov.io/gh/niess/pumas/branch/master/graph/badge.svg)](https://codecov.io/gh/niess/pumas)

# PUMAS
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
can be found on the [wiki](https://github.com/niess/pumas/wiki/Installation).

## Materials

A compilation of materials for PUMAS is hosted on GitHub as a separate project,
[pumas-materials](https://github.com/niess/pumas-materials). These tables have
been generated with the PUMAS library. For muons one can also use the whole set
of tabulations provided online by the
[PDG](http://pdg.lbl.gov/2019/AtomicNuclearProperties/index.html).  Note however
that those are slightly less accurate than PUMAS ones above 100 TeV.

## Documentation

The API documentation can be found
[here](https://niess.github.io/pumas-docs). For tutorials one can check the
[wiki](https://github.com/niess/pumas/wiki/Tutorials). You might also directly
browse the provided [examples](examples).

_**Note** that for the examples to work you need the corresponding MDF and
energy loss tables. Those can be downloaded with git, as following:_
```bash
git clone https://gitub.com/niess/pumas-materials materials
```

## License
The PUMAS library is  under the **GNU LGPLv3** license. See the provided
[LICENSE](LICENSE) and [COPYING.LESSER](COPYING.LESSER) files. The
[examples](examples) however have a separate public domain license allowing them
to be copied without any restriction.
