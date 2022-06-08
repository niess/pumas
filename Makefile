# Installation prefix
PREFIX= $(PWD)


# PUMAS library
CFLAGS= -O3 -std=c99 -pedantic -Wall -Wfatal-errors -Iinclude
LIBS=   -lm

VERSION_MAJOR= $(shell grep PUMAS_VERSION_MAJOR include/pumas.h | cut -d' ' -f3)
VERSION_MINOR= $(shell grep PUMAS_VERSION_MINOR include/pumas.h | cut -d' ' -f3)
VERSION_PATCH= $(shell grep PUMAS_VERSION_PATCH include/pumas.h | cut -d' ' -f3)

LIB_NAME=      libpumas.so
LIB_SHORTNAME= $(LIB_NAME).$(VERSION_MAJOR)
LIB_FULLNAME=  $(LIB_SHORTNAME).$(VERSION_MINOR).$(VERSION_PATCH)

.PHONY: lib
lib: $(PREFIX)/lib/$(LIB_FULLNAME) \
     $(PREFIX)/lib/$(LIB_SHORTNAME) \
     $(PREFIX)/lib/$(LIB_NAME)

$(PREFIX)/lib/$(LIB_FULLNAME): src/pumas.c include/pumas.h | $(PREFIX)/lib
	@$(CC) -o $@ $(CFLAGS) -shared -fPIC $< $(LIBS)

$(PREFIX)/lib/$(LIB_SHORTNAME): $(PREFIX)/lib/$(LIB_FULLNAME)
	@ln -fs $(LIB_FULLNAME) $@

$(PREFIX)/lib/$(LIB_NAME): $(PREFIX)/lib/$(LIB_SHORTNAME)
	@ln -fs $(LIB_SHORTNAME) $@

$(PREFIX)/lib:
	@mkdir -p $@


# Cleanup
clean:
	@rm -rf $(PREFIX)/lib $(PREFIX)/bin


# PUMAS examples
LIBS_EXAMPLES= -L$(PREFIX)/lib -lpumas -Wl,-rpath,$(PREFIX)/lib $(LIBS)

.PHONY: examples
examples: $(PREFIX)/bin/example-dump \
          $(PREFIX)/bin/example-geometry \
          $(PREFIX)/bin/example-loader \
          $(PREFIX)/bin/example-geometry \
          $(PREFIX)/bin/example-straight \
          $(PREFIX)/bin/example-tabulate

$(PREFIX)/bin/example-geometry: examples/pumas/geometry.c \
                                examples/pumas/flux.c \
                                examples/pumas/flux.h \
                                $(PREFIX)/lib/libpumas.so | \
                                $(PREFIX)/bin
	@$(CC) $(CFLAGS) -o $@ examples/pumas/geometry.c \
	    examples/pumas/flux.c $(LIBS_EXAMPLES)

$(PREFIX)/bin/example-straight: examples/pumas/straight.c \
                                examples/pumas/flux.c \
                                examples/pumas/flux.h \
                                $(PREFIX)/lib/libpumas.so | \
                                $(PREFIX)/bin
	@$(CC) $(CFLAGS) -o $@ examples/pumas/straight.c \
	    examples/pumas/flux.c $(LIBS_EXAMPLES)

$(PREFIX)/bin/example-%: examples/pumas/%.c \
                         $(PREFIX)/lib/libpumas.so | \
                         $(PREFIX)/bin
	@$(CC) $(CFLAGS) -o $@ $< $(LIBS_EXAMPLES)

$(PREFIX)/bin:
	@mkdir -p $@


# TURTLE examples
PREFIX_TURTLE= $(PREFIX)
LIBS_TURTLE=   -L$(PREFIX_TURTLE)/lib -lturtle -Wl,-rpath,$(PREFIX_TURTLE)/lib

examples-turtle: $(PREFIX)/bin/example-turtle-earth

$(PREFIX)/bin/example-turtle-earth: examples/turtle/earth.c \
                                    $(PREFIX)/lib/libpumas.so \
                                    $(PREFIX_TURTLE)/lib/libturtle.so | \
                                    $(PREFIX)/bin
	@$(CC) $(CFLAGS) -o $@ $< $(LIBS_TURTLE) $(LIBS_EXAMPLES)


# Geant4 examples
CFLAGS_G4= -O2 -Wall $(shell geant4-config --cflags) -Iexamples/geant4
LIBS_G4=   $(shell geant4-config --libs) -lxerces-c

examples-geant4: $(PREFIX)/bin/example-geant4-generate \
                 $(PREFIX)/bin/example-geant4-run

$(PREFIX)/bin/example-geant4-generate: examples/geant4/generate.cpp | \
                                       $(PREFIX)/bin
	@$(CXX) $(CFLAGS_G4) -o $@ $< $(LIBS_G4)

$(PREFIX)/bin/example-geant4-run: examples/geant4/run.cpp \
                                  examples/geant4/g4pumas.cpp \
                                  examples/geant4/g4pumas.h \
                                  $(PREFIX)/lib/libpumas.so | \
                                  $(PREFIX)/bin
	@$(CXX) $(CFLAGS_G4) -Iinclude -o $@ examples/geant4/run.cpp \
	    examples/geant4/g4pumas.cpp $(LIBS_G4) $(LIBS_EXAMPLES)
