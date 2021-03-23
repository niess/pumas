CFLAGS := -O2 -std=c99 -pedantic -Wall -Wfatal-errors
INCLUDES := -Iinclude
LIBS := -Llib -lpumas -Wl,-rpath,$(PWD)/lib -lm

.PHONY: lib clean examples examples-turtle examples-geant4


lib: lib/libpumas.so
	@rm -f *.o

lib/libpumas.so: src/pumas.c include/pumas.h
	@mkdir -p lib
	@$(CC) -o $@ $(CFLAGS) -fPIC -shared $(INCLUDES) $<

clean:
	@rm -rf lib bin


examples: bin/example-straight bin/example-load bin/example-geometry

bin/example-%: examples/%.c lib/libpumas.so
	@mkdir -p bin
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)


examples-turtle: bin/example-earth

bin/example-earth: examples/earth.c lib/libpumas.so lib/libturtle.so
	@mkdir -p bin
	@$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS) -lturtle


G4FLAGS = -O2 -Wall $(shell geant4-config --cflags) -Iexamples/geant4
G4LIBS = $(shell geant4-config --libs) -lxerces-c

examples-geant4: bin/example-geant4-generate bin/example-geant4-run

bin/example-geant4-generate: examples/geant4/generate.cpp
	@mkdir -p bin
	@$(CXX) $(G4FLAGS) -o $@ examples/geant4/generate.cpp $(G4LIBS)

bin/example-geant4-run: examples/geant4/run.cpp examples/geant4/g4pumas.cpp examples/geant4/g4pumas.h
	@mkdir -p bin
	@$(CXX) $(G4FLAGS) $(INCLUDES) -o $@ examples/geant4/run.cpp \
	    examples/geant4/g4pumas.cpp $(G4LIBS) $(LIBS)
