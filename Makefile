CFLAGS := -O0 -g -std=c99 -pedantic -Wall -Wfatal-errors
INCLUDES := -Iinclude
LIBS := -Llib -lpumas -Wl,-rpath,$(PWD)/lib -lm

.PHONY: lib clean examples

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
