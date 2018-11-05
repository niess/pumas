CFLAGS := -O3 -std=c99 -pedantic -Wall -D_BUILD_TABULATE -Wfatal-errors
INCLUDES := -Iinclude
LIBS := -Llib -lpumas -Wl,-rpath,$(PWD)/lib -lm

.PHONY: lib clean examples tabulate

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

tabulate: bin/pumas-tabulate

bin/pumas-%: src/pumas-%.c lib
	@mkdir -p bin
	@$(CC) $(CFLAGS) -Wno-unused-function $(INCLUDES) -o $@ $< $(LIBS)
