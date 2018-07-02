CFLAGS := -O2 -std=c99 -pedantic -Wall -D_BUILD_TABULATE
INC := -Iinclude
LIBS := -Llib -lpumas -Wl,-rpath,$(PWD)/lib -lm

.PHONY: lib clean examples tabulate

lib: lib/libpumas.so
	@rm -f *.o

lib/libpumas.so: src/pumas.c include/pumas.h
	@mkdir -p lib
	@$(CC) -o $@ $(CFLAGS) -fPIC -shared $(INC) $<

clean:
	@rm -rf lib bin

examples: bin/example-straight bin/example-load bin/example-geometry

bin/example-%: examples/%.c lib/libpumas.so
	@mkdir -p bin
	@$(CC) $(CFLAGS) $(INC) -o $@ $< $(LIBS)

tabulate: bin/pumas-tabulate

bin/pumas-%: src/pumas-%.c lib
	@mkdir -p bin
	@$(CC) $(CFLAGS) -Wno-unused-function $(INC) -o $@ $< $(LIBS)
