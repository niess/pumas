CC := gcc
CFLAGS := -O2 -std=c99 -pedantic -fPIC -Wall -D_BUILD_TABULATE
INC := -Iinclude

.PHONY: lib clean tabulate

lib: lib/libpumas.so
	@rm -f *.o

lib/lib%.so: %.o
	@mkdir -p lib
	@$(CC) -o $@ $(CFLAGS) -shared $(INC) $<

%.o: src/%.c include/%.h
	@$(CC) $(CFLAGS) $(INC) -o $@ -c $<

clean:
	@rm -rf lib bin *.o

tabulate: bin/pumas-tabulate
	@rm -f pumas.o

bin/pumas-%: src/pumas-%.c lib
	@mkdir -p bin
	@$(CC) $(CFLAGS) $(INC) -o $@ $< -lm -Llib -lpumas
