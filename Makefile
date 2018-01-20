CC := gcc -fopenmp
CFLAGS := -std=c99 -Wall -Wextra -Wpedantic -Wno-unused-parameter -mtune=native -march=native -O2

.PHONY: dirs examples clean

dirs:
	@mkdir -p bin obj

obj/crntk.o: dirs src/crntk.c include/crntk.h
	$(CC) $(CFLAGS) -c src/crntk.c -o obj/crntk.o

bin/nullspace: examples/nullspace.c obj/crntk.o
	$(CC) $(CFLAGS) examples/nullspace.c obj/crntk.o -lm -o bin/nullspace

bin/extinction: examples/extinction.c obj/crntk.o
	$(CC) $(CFLAGS) examples/extinction.c obj/crntk.o -lm -o bin/extinction

bin/birth: examples/extinction.c obj/crntk.o
	$(CC) $(CFLAGS) examples/birth.c obj/crntk.o -lm -o bin/birth

examples: bin/nullspace bin/extinction bin/birth

clean:
	rm -rf obj && rm -rf bin

