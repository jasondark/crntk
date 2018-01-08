CC := gcc
CFLAGS := -fopenmp -std=c99 -pedantic -Wall -mtune=native -march=native -O2

obj/crntk.o: src/crntk.c include/crntk.h
	@mkdir -p obj
	$(CC) $(CFLAGS) -c src/crntk.c -o obj/crntk.o

example: src/example.c obj/crntk.o
	@mkdir -p bin
	$(CC) $(CFLAGS) src/example.c obj/crntk.o -lm -o bin/example

clean:
	rm -f obj/* && rm -f bin/*

